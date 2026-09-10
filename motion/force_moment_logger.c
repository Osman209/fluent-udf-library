/*
 * force_moment_logger.c
 *
 * Writes the pressure force, the viscous force and the moment on a wall
 * zone to a CSV file at the end of every time step. Parallel safe: only
 * principal faces are summed (no double counting across partition
 * boundaries), partial sums are reduced with PRF_GRSUM1, and one
 * process writes the file.
 *
 * Hook: User-Defined > Function Hooks > Execute at End.
 * Set WALL_ZONE_ID to the zone ID of the body wall (shown in the
 * Boundary Conditions panel). MOMENT_CENTER is the point moments are
 * taken about, in the global frame.
 *
 * Pressure force
 *   F_AREA on a boundary face points out of the fluid, into the wall,
 *   so the pressure force of the fluid on the wall is p * A with
 *   p = F_P(f,t), the gauge pressure. On a closed body the operating
 *   pressure cancels over the surface. On an open wall (a plate wetted
 *   on one side only) it does not, so add it yourself if that matters.
 *
 * Viscous force: two methods, chosen by VISC_METHOD
 *   1 = velocity gradient (default). The viscous stress at the wall is
 *       built from the effective viscosity and the velocity gradient in
 *       the cell next to the wall, using documented macros (C_MU_EFF,
 *       C_U_G, C_V_G, C_W_G). The sign follows from the formula, and
 *       the code works in any Fluent version that has those macros.
 *       With wall functions and a coarse near-wall mesh this
 *       underestimates the shear, exactly as any cell-gradient estimate
 *       does; on a y+ ~ 1 mesh it is accurate.
 *   2 = SV_WALL_SHEAR. Reads Fluent's own stored wall shear vector, so
 *       it matches Report > Forces on any mesh. This macro is NOT in
 *       the Ansys UDF manual: its sign convention is not published and
 *       has been reported to differ from the contour plot, and it has
 *       been reported to crash in some versions. If you use it, run
 *       calibrate_viscous_sign below once and compare with
 *       Report > Forces, then set VISC_SIGN.
 *
 * Method 1 is the default because it is documented and its sign is
 * derivable. Use method 2 when you need to match the Fluent report
 * exactly on a wall-function mesh.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define WALL_ZONE_ID   7
#define MOMENT_CENTER  {0.0, 0.0, 0.0}
#define VISC_METHOD    1      /* 1 = velocity gradient, 2 = SV_WALL_SHEAR */
#define VISC_SIGN     -1.0    /* only used by method 2                    */
#define LOG_FILE      "forces.csv"
/* --------------------------------------------------------------------- */

static int first = 1;

/* Viscous force on one boundary face, from the velocity gradient in the
 * adjacent cell. tau = mu_eff (grad u + grad u^T); the traction on the
 * face is tau . n_hat; the force on the wall is traction * area, with
 * the normal component removed so only shear is left. */
#if !RP_HOST
static void visc_force_gradient(face_t f, Thread *t, real fv[3])
{
    cell_t c0 = F_C0(f, t);
    Thread *t0 = THREAD_T0(t);
    real A[ND_ND], nhat[3] = {0.0, 0.0, 0.0}, area, mu;
    real g[3][3], tau[3][3], trac[3] = {0.0, 0.0, 0.0}, tn = 0.0;
    int i, j;

    F_AREA(A, f, t);
    area = NV_MAG(A);
    if (area <= 0.0) { fv[0] = fv[1] = fv[2] = 0.0; return; }
    for (i = 0; i < ND_ND; i++) nhat[i] = A[i] / area;

    mu = C_MU_EFF(c0, t0);

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) g[i][j] = 0.0;
    for (j = 0; j < ND_ND; j++)
    {
        g[0][j] = C_U_G(c0, t0)[j];
        g[1][j] = C_V_G(c0, t0)[j];
#if ND_ND == 3
        g[2][j] = C_W_G(c0, t0)[j];
#endif
    }

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            tau[i][j] = mu * (g[i][j] + g[j][i]);

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            trac[i] += tau[i][j] * nhat[j];

    for (i = 0; i < 3; i++) tn += trac[i] * nhat[i];
    for (i = 0; i < 3; i++) fv[i] = (trac[i] - tn * nhat[i]) * area;
}
#endif

DEFINE_EXECUTE_AT_END(log_forces)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t = Lookup_Thread(d, WALL_ZONE_ID);
    face_t f;
    real A[ND_ND], xc[ND_ND], r[3];
    real Fp[3] = {0.0, 0.0, 0.0};
    real Fv[3] = {0.0, 0.0, 0.0};
    real M[3]  = {0.0, 0.0, 0.0};
    real xm[3] = MOMENT_CENTER;
    real p, fp3[3], fv3[3], ftot[3];
    int i;

    if (t == NULL)
    {
        if (UDF_IS_WRITER)
            Message("log_forces: zone %d not found, nothing written\n", WALL_ZONE_ID);
        return;
    }

    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            F_AREA(A, f, t);
            F_CENTROID(xc, f, t);
            p = F_P(f, t);

            for (i = 0; i < 3; i++) { fp3[i] = 0.0; fv3[i] = 0.0; r[i] = 0.0; }
            for (i = 0; i < ND_ND; i++)
            {
                fp3[i] = p * A[i];
                r[i]   = xc[i] - xm[i];
            }
#if VISC_METHOD == 1
            visc_force_gradient(f, t, fv3);
#else
            for (i = 0; i < ND_ND; i++)
                fv3[i] = VISC_SIGN * F_STORAGE_R_N3V(f, t, SV_WALL_SHEAR)[i];
#endif
            for (i = 0; i < 3; i++)
            {
                ftot[i] = fp3[i] + fv3[i];
                Fp[i]  += fp3[i];
                Fv[i]  += fv3[i];
            }
            M[0] += r[1] * ftot[2] - r[2] * ftot[1];
            M[1] += r[2] * ftot[0] - r[0] * ftot[2];
            M[2] += r[0] * ftot[1] - r[1] * ftot[0];
        }
    }
    end_f_loop(f, t)

    for (i = 0; i < 3; i++)
    {
        Fp[i] = PRF_GRSUM1(Fp[i]);
        Fv[i] = PRF_GRSUM1(Fv[i]);
        M[i]  = PRF_GRSUM1(M[i]);
    }

    if (UDF_IS_WRITER)
    {
        FILE *fpf;
        if (first)
        {
            fpf = fopen(LOG_FILE, "w");
            if (fpf)
            {
                fprintf(fpf, "t,Fpx,Fpy,Fpz,Fvx,Fvy,Fvz,Mx,My,Mz\n");
                fclose(fpf);
            }
            first = 0;
        }
        fpf = fopen(LOG_FILE, "a");
        if (fpf)
        {
            fprintf(fpf, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                    CURRENT_TIME, Fp[0], Fp[1], Fp[2],
                    Fv[0], Fv[1], Fv[2], M[0], M[1], M[2]);
            fclose(fpf);
        }
    }
#endif
}

/*
 * Run once from User-Defined > Execute on Demand after a converged
 * iteration. Prints the pressure force and the viscous force computed
 * both ways, so you can compare with Report > Forces on the same zone
 * and pick VISC_METHOD (and VISC_SIGN, if you choose method 2).
 */
DEFINE_ON_DEMAND(calibrate_viscous_sign)
{
#if !RP_HOST
    Domain *d = Get_Domain(1);
    Thread *t = Lookup_Thread(d, WALL_ZONE_ID);
    face_t f;
    real A[ND_ND], p, fv3[3];
    real Fp[3] = {0.0, 0.0, 0.0};
    real Fg[3] = {0.0, 0.0, 0.0};
    real Fs[3] = {0.0, 0.0, 0.0};
    int i;

    if (t == NULL) { if (UDF_IS_WRITER) Message("zone %d not found\n", WALL_ZONE_ID); return; }

    begin_f_loop(f, t)
    {
        if (PRINCIPAL_FACE_P(f, t))
        {
            F_AREA(A, f, t);
            p = F_P(f, t);
            for (i = 0; i < ND_ND; i++) Fp[i] += p * A[i];

            visc_force_gradient(f, t, fv3);
            for (i = 0; i < 3; i++) Fg[i] += fv3[i];

            for (i = 0; i < ND_ND; i++)
                Fs[i] += F_STORAGE_R_N3V(f, t, SV_WALL_SHEAR)[i];
        }
    }
    end_f_loop(f, t)

    for (i = 0; i < 3; i++)
    {
        Fp[i] = PRF_GRSUM1(Fp[i]);
        Fg[i] = PRF_GRSUM1(Fg[i]);
        Fs[i] = PRF_GRSUM1(Fs[i]);
    }

    if (UDF_IS_WRITER)
    {
        Message("\n--- force calibration on zone %d ---\n", WALL_ZONE_ID);
        Message("pressure force            : %g %g %g\n", Fp[0], Fp[1], Fp[2]);
        Message("viscous, gradient method  : %g %g %g\n", Fg[0], Fg[1], Fg[2]);
        Message("viscous, SV_WALL_SHEAR raw: %g %g %g\n", Fs[0], Fs[1], Fs[2]);
        Message("Compare with Report > Forces on this zone.\n");
        Message("If the raw SV_WALL_SHEAR row has the opposite sign to the\n");
        Message("report, set VISC_SIGN to -1; if it matches, set it to +1.\n\n");
    }
#endif
}
