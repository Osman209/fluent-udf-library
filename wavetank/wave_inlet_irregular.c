/*
 * wave_inlet_irregular.c
 *
 * Irregular (random sea) wave generation at a velocity inlet: linear
 * superposition of N_COMP components drawn from a JONSWAP spectrum
 * (significant height HS, peak period TP, peakedness GAMMA).
 *
 * Profiles: inlet_u_irr, inlet_v_irr, inlet_vof_irr, hooked exactly as
 * in wave_inlet_regular.c.
 *
 * The component phases come from a small deterministic generator with
 * a fixed SEED, so the same input always gives the same sea state and
 * the Python script validation/jonswap_check.py reproduces the exact
 * component list (amplitudes, wavenumbers, phases) written by this UDF
 * to wave_components.csv on the first call. That file is the record of
 * which sea state the run used.
 *
 * Checks printed at start-up: Hs recovered from the discretised
 * spectrum as 4 sqrt(m0); it should be within a few percent of HS.
 * If it is far off, widen [F_MIN, F_MAX] or increase N_COMP.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"
#include "wave_theory.h"

/* ---- user parameters ------------------------------------------------ */
#define HS        0.10    /* m, significant wave height                  */
#define TP        1.50    /* s, peak period                              */
#define GAMMA     3.3     /* JONSWAP peakedness (3.3 = North Sea mean)   */
#define N_COMP    100     /* number of components                        */
#define F_MIN     (0.4/TP)/* Hz, lower cut of the spectrum               */
#define F_MAX     (3.0/TP)/* Hz, upper cut                               */
#define SEED      12345u
#define DEPTH_D   1.00
#define Y_SWL     0.00
#define GRAV      9.81
#define RAMP_S    3.0     /* s, amplitude ramp                           */
#define VOF_SMOOTH 0.01
#define COMP_FILE "wave_components.csv"
/* --------------------------------------------------------------------- */

static int    init_done = 0;
static double ca[N_COMP], ck[N_COMP], cw[N_COMP], cp[N_COMP];

static void irr_init(void)
{
    double m0;
    int i;

    if (init_done) return;
    m0 = wt_build_components(N_COMP, HS, TP, GAMMA, F_MIN, F_MAX, DEPTH_D,
                             GRAV, SEED, ca, ck, cw, cp);
    init_done = 1;

    if (UDF_IS_WRITER)
    {
        FILE *fp;
        Message("wave_inlet_irregular: Hs=%g Tp=%g gamma=%g N=%d -> "
                "4 sqrt(m0) = %g\n", HS, TP, GAMMA, N_COMP, 4.0 * sqrt(m0));
        fp = fopen(COMP_FILE, "w");
        if (fp)
        {
            fprintf(fp, "# Hs=%g Tp=%g gamma=%g N=%d fmin=%g fmax=%g d=%g seed=%u\n",
                    HS, TP, GAMMA, N_COMP, (double)F_MIN, (double)F_MAX, DEPTH_D, SEED);
            fprintf(fp, "i,amplitude,wavenumber,omega,phase\n");
            for (i = 0; i < N_COMP; i++)
                fprintf(fp, "%d,%.10g,%.10g,%.10g,%.10g\n", i, ca[i], ck[i], cw[i], cp[i]);
            fclose(fp);
        }
    }
}

static double ramp(double t)
{
    if (RAMP_S <= 0.0 || t >= RAMP_S) return 1.0;
    if (t <= 0.0) return 0.0;
    return 0.5 * (1.0 - cos(M_PI * t / RAMP_S));
}

static double irr_at(double x, double y, double t, double *u, double *v)
{
    double yy = y - Y_SWL, rr = ramp(t), eta;
    eta = wt_irregular(N_COMP, ca, ck, cw, cp, DEPTH_D, x, yy, t, u, v);
    eta *= rr; *u *= rr; *v *= rr;
    if (yy > eta) { *u = 0.0; *v = 0.0; }
    return eta;
}

DEFINE_PROFILE(inlet_u_irr, t, i)
{
    face_t f; real x[ND_ND]; double u, v, tt = CURRENT_TIME;
    irr_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        irr_at(x[0], x[1], tt, &u, &v);
        F_PROFILE(f, t, i) = u;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_v_irr, t, i)
{
    face_t f; real x[ND_ND]; double u, v, tt = CURRENT_TIME;
    irr_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        irr_at(x[0], x[1], tt, &u, &v);
        F_PROFILE(f, t, i) = v;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_vof_irr, t, i)
{
    face_t f; real x[ND_ND]; double u, v, eta, a, tt = CURRENT_TIME;
    irr_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        eta = irr_at(x[0], x[1], tt, &u, &v);
        a = 0.5 + (eta - (x[1] - Y_SWL)) / (VOF_SMOOTH > 0.0 ? VOF_SMOOTH : 1.0e-9);
        if (a > 1.0) a = 1.0;
        if (a < 0.0) a = 0.0;
        F_PROFILE(f, t, i) = a;
    }
    end_f_loop(f, t)
}
