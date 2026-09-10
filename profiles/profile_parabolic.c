/*
 * profile_parabolic.c
 *
 * Fully developed laminar inlet profiles, for pipes and channels where
 * you do not want to waste an entry length on the mesh.
 *
 *   parabolic_pipe    round pipe, radius R about axis (AXIS_X0, AXIS_Y0)
 *                     U(r) = 2 U_mean (1 - (r/R)^2)
 *   parabolic_channel plane channel of half-height H
 *                     U(y) = 1.5 U_mean (1 - (y/H)^2)
 *
 * Both are normalised so that the area-average velocity is U_MEAN
 * exactly, which means the mass flow through the inlet is
 * rho * U_MEAN * area, whatever the mesh. That is what you check.
 *
 * Hook on the velocity inlet, the component normal to the inlet.
 * FLOW_DIR selects which velocity component the profile is written to
 * and which coordinates form the cross-section.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define U_MEAN     1.0    /* m/s, area-averaged velocity                */
#define PIPE_R     0.05   /* m, pipe radius (parabolic_pipe)            */
#define CHAN_H     0.05   /* m, channel half-height (parabolic_channel) */
#define AXIS_A0    0.0    /* m, pipe axis position, first cross coord   */
#define AXIS_B0    0.0    /* m, pipe axis position, second cross coord  */
#define WALL_MID   0.0    /* m, channel mid-plane position              */
#define FLOW_DIR   0      /* 0 = x, 1 = y, 2 = z                        */
#define CROSS_A    1      /* first cross-section coordinate             */
#define CROSS_B    2      /* second cross-section coordinate (3D only)  */
#define CHAN_DIR   1      /* coordinate across the channel              */
/* --------------------------------------------------------------------- */

DEFINE_PROFILE(parabolic_pipe, t, i)
{
    face_t f;
    real x[ND_ND], da, db, r2, R2 = PIPE_R * PIPE_R, u;

    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        da = x[CROSS_A] - AXIS_A0;
#if ND_ND == 3
        db = x[CROSS_B] - AXIS_B0;
#else
        db = 0.0;
#endif
        r2 = da * da + db * db;
        u  = 2.0 * U_MEAN * (1.0 - r2 / R2);
        F_PROFILE(f, t, i) = (u > 0.0) ? u : 0.0;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(parabolic_channel, t, i)
{
    face_t f;
    real x[ND_ND], y, u;

    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        y = (x[CHAN_DIR] - WALL_MID) / CHAN_H;
        u = 1.5 * U_MEAN * (1.0 - y * y);
        F_PROFILE(f, t, i) = (u > 0.0) ? u : 0.0;
    }
    end_f_loop(f, t)
}

/*
 * Turbulent pipe inlet, 1/7 power law:
 *   U(r) = U_max (1 - r/R)^(1/n),  U_mean/U_max = 2 n^2 / ((n+1)(2n+1))
 * with n = 7 this gives U_mean = 0.8167 U_max. Normalised to U_MEAN.
 */
#define POW_N 7.0

DEFINE_PROFILE(power_law_pipe, t, i)
{
    face_t f;
    real x[ND_ND], da, db, r, u, umax;

    umax = U_MEAN * ((POW_N + 1.0) * (2.0 * POW_N + 1.0)) / (2.0 * POW_N * POW_N);

    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        da = x[CROSS_A] - AXIS_A0;
#if ND_ND == 3
        db = x[CROSS_B] - AXIS_B0;
#else
        db = 0.0;
#endif
        r = sqrt(da * da + db * db);
        if (r > PIPE_R) r = PIPE_R;
        u = umax * pow(1.0 - r / PIPE_R, 1.0 / POW_N);
        F_PROFILE(f, t, i) = (u > 0.0) ? u : 0.0;
    }
    end_f_loop(f, t)
}
