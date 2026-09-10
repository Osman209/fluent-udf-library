/*
 * wave_inlet_regular.c
 *
 * Regular wave generation at a velocity inlet for a VOF numerical wave
 * tank (inlet-velocity method). Provides three profiles:
 *
 *   inlet_u     horizontal velocity   -> velocity inlet, X-Velocity
 *   inlet_v     vertical velocity     -> velocity inlet, Y-Velocity
 *   inlet_vof   water volume fraction -> velocity inlet, phase "water",
 *                                        Volume Fraction
 *
 * Theory: WAVE_THEORY 1 = Airy (linear), 2 = Stokes second order.
 * The maths lives in common/wave_theory.h and is checked outside Fluent
 * against validation/wave_theory_check.py.
 *
 * Coordinates: x along the tank, y vertical, still water level at
 * y = Y_SWL, depth D (bed at Y_SWL - D). Air above the free surface gets
 * zero velocity and zero water fraction.
 *
 * Start-up: the wave amplitude is ramped over RAMP_PERIODS periods.
 * Without this the first time steps see a full-amplitude wave hitting
 * still water and the run often diverges or the surface splashes at the
 * inlet.
 *
 * Volume fraction: sharp 0/1 with a linear smoothing band of width
 * VOF_SMOOTH (m) around eta. Set it to about one inlet cell height.
 *
 * Reflections: use damping_zone_source.c at the far end of the tank.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"
#include "wave_theory.h"

/* ---- user parameters ------------------------------------------------ */
#define WAVE_THEORY   2       /* 1 = Airy, 2 = Stokes 2nd order         */
#define WAVE_H        0.10    /* m, wave height (crest to trough)        */
#define WAVE_T        1.50    /* s, wave period                          */
#define DEPTH_D       1.00    /* m, still water depth                    */
#define Y_SWL         0.00    /* m, y coordinate of still water level    */
#define GRAV          9.81    /* m/s^2                                   */
#define RAMP_PERIODS  2.0     /* periods over which amplitude ramps in   */
#define VOF_SMOOTH    0.01    /* m, smoothing band for volume fraction   */
/* --------------------------------------------------------------------- */

static int    init_done = 0;
static double k_wave, w_wave, A_wave;

static void wave_init(void)
{
    if (init_done) return;
    w_wave = 2.0 * M_PI / WAVE_T;
    k_wave = wt_wavenumber(w_wave, DEPTH_D, GRAV);
    A_wave = 0.5 * WAVE_H;
    init_done = 1;
    if (UDF_IS_WRITER)
        Message("wave_inlet_regular: H=%g T=%g d=%g -> k=%g L=%g Ur=%g\n",
                WAVE_H, WAVE_T, DEPTH_D, k_wave, 2.0 * M_PI / k_wave,
                wt_ursell(WAVE_H, k_wave, DEPTH_D));
}

static double ramp(double t)
{
    double tr = RAMP_PERIODS * WAVE_T;
    if (tr <= 0.0 || t >= tr) return 1.0;
    if (t <= 0.0) return 0.0;
    return 0.5 * (1.0 - cos(M_PI * t / tr));
}

/* evaluate eta, u, v at (x, y) at the current time, ramp applied */
static double wave_at(double x, double y, double t, double *u, double *v)
{
    double eta, rr = ramp(t);
    double yy = y - Y_SWL;
    double A  = rr * A_wave;

#if WAVE_THEORY == 1
    eta = wt_airy(A, k_wave, w_wave, DEPTH_D, x, yy, t, u, v);
#else
    eta = wt_stokes2(A, k_wave, w_wave, DEPTH_D, x, yy, t, u, v);
#endif
    if (yy > eta) { *u = 0.0; *v = 0.0; }
    return eta;
}

DEFINE_PROFILE(inlet_u, t, i)
{
    face_t f;
    real x[ND_ND];
    double u, v, tt = CURRENT_TIME;

    wave_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        wave_at(x[0], x[1], tt, &u, &v);
        F_PROFILE(f, t, i) = u;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_v, t, i)
{
    face_t f;
    real x[ND_ND];
    double u, v, tt = CURRENT_TIME;

    wave_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        wave_at(x[0], x[1], tt, &u, &v);
        F_PROFILE(f, t, i) = v;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_vof, t, i)
{
    face_t f;
    real x[ND_ND];
    double u, v, eta, a, tt = CURRENT_TIME;

    wave_init();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        eta = wave_at(x[0], x[1], tt, &u, &v);
        /* fraction of water: 1 well below eta, 0 well above, linear in between */
        a = 0.5 + (eta - (x[1] - Y_SWL)) / (VOF_SMOOTH > 0.0 ? VOF_SMOOTH : 1.0e-9);
        if (a > 1.0) a = 1.0;
        if (a < 0.0) a = 0.0;
        F_PROFILE(f, t, i) = a;
    }
    end_f_loop(f, t)
}
