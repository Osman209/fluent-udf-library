/*
 * profile_pulsatile.c
 *
 * Pulsatile inlet for cardiovascular and other periodic internal flows.
 *
 *   pulsatile_flat      spatially uniform, waveform in time
 *   pulsatile_parabolic parabolic in space, waveform in time
 *
 * The waveform is a Fourier series, which is how measured flow-rate
 * waveforms are normally reported:
 *
 *   Q(t) = Q0 + sum_{n=1}^{N} [ An cos(n w t) + Bn sin(n w t) ],  w = 2 pi / T
 *
 * Fill FOUR_A and FOUR_B with your own harmonics. The default set is a
 * simple smooth pulse with a systolic peak; replace it with the
 * coefficients of your measured waveform.
 *
 * The profile returns a VELOCITY. If your coefficients are a flow rate
 * in m^3/s, set INLET_AREA and the code divides by it, so the mass flow
 * through the inlet follows the waveform exactly. If they are already
 * velocities, set INLET_AREA to 1.0.
 *
 * A note on physics: a parabolic profile in a pulsatile flow is only
 * correct at low Womersley number (alpha = R sqrt(w rho / mu) below
 * about 1). For a human aorta alpha is roughly 10 to 20, where the real
 * profile is much flatter with a thin oscillating layer near the wall.
 * At high alpha the flat profile is the closer approximation of the
 * two, or use a long enough entry length and let the solver develop it.
 * The code prints alpha at start-up so you can see where you are.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define PERIOD      0.8      /* s, cycle period (0.8 s = 75 bpm)        */
#define N_HARM      4        /* number of harmonics used                */
#define Q_MEAN      1.0      /* mean of the waveform (same units as A,B)*/
#define INLET_AREA  1.0      /* m^2, set to 1.0 if A,B are velocities   */
#define PIPE_R      0.01     /* m, radius (parabolic version, and alpha) */
#define AXIS_A0     0.0
#define AXIS_B0     0.0
#define CROSS_A     1
#define CROSS_B     2
#define RHO_F    1060.0      /* kg/m^3, for the Womersley number print  */
#define MU_F        0.0035   /* Pa s                                     */
/* --------------------------------------------------------------------- */

static const real FOUR_A[N_HARM] = { 0.33, -0.13,  0.02,  0.01};
static const real FOUR_B[N_HARM] = { 0.42,  0.05, -0.06,  0.01};

static int printed = 0;

static real waveform(real t)
{
    real w = 2.0 * M_PI / PERIOD, q = Q_MEAN;
    int n;
    for (n = 1; n <= N_HARM; n++)
        q += FOUR_A[n - 1] * cos(n * w * t) + FOUR_B[n - 1] * sin(n * w * t);
    return q;
}

static void print_alpha(void)
{
    real w, alpha;
    if (printed) return;
    printed = 1;
    w = 2.0 * M_PI / PERIOD;
    alpha = PIPE_R * sqrt(w * RHO_F / MU_F);
    Message("profile_pulsatile: period %g s, Womersley alpha = %.2f\n",
            (double)PERIOD, (double)alpha);
    if (alpha > 2.0)
        Message("profile_pulsatile: alpha > 2, the true profile is flatter "
                "than parabolic; consider pulsatile_flat.\n");
}

DEFINE_PROFILE(pulsatile_flat, t, i)
{
    face_t f;
    real u;
    print_alpha();
    u = waveform(CURRENT_TIME) / INLET_AREA;
    begin_f_loop(f, t) { F_PROFILE(f, t, i) = u; } end_f_loop(f, t)
}

DEFINE_PROFILE(pulsatile_parabolic, t, i)
{
    face_t f;
    real x[ND_ND], da, db, r2, R2 = PIPE_R * PIPE_R, umean, u;

    print_alpha();
    umean = waveform(CURRENT_TIME) / INLET_AREA;

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
        u  = 2.0 * umean * (1.0 - r2 / R2);
        F_PROFILE(f, t, i) = u;
    }
    end_f_loop(f, t)
}
