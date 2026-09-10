/*
 * cg_motion_oscillate.c
 *
 * Prescribed harmonic heave and pitch, the standard oscillating foil /
 * cylinder / bridge-deck forced-motion case.
 *
 *   heave:  h(t)     = H0 * sin(2 pi f t + PHASE_H)
 *   pitch:  theta(t) = TH0 * sin(2 pi f t + PHASE_TH)
 *
 * DEFINE_CG_MOTION must return velocities, not positions, so the UDF
 * returns the time derivatives. Fluent integrates them. The pitch axis
 * is the "Center of Gravity Location" you enter in the Dynamic Mesh
 * Zones panel for this zone, and the rotation axis is chosen by which
 * omega component is nonzero (here z, i.e. the 2D case).
 *
 * Two functions are provided so you can hook heave only, pitch only,
 * or both (both = the combined function).
 *
 * The transient start is ramped over T_RAMP seconds to avoid the
 * impulsive first step. Set T_RAMP = 0 to disable.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define FREQ       1.0        /* Hz, oscillation frequency               */
#define H0         0.05       /* m,  heave amplitude                     */
#define PHASE_H    0.0        /* rad, heave phase                        */
#define TH0        (10.0*M_PI/180.0) /* rad, pitch amplitude (10 deg)    */
#define PHASE_TH   (M_PI/2.0) /* rad, pitch phase (90 deg lead is usual) */
#define HEAVE_DIR  1          /* 0 = x, 1 = y, 2 = z                      */
#define PITCH_AXIS 2          /* 0 = x, 1 = y, 2 = z                      */
#define T_RAMP     0.5        /* s, ramp-in of the amplitudes             */
/* --------------------------------------------------------------------- */

/* The ramp multiplies the AMPLITUDE, h(t) = r(t) H0 sin(w t + phi), and
 * the UDF returns the exact derivative  r' H0 sin + r H0 w cos.  This
 * way the motion after the ramp is exactly H0 sin(w t + phi) with no
 * drift of the mean position. (Ramping the velocity alone leaves a
 * permanent offset; validation/motion_reference.py shows the size.)  */
static real ramp_r(real t)
{
    if (T_RAMP <= 0.0) return 1.0;
    if (t >= T_RAMP)   return 1.0;
    if (t <= 0.0)      return 0.0;
    return 0.5 * (1.0 - cos(M_PI * t / T_RAMP));
}

static real ramp_dr(real t)
{
    if (T_RAMP <= 0.0) return 0.0;
    if (t >= T_RAMP)   return 0.0;
    if (t <= 0.0)      return 0.0;
    return 0.5 * (M_PI / T_RAMP) * sin(M_PI * t / T_RAMP);
}

static real harmonic_rate(real amp, real w, real phase, real t)
{
    return ramp_dr(t) * amp * sin(w * t + phase)
         + ramp_r(t)  * amp * w * cos(w * t + phase);
}

DEFINE_CG_MOTION(heave_pitch, dt, vel, omega, time, dtime)
{
    real w = 2.0 * M_PI * FREQ;

    NV_S(vel,   =, 0.0);
    NV_S(omega, =, 0.0);

    vel[HEAVE_DIR]    = harmonic_rate(H0,  w, PHASE_H,  time);
    omega[PITCH_AXIS] = harmonic_rate(TH0, w, PHASE_TH, time);
}

DEFINE_CG_MOTION(heave_only, dt, vel, omega, time, dtime)
{
    real w = 2.0 * M_PI * FREQ;

    NV_S(vel,   =, 0.0);
    NV_S(omega, =, 0.0);
    vel[HEAVE_DIR] = harmonic_rate(H0, w, PHASE_H, time);
}

DEFINE_CG_MOTION(pitch_only, dt, vel, omega, time, dtime)
{
    real w = 2.0 * M_PI * FREQ;

    NV_S(vel,   =, 0.0);
    NV_S(omega, =, 0.0);
    omega[PITCH_AXIS] = harmonic_rate(TH0, w, PHASE_TH, time);
}
