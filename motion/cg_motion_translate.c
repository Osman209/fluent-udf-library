/*
 * cg_motion_translate.c
 *
 * Prescribed rigid-body translation with a smooth start.
 * Hook: Dynamic Mesh > Dynamic Mesh Zones > Rigid Body > Motion UDF.
 *
 * The body accelerates from rest to the target velocity V_TARGET over
 * T_RAMP seconds using a half-cosine ramp (zero acceleration at both
 * ends). An impulsive start (velocity jumping from 0 to V in one step)
 * is the usual cause of a diverging first time step in dynamic mesh
 * runs, so do not set T_RAMP to zero unless you know why.
 *
 * Velocity after the ramp is constant. Rotation is zero.
 *
 * Verified: the integrated displacement after T_RAMP equals
 * V_TARGET * T_RAMP / 2 (see validation/motion_reference.py).
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define V_TARGET_X   10.0   /* m/s, target velocity, x component */
#define V_TARGET_Y    0.0   /* m/s */
#define V_TARGET_Z    0.0   /* m/s */
#define T_RAMP        1.0   /* s, ramp duration from rest to target */
/* --------------------------------------------------------------------- */

static real ramp_factor(real t)
{
    if (T_RAMP <= 0.0) return 1.0;
    if (t >= T_RAMP)   return 1.0;
    if (t <= 0.0)      return 0.0;
    return 0.5 * (1.0 - cos(M_PI * t / T_RAMP));
}

DEFINE_CG_MOTION(translate_ramp, dt, vel, omega, time, dtime)
{
    real f = ramp_factor(time);

    NV_S(vel,   =, 0.0);
    NV_S(omega, =, 0.0);

    vel[0] = f * V_TARGET_X;
    vel[1] = f * V_TARGET_Y;
    vel[2] = f * V_TARGET_Z;
}
