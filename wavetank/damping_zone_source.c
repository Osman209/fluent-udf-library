/*
 * damping_zone_source.c
 *
 * Numerical beach: a momentum sink that grows smoothly from zero at
 * X_START to its full value at X_END, absorbing outgoing waves so they
 * do not reflect from the outlet.
 *
 *   S_i = - rho * C(x) * u_i,     C(x) = C0 * ((x - X_START)/(X_END - X_START))^POW
 *
 * Hook: Cell Zone Conditions > fluid zone > Source Terms > X Momentum
 * and Y Momentum (mixture phase for VOF), pick x_mom_damp / y_mom_damp.
 * The implicit part dS/du = -rho C is returned so the solver stays
 * stable at large C0.
 *
 * Sizing rules that work in practice (from the wave-tank literature
 * and our own runs):
 *   - beach length 1.5 to 2 wavelengths
 *   - C0 of order 5 to 20 s^-1 for waves with T of about 1 to 2 s;
 *     too strong a beach reflects like a wall, too weak lets the wave
 *     through. Check the reflection with two probes in the working
 *     section (validation/wave_probe_analysis.py).
 *   - POW 2 (quadratic ramp) is the standard choice.
 *
 * Apply it to the mixture momentum equations. It damps both phases;
 * that is what you want (air velocity above a decaying wave should die
 * too). To damp only the vertical velocity, hook y_mom_damp alone.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define X_START  20.0   /* m, beach begins                              */
#define X_END    25.0   /* m, beach ends (outlet)                       */
#define C0       10.0   /* 1/s, maximum damping coefficient             */
#define POW       2.0   /* ramp exponent                                */
/* --------------------------------------------------------------------- */

static real damp_coeff(real x)
{
    real s;
    if (x <= X_START) return 0.0;
    if (x >= X_END)   return C0;
    s = (x - X_START) / (X_END - X_START);
    return C0 * pow(s, POW);
}

DEFINE_SOURCE(x_mom_damp, c, t, dS, eqn)
{
    real xc[ND_ND], C, rho, src;
    C_CENTROID(xc, c, t);
    C   = damp_coeff(xc[0]);
    rho = C_R(c, t);
    src = -rho * C * C_U(c, t);
    dS[eqn] = -rho * C;
    return src;
}

DEFINE_SOURCE(y_mom_damp, c, t, dS, eqn)
{
    real xc[ND_ND], C, rho, src;
    C_CENTROID(xc, c, t);
    C   = damp_coeff(xc[0]);
    rho = C_R(c, t);
    src = -rho * C * C_V(c, t);
    dS[eqn] = -rho * C;
    return src;
}

#if ND_ND == 3
DEFINE_SOURCE(z_mom_damp, c, t, dS, eqn)
{
    real xc[ND_ND], C, rho, src;
    C_CENTROID(xc, c, t);
    C   = damp_coeff(xc[0]);
    rho = C_R(c, t);
    src = -rho * C * C_W(c, t);
    dS[eqn] = -rho * C;
    return src;
}
#endif
