/*
 * wave_theory.h
 *
 * Pure C wave theory. No Fluent macros here, so this header can be
 * compiled and tested with plain gcc (see tests/test_wave_theory.c)
 * and then included unchanged by the Fluent UDFs in wavetank/.
 *
 * Conventions
 *   - y is vertical, positive upward, y = 0 at the still water level (SWL).
 *   - d is the still water depth (bed at y = -d).
 *   - x is the wave propagation direction.
 *   - theta = k*x - omega*t is the wave phase.
 *   - A = H/2 is the wave amplitude.
 *
 * References
 *   Dean & Dalrymple, "Water Wave Mechanics for Engineers and Scientists",
 *   chapters 3 (linear theory) and 11 (Stokes second order).
 *   Goda, "Random Seas and Design of Maritime Structures", JONSWAP form.
 */

#ifndef WAVE_THEORY_H
#define WAVE_THEORY_H

#include <math.h>

#ifndef WT_PI
#define WT_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Dispersion relation: omega^2 = g k tanh(k d). Solve for k by Newton. */
/* ------------------------------------------------------------------ */
static double wt_wavenumber(double omega, double d, double g)
{
    double k, f, df, k0;
    int it;

    /* deep water start guess, then Newton */
    k = omega * omega / g;
    if (k * d < 1.0e-6)
        k = omega / sqrt(g * d); /* shallow guess if very shallow */

    for (it = 0; it < 100; it++)
    {
        double th = tanh(k * d);
        double ch = cosh(k * d);
        f  = g * k * th - omega * omega;
        df = g * th + g * k * d / (ch * ch);
        k0 = k;
        k  = k - f / df;
        if (k <= 0.0)
            k = 0.5 * k0;
        if (fabs(k - k0) < 1.0e-12 * k)
            break;
    }
    return k;
}

/* ------------------------------------------------------------------ */
/* Linear (Airy) wave.                                                 */
/* Returns surface elevation eta at (x, t). u, v are horizontal and     */
/* vertical velocity at (x, y, t). y is clamped to the free surface     */
/* if y > eta so that a value just above the surface is not blown up.  */
/* ------------------------------------------------------------------ */
static double wt_airy(double A, double k, double omega, double d,
                      double x, double y, double t,
                      double *u, double *v)
{
    double theta = k * x - omega * t;
    double eta   = A * cos(theta);
    double yy    = (y > eta) ? eta : y;
    double sh    = sinh(k * d);

    *u = A * omega * cosh(k * (yy + d)) / sh * cos(theta);
    *v = A * omega * sinh(k * (yy + d)) / sh * sin(theta);
    return eta;
}

/* ------------------------------------------------------------------ */
/* Stokes second order wave (Dean & Dalrymple eq. 11.35, 11.37).        */
/* ------------------------------------------------------------------ */
static double wt_stokes2(double A, double k, double omega, double d,
                         double x, double y, double t,
                         double *u, double *v)
{
    double theta = k * x - omega * t;
    double sh    = sinh(k * d);
    double ch    = cosh(k * d);
    double sh4   = sh * sh * sh * sh;

    double eta = A * cos(theta)
               + (k * A * A / 4.0) * ch / (sh * sh * sh)
                 * (2.0 + cosh(2.0 * k * d)) * cos(2.0 * theta);

    double yy = (y > eta) ? eta : y;

    *u = A * omega * cosh(k * (yy + d)) / sh * cos(theta)
       + 0.75 * A * A * omega * k * cosh(2.0 * k * (yy + d)) / sh4
         * cos(2.0 * theta);

    *v = A * omega * sinh(k * (yy + d)) / sh * sin(theta)
       + 0.75 * A * A * omega * k * sinh(2.0 * k * (yy + d)) / sh4
         * sin(2.0 * theta);

    return eta;
}

/* ------------------------------------------------------------------ */
/* Ursell number H L^2 / d^3. Stokes theory is usually considered      */
/* valid for Ur < about 26 (deep/intermediate water). Report it.        */
/* ------------------------------------------------------------------ */
static double wt_ursell(double H, double k, double d)
{
    double L = 2.0 * WT_PI / k;
    return H * L * L / (d * d * d);
}

/* ------------------------------------------------------------------ */
/* JONSWAP spectrum in the Goda form (Goda 2000, eq. 2.12).             */
/* S(f) = beta_J Hs^2 fp^4 f^-5 exp(-1.25 (fp/f)^4) gamma^r             */
/* r = exp( -(f/fp - 1)^2 / (2 sigma^2) ), sigma = 0.07 (f<=fp), 0.09   */
/* beta_J = 0.0624 (1.094 - 0.01915 ln gamma)                          */
/*          / (0.230 + 0.0336 gamma - 0.185 (1.9 + gamma)^-1)          */
/* With this normalisation 4 sqrt(m0) reproduces Hs closely.           */
/* ------------------------------------------------------------------ */
static double wt_jonswap(double f, double Hs, double fp, double gamma)
{
    double sigma, r, betaJ, S;

    if (f <= 0.0)
        return 0.0;

    sigma = (f <= fp) ? 0.07 : 0.09;
    r     = exp(-(f / fp - 1.0) * (f / fp - 1.0) / (2.0 * sigma * sigma));
    betaJ = 0.0624 * (1.094 - 0.01915 * log(gamma))
          / (0.230 + 0.0336 * gamma - 0.185 / (1.9 + gamma));

    S = betaJ * Hs * Hs * pow(fp, 4.0) * pow(f, -5.0)
      * exp(-1.25 * pow(fp / f, 4.0)) * pow(gamma, r);
    return S;
}

/* ------------------------------------------------------------------ */
/* Tiny deterministic LCG so C and Python build identical phases.      */
/* Numerical Recipes constants. Returns a value in [0, 1).             */
/* ------------------------------------------------------------------ */
static unsigned int wt_lcg_state = 12345u;

static void wt_lcg_seed(unsigned int seed)
{
    wt_lcg_state = seed;
}

static double wt_lcg_uniform(void)
{
    wt_lcg_state = 1664525u * wt_lcg_state + 1013904223u; /* mod 2^32 */
    return (double)wt_lcg_state / 4294967296.0;
}

/* ------------------------------------------------------------------ */
/* Build N linear components from the JONSWAP spectrum between fmin    */
/* and fmax (equal frequency spacing). Fills arrays of length N:       */
/* amplitude a[i], wavenumber k[i], angular frequency w[i], phase p[i].*/
/* Returns m0 = sum S df, so the caller can check Hs = 4 sqrt(m0).      */
/* ------------------------------------------------------------------ */
static double wt_build_components(int N, double Hs, double Tp, double gamma,
                                  double fmin, double fmax, double d, double g,
                                  unsigned int seed,
                                  double *a, double *k, double *w, double *p)
{
    int i;
    double fp = 1.0 / Tp;
    double df = (fmax - fmin) / (double)N;
    double m0 = 0.0;

    wt_lcg_seed(seed);
    for (i = 0; i < N; i++)
    {
        double f = fmin + (i + 0.5) * df;
        double S = wt_jonswap(f, Hs, fp, gamma);
        a[i] = sqrt(2.0 * S * df);
        w[i] = 2.0 * WT_PI * f;
        k[i] = wt_wavenumber(w[i], d, g);
        p[i] = 2.0 * WT_PI * wt_lcg_uniform();
        m0  += S * df;
    }
    return m0;
}

/* ------------------------------------------------------------------ */
/* Irregular wave: superposition of N linear components.               */
/* ------------------------------------------------------------------ */
static double wt_irregular(int N, const double *a, const double *k,
                           const double *w, const double *p, double d,
                           double x, double y, double t,
                           double *u, double *v)
{
    int i;
    double eta = 0.0, uu = 0.0, vv = 0.0, yy;

    for (i = 0; i < N; i++)
        eta += a[i] * cos(k[i] * x - w[i] * t + p[i]);

    yy = (y > eta) ? eta : y;

    for (i = 0; i < N; i++)
    {
        double th = k[i] * x - w[i] * t + p[i];
        double sh = sinh(k[i] * d);
        uu += a[i] * w[i] * cosh(k[i] * (yy + d)) / sh * cos(th);
        vv += a[i] * w[i] * sinh(k[i] * (yy + d)) / sh * sin(th);
    }
    *u = uu;
    *v = vv;
    return eta;
}

#endif /* WAVE_THEORY_H */
