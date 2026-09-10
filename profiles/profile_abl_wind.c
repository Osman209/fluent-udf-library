/*
 * profile_abl_wind.c
 *
 * Atmospheric boundary layer inlet: wind speed, turbulence kinetic
 * energy and dissipation (or specific dissipation) as functions of
 * height. This is the standard inlet for building aerodynamics, wind
 * comfort, pollutant dispersion and wind-turbine siting.
 *
 * Two velocity laws, chosen by ABL_LAW:
 *   1 = log law     U(z) = (ustar / kappa) ln((z + z0)/z0)
 *   2 = power law   U(z) = U_ref (z/Z_REF)^ALPHA
 *
 * The log law is the one that is consistent with the k-epsilon
 * turbulence profiles below (Richards & Hoxey 1993):
 *   k(z)       = ustar^2 / sqrt(C_mu)                  (constant with height)
 *   epsilon(z) = ustar^3 / (kappa (z + z0))
 *   omega(z)   = epsilon / (C_mu k)                 (for k-omega / SST)
 *
 * ustar is computed from the reference speed at the reference height:
 *   ustar = kappa U_REF / ln((Z_REF + z0)/z0)
 *
 * Wind direction: WIND_DIR_DEG is the direction the wind blows TOWARD,
 * measured counter-clockwise from the +x axis. The two horizontal
 * components are set by inlet_u_abl and inlet_w_abl (or inlet_v_abl in
 * a Z-up case). Change VERT_DIR if your vertical axis is not y.
 *
 * Ground level is Z_GROUND, so height above ground is (coordinate -
 * Z_GROUND). Roughness z0 by terrain: 0.0002 m open sea, 0.03 m open
 * farmland, 0.3 m suburban, 1.0 m dense urban (Eurocode terrain
 * categories 0, II, III, IV).
 *
 * Hook the profiles on the velocity inlet: X-Velocity -> inlet_u_abl,
 * Z-Velocity -> inlet_w_abl, Turbulent Kinetic Energy -> inlet_k_abl,
 * Turbulent Dissipation Rate -> inlet_eps_abl (or Specific Dissipation
 * Rate -> inlet_omega_abl for k-omega / SST).
 *
 * To keep the profile from decaying along the domain you also need a
 * matching ground wall roughness (Ks about 9.793 z0 / Cs with Cs the
 * roughness constant) and a top boundary that does not fight the
 * profile. That is a case-setup matter, not a UDF matter.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define ABL_LAW       1        /* 1 = log law, 2 = power law            */
#define U_REF        10.0      /* m/s, reference wind speed             */
#define Z_REF        10.0      /* m, reference height above ground      */
#define Z0            0.03     /* m, aerodynamic roughness length       */
#define ALPHA         0.16     /* power-law exponent (law 2 only)       */
#define WIND_DIR_DEG  0.0      /* deg, direction wind blows toward      */
#define VERT_DIR      1        /* 0 = x, 1 = y, 2 = z: vertical axis    */
#define Z_GROUND      0.0      /* m, ground level                       */
#define KAPPA         0.42     /* von Karman constant                   */
#define C_MU          0.09
#define U_MIN         0.01     /* m/s, floor to avoid zero at the ground */
/* --------------------------------------------------------------------- */

static real u_star(void)
{
    return KAPPA * U_REF / log((Z_REF + Z0) / Z0);
}

static real abl_speed(real z)
{
    real h = z - Z_GROUND;
    if (h < 0.0) h = 0.0;
#if ABL_LAW == 1
    {
        real U = (u_star() / KAPPA) * log((h + Z0) / Z0);
        return (U < U_MIN) ? U_MIN : U;
    }
#else
    {
        real U = U_REF * pow((h < 1.0e-6 ? 1.0e-6 : h) / Z_REF, ALPHA);
        return (U < U_MIN) ? U_MIN : U;
    }
#endif
}

/* horizontal component indices: everything that is not the vertical */
static void horiz_dirs(int *a, int *b)
{
    if      (VERT_DIR == 1) { *a = 0; *b = 2; }
    else if (VERT_DIR == 2) { *a = 0; *b = 1; }
    else                    { *a = 1; *b = 2; }
}

DEFINE_PROFILE(inlet_u_abl, t, i)
{
    face_t f; real x[ND_ND]; int a, b;
    horiz_dirs(&a, &b);
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        F_PROFILE(f, t, i) = abl_speed(x[VERT_DIR]) * cos(WIND_DIR_DEG * M_PI / 180.0);
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_w_abl, t, i)
{
    face_t f; real x[ND_ND]; int a, b;
    horiz_dirs(&a, &b);
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        F_PROFILE(f, t, i) = abl_speed(x[VERT_DIR]) * sin(WIND_DIR_DEG * M_PI / 180.0);
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_k_abl, t, i)
{
    face_t f;
    real us = u_star();
    real k  = us * us / sqrt(C_MU);
    begin_f_loop(f, t)
    {
        F_PROFILE(f, t, i) = k;
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_eps_abl, t, i)
{
    face_t f; real x[ND_ND], h, us = u_star();
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        h = x[VERT_DIR] - Z_GROUND;
        if (h < 0.0) h = 0.0;
        F_PROFILE(f, t, i) = us * us * us / (KAPPA * (h + Z0));
    }
    end_f_loop(f, t)
}

DEFINE_PROFILE(inlet_omega_abl, t, i)
{
    face_t f; real x[ND_ND], h, us = u_star(), k, eps;
    k = us * us / sqrt(C_MU);
    begin_f_loop(f, t)
    {
        F_CENTROID(x, f, t);
        h = x[VERT_DIR] - Z_GROUND;
        if (h < 0.0) h = 0.0;
        eps = us * us * us / (KAPPA * (h + Z0));
        F_PROFILE(f, t, i) = eps / (C_MU * k);
    }
    end_f_loop(f, t)
}
