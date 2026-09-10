/*
 * test_profiles.c
 *
 * Numerical checks of the profile formulas, compiled with plain gcc.
 * The formulas are duplicated here from the UDFs (they are inside
 * DEFINE_PROFILE bodies and cannot be called directly), so this test
 * guards the maths, and profiles_check.py guards it again with numpy.
 *
 *   1. ABL log law returns U_REF at Z_REF exactly
 *   2. ABL k-epsilon profile satisfies the Richards-Hoxey consistency
 *      relation:  epsilon = ustar^3 / (kappa (z+z0))  and  k = ustar^2/sqrt(Cmu)
 *      give a constant turbulent viscosity gradient consistent with
 *      the log law, i.e. nu_t = Cmu k^2/eps = kappa ustar (z+z0)
 *   3. parabolic pipe profile: area average equals U_MEAN
 *   4. power-law pipe profile: area average equals U_MEAN
 *   5. parabolic channel profile: average across the gap equals U_MEAN
 *   6. cp table enthalpy: h(T) from the trapezium integral matches the
 *      analytic integral for a linear cp, and h(Tref) = 0
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "udf_common.h"

static int fails = 0;
#define CHECK(c, m) do { if (!(c)) { printf("FAIL: %s\n", m); fails++; } else printf("ok:   %s\n", m); } while (0)

int main(void)
{
    const double kappa = 0.42, Cmu = 0.09, Uref = 10.0, Zref = 10.0, z0 = 0.03;
    double ustar, U, nut, z, k, eps, sum, area, r, dr, R = 0.05, Umean = 1.0, n = 7.0, umax;
    int i, N = 200000;

    /* 1. log law */
    ustar = kappa * Uref / log((Zref + z0) / z0);
    U = (ustar / kappa) * log((Zref + z0) / z0);
    printf("ustar = %.6f, U(Zref) = %.10f\n", ustar, U);
    CHECK(fabs(U - Uref) < 1e-10, "ABL log law returns U_REF at Z_REF");

    /* 2. Richards-Hoxey consistency: Cmu k^2/eps = kappa ustar (z+z0) */
    k = ustar * ustar / sqrt(Cmu);
    {
        double worst = 0.0;
        for (z = 0.0; z < 200.0; z += 0.5)
        {
            eps = ustar * ustar * ustar / (kappa * (z + z0));
            nut = Cmu * k * k / eps;
            worst = fmax(worst, fabs(nut - kappa * ustar * (z + z0)) / (kappa * ustar * (z + z0)));
        }
        printf("worst relative error in nu_t = Cmu k^2/eps vs kappa ustar (z+z0): %.2e\n", worst);
        CHECK(worst < 1e-12, "ABL k-epsilon profiles are Richards-Hoxey consistent");
    }

    /* 3. parabolic pipe: (1/A) int u dA = Umean */
    sum = 0.0; area = 0.0; dr = R / N;
    for (i = 0; i < N; i++)
    {
        r = (i + 0.5) * dr;
        U = 2.0 * Umean * (1.0 - r * r / (R * R));
        sum  += U * 2.0 * M_PI * r * dr;
        area += 2.0 * M_PI * r * dr;
    }
    printf("parabolic pipe area-average = %.10f\n", sum / area);
    CHECK(fabs(sum / area - Umean) < 1e-6, "parabolic pipe area-average = U_MEAN");

    /* 4. power law pipe */
    umax = Umean * ((n + 1.0) * (2.0 * n + 1.0)) / (2.0 * n * n);
    sum = 0.0; area = 0.0;
    for (i = 0; i < N; i++)
    {
        r = (i + 0.5) * dr;
        U = umax * pow(1.0 - r / R, 1.0 / n);
        sum  += U * 2.0 * M_PI * r * dr;
        area += 2.0 * M_PI * r * dr;
    }
    printf("power-law pipe area-average = %.10f (umax = %.6f)\n", sum / area, umax);
    CHECK(fabs(sum / area - Umean) < 1e-5, "power-law pipe area-average = U_MEAN");

    /* 5. parabolic channel: mean across -H..H */
    {
        double H = 0.05, dy = 2.0 * H / N, y;
        sum = 0.0;
        for (i = 0; i < N; i++)
        {
            y = -H + (i + 0.5) * dy;
            sum += 1.5 * Umean * (1.0 - y * y / (H * H)) * dy;
        }
        printf("parabolic channel average = %.10f\n", sum / (2.0 * H));
        CHECK(fabs(sum / (2.0 * H) - Umean) < 1e-6, "parabolic channel average = U_MEAN");
    }

    /* 6. cp enthalpy: cp(T) = a + b T, so h(T) = a(T-Tr) + b(T^2-Tr^2)/2 */
    {
        int M = 401, j;
        double a = 4000.0, b = 0.6, T0 = 273.15, T1 = 473.15, Tref = 298.15;
        double *x = malloc(M * sizeof(double)), *y = malloc(M * sizeof(double));
        double *c = malloc(M * sizeof(double)), href, T, h_num, h_ana, worst = 0.0;
        for (j = 0; j < M; j++) { x[j] = T0 + (T1 - T0) * j / (M - 1.0); y[j] = a + b * x[j]; }
        c[0] = 0.0;
        for (j = 1; j < M; j++) c[j] = c[j-1] + 0.5 * (y[j] + y[j-1]) * (x[j] - x[j-1]);
        href = udf_interp1(x, c, M, Tref);
        for (j = 0; j < M; j++) c[j] -= href;
        CHECK(fabs(udf_interp1(x, c, M, Tref)) < 1e-6, "cp enthalpy is zero at T_REF");
        for (T = T0; T <= T1; T += 5.0)
        {
            h_num = udf_interp1(x, c, M, T);
            h_ana = a * (T - Tref) + 0.5 * b * (T * T - Tref * Tref);
            worst = fmax(worst, fabs(h_num - h_ana) / (fabs(h_ana) + 1.0));
        }
        printf("cp enthalpy worst relative error vs analytic: %.2e\n", worst);
        CHECK(worst < 1e-6, "cp enthalpy integral matches the analytic integral");
        free(x); free(y); free(c);
    }

    printf("%s\n", fails ? "TESTS FAILED" : "all profile checks passed");
    return fails ? 1 : 0;
}
