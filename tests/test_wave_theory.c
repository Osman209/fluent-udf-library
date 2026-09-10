/*
 * test_wave_theory.c
 * Compiles common/wave_theory.h with plain gcc and prints reference
 * values that validation/wave_theory_check.py recomputes independently
 * with numpy. Run via tests/run_tests.sh.
 *
 * Checks done here (exit code 1 on failure):
 *   1. dispersion: g k tanh(kd) - w^2 = 0 to 1e-9 relative
 *   2. deep and shallow limits of k
 *   3. Stokes 2nd order -> Airy as H -> 0
 *   4. JONSWAP: 4 sqrt(m0) within 3 % of Hs
 * Then it writes wave_ref.csv with (x, y, t, eta, u, v) samples and
 * jonswap_ref.csv with the component list for the Python comparison.
 */
#include <stdio.h>
#include <stdlib.h>
#include "wave_theory.h"

static int fails = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); fails++; } else printf("ok:   %s\n", msg); } while (0)

int main(void)
{
    double g = 9.81, d = 1.0, H = 0.10, T = 1.5;
    double w = 2.0 * WT_PI / T, k, res, kdeep, kshal;
    double u1, v1, u2, v2, e1, e2, x = 3.7, y = -0.3, t = 0.83;
    double A = 0.5 * H;
    int i;

    /* 1. dispersion */
    k = wt_wavenumber(w, d, g);
    res = fabs(g * k * tanh(k * d) - w * w) / (w * w);
    printf("k = %.12g  L = %.12g  residual = %.3g\n", k, 2 * WT_PI / k, res);
    CHECK(res < 1e-9, "dispersion residual < 1e-9");

    /* 2. limits */
    kdeep = wt_wavenumber(w, 1000.0, g);
    CHECK(fabs(kdeep - w * w / g) / kdeep < 1e-6, "deep water k = w^2/g");
    kshal = wt_wavenumber(2 * WT_PI / 60.0, 0.5, g);
    CHECK(fabs(kshal - (2 * WT_PI / 60.0) / sqrt(g * 0.5)) / kshal < 1e-2,
          "shallow water k = w/sqrt(gd) within 1 %");

    /* 3. Stokes -> Airy: the correction must scale as A^2, so halving
     * the amplitude must divide (stokes - airy) by 4 in eta, u and v. */
    {
        double dA1, dA2, ua, va, us, vs, ea, es;
        ea = wt_airy(A,   k, w, d, x, y, t, &ua, &va);
        es = wt_stokes2(A, k, w, d, x, y, t, &us, &vs);
        dA1 = fabs(es - ea) + fabs(us - ua) + fabs(vs - va);
        ea = wt_airy(A/2,   k, w, d, x, y, t, &ua, &va);
        es = wt_stokes2(A/2, k, w, d, x, y, t, &us, &vs);
        dA2 = fabs(es - ea) + fabs(us - ua) + fabs(vs - va);
        printf("Stokes-Airy difference: A -> %.3e, A/2 -> %.3e, ratio %.4f (expect 4)\n",
               dA1, dA2, dA1 / dA2);
        CHECK(fabs(dA1 / dA2 - 4.0) < 1e-3,
              "Stokes 2nd order correction scales as A^2 (reduces to Airy)");
    }

    /* reference samples */
    {
        FILE *fp = fopen("wave_ref.csv", "w");
        double xs[3] = {0.0, 1.3, 3.7}, ys[3] = {-0.9, -0.3, 0.0}, ts[3] = {0.0, 0.83, 1.2};
        int a, b, c;
        fprintf(fp, "# H=%g T=%g d=%g g=%g k=%.15g\n", H, T, d, g, k);
        fprintf(fp, "theory,x,y,t,eta,u,v\n");
        for (a = 0; a < 3; a++) for (b = 0; b < 3; b++) for (c = 0; c < 3; c++)
        {
            e1 = wt_airy(A, k, w, d, xs[a], ys[b], ts[c], &u1, &v1);
            fprintf(fp, "airy,%g,%g,%g,%.15g,%.15g,%.15g\n", xs[a], ys[b], ts[c], e1, u1, v1);
            e2 = wt_stokes2(A, k, w, d, xs[a], ys[b], ts[c], &u2, &v2);
            fprintf(fp, "stokes2,%g,%g,%g,%.15g,%.15g,%.15g\n", xs[a], ys[b], ts[c], e2, u2, v2);
        }
        fclose(fp);
    }
    printf("Ursell number for H=%g T=%g d=%g: %.3f\n", H, T, d, wt_ursell(H, k, d));

    /* 4. JONSWAP */
    {
        int N = 100;
        double Hs = 0.10, Tp = 1.5, gam = 3.3, fmin = 0.4 / Tp, fmax = 3.0 / Tp, m0;
        double *a = malloc(N * sizeof(double)), *kk = malloc(N * sizeof(double));
        double *ww = malloc(N * sizeof(double)), *pp = malloc(N * sizeof(double));
        FILE *fp;
        m0 = wt_build_components(N, Hs, Tp, gam, fmin, fmax, d, g, 12345u, a, kk, ww, pp);
        printf("JONSWAP: 4 sqrt(m0) = %.5f (Hs = %.5f)\n", 4 * sqrt(m0), Hs);
        CHECK(fabs(4 * sqrt(m0) - Hs) / Hs < 0.03, "4 sqrt(m0) within 3 % of Hs");
        fp = fopen("jonswap_ref.csv", "w");
        fprintf(fp, "# Hs=%g Tp=%g gamma=%g N=%d fmin=%.15g fmax=%.15g d=%g seed=12345\n",
                Hs, Tp, gam, N, fmin, fmax, d);
        fprintf(fp, "i,amplitude,wavenumber,omega,phase\n");
        for (i = 0; i < N; i++)
            fprintf(fp, "%d,%.15g,%.15g,%.15g,%.15g\n", i, a[i], kk[i], ww[i], pp[i]);
        fclose(fp);
        e1 = wt_irregular(N, a, kk, ww, pp, d, x, y, t, &u1, &v1);
        fp = fopen("irregular_ref.csv", "w");
        fprintf(fp, "x,y,t,eta,u,v\n%g,%g,%g,%.15g,%.15g,%.15g\n", x, y, t, e1, u1, v1);
        fclose(fp);
        free(a); free(kk); free(ww); free(pp);
    }

    printf("%s\n", fails ? "TESTS FAILED" : "all C wave-theory checks passed");
    return fails ? 1 : 0;
}
