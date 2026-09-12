/*
 * property_from_table.c
 *
 * Material properties interpolated from a data file, instead of a
 * polynomial typed into the Materials panel. This covers the case that
 * comes up most often: you have measured or tabulated data for a fluid
 * and you want Fluent to use it directly.
 *
 * Two forms:
 *
 *   1D, property against temperature. File "prop_T.txt":
 *       # T [K]   value
 *       273.15    0.001792
 *       283.15    0.001307
 *       ...
 *   Functions: density_T, viscosity_T, conductivity_T, cp_T.
 *
 *   2D, property against temperature and pressure. File "prop_TP.txt",
 *   a regular grid written row by row:
 *       # nT nP
 *       # T values, then P values, then nT*nP property values
 *   Function: density_TP. Bilinear interpolation, clamped at the edges.
 *
 * The file names are set per property below, so you can use several at
 * once (one file for viscosity, one for conductivity, and so on) in a
 * single compiled library.
 *
 * DEFINE_SPECIFIC_HEAT is different from the other property macros: it
 * must also return the sensible enthalpy h(T) = integral of cp dT from
 * the reference temperature, and Fluent uses that enthalpy in the
 * energy equation. The code integrates the table with the trapezium
 * rule once at load time, so h is consistent with the cp you supply.
 * Getting this wrong is the usual reason a cp UDF gives a wrong
 * temperature field while cp itself looks right.
 *
 * Hook: Materials panel, set the property to "user-defined" and pick
 * the function.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define FILE_RHO    "rho_T.txt"
#define FILE_MU     "mu_T.txt"
#define FILE_K      "k_T.txt"
#define FILE_CP     "cp_T.txt"
#define FILE_RHO_TP "rho_TP.txt"
#define T_REF       298.15     /* K, reference temperature for enthalpy */
/* --------------------------------------------------------------------- */

typedef struct {
    double *x, *y, *h;   /* h used only for cp: integral of cp dT       */
    int n, loaded;
    const char *fname;
} Table1D;

static Table1D tb_rho = {NULL, NULL, NULL, 0, 0, FILE_RHO};
static Table1D tb_mu  = {NULL, NULL, NULL, 0, 0, FILE_MU};
static Table1D tb_k   = {NULL, NULL, NULL, 0, 0, FILE_K};
static Table1D tb_cp  = {NULL, NULL, NULL, 0, 0, FILE_CP};

static void t1_load(Table1D *tb, int build_enthalpy)
{
    double *raw = NULL;
    int r, n;

    if (tb->loaded) return;
    tb->loaded = 1;

    n = udf_read_table(tb->fname, 2, &raw);
    if (n <= 0)
    {
        if (UDF_IS_WRITER)
            Message("property_from_table: cannot read %s\n", tb->fname);
        tb->n = 0;
        return;
    }
    tb->x = (double *)malloc(sizeof(double) * n);
    tb->y = (double *)malloc(sizeof(double) * n);
    for (r = 0; r < n; r++) { tb->x[r] = raw[2 * r]; tb->y[r] = raw[2 * r + 1]; }
    tb->n = n;
    free(raw);

    if (build_enthalpy)
    {
        /* h(x_i) = integral from T_REF to x_i of y dx, trapezium rule.
         * Built by cumulating from the first point then shifting so
         * that h(T_REF) = 0. */
        double *c = (double *)malloc(sizeof(double) * n);
        double h_ref;
        c[0] = 0.0;
        for (r = 1; r < n; r++)
            c[r] = c[r - 1] + 0.5 * (tb->y[r] + tb->y[r - 1]) * (tb->x[r] - tb->x[r - 1]);
        h_ref = udf_interp1(tb->x, c, n, T_REF);
        for (r = 0; r < n; r++) c[r] -= h_ref;
        tb->h = c;
    }

    if (UDF_IS_WRITER)
        Message("property_from_table: %s, %d rows, x = %g .. %g, y = %g .. %g\n",
                tb->fname, n, tb->x[0], tb->x[n - 1], tb->y[0], tb->y[n - 1]);
}

static double t1_value(Table1D *tb, double x, int build_enthalpy, double fallback)
{
    t1_load(tb, build_enthalpy);
    if (tb->n == 0) return fallback;
    return udf_interp1(tb->x, tb->y, tb->n, x);
}

DEFINE_PROPERTY(density_T, c, t)
{
    return (real)t1_value(&tb_rho, C_T(c, t), 0, 998.2);
}

DEFINE_PROPERTY(viscosity_T, c, t)
{
    return (real)t1_value(&tb_mu, C_T(c, t), 0, 1.003e-3);
}

DEFINE_PROPERTY(conductivity_T, c, t)
{
    return (real)t1_value(&tb_k, C_T(c, t), 0, 0.6);
}

/*
 * Specific heat. Fluent asks for cp(T) AND the sensible enthalpy h(T)
 * measured from Tref. Both come from the same table, so they cannot
 * disagree.
 */
DEFINE_SPECIFIC_HEAT(cp_T, T, Tref, h, yi)
{
    real cp;
    t1_load(&tb_cp, 1);
    if (tb_cp.n == 0) { *h = 4182.0 * (T - Tref); return 4182.0; }
    cp = (real)udf_interp1(tb_cp.x, tb_cp.y, tb_cp.n, T);
    *h = (real)(udf_interp1(tb_cp.x, tb_cp.h, tb_cp.n, T)
              - udf_interp1(tb_cp.x, tb_cp.h, tb_cp.n, Tref));
    return cp;
}

/* ---- 2D table: property(T, P) --------------------------------------- */

static double *g_T = NULL, *g_P = NULL, *g_V = NULL;
static int g_nT = 0, g_nP = 0, g_loaded = 0;

static void t2_load(void)
{
    FILE *fp;
    int i, j;

    if (g_loaded) return;
    g_loaded = 1;

    fp = fopen(FILE_RHO_TP, "r");
    if (!fp)
    {
        if (UDF_IS_WRITER)
            Message("property_from_table: cannot read %s\n", FILE_RHO_TP);
        return;
    }
    if (fscanf(fp, "%d %d", &g_nT, &g_nP) != 2 || g_nT < 2 || g_nP < 2)
    {
        if (UDF_IS_WRITER)
            Message("property_from_table: bad header in %s "
                    "(expected: nT nP)\n", FILE_RHO_TP);
        g_nT = g_nP = 0; fclose(fp); return;
    }
    g_T = (double *)malloc(sizeof(double) * g_nT);
    g_P = (double *)malloc(sizeof(double) * g_nP);
    g_V = (double *)malloc(sizeof(double) * g_nT * g_nP);
    for (i = 0; i < g_nT; i++) if (fscanf(fp, "%lf", &g_T[i]) != 1) goto bad;
    for (j = 0; j < g_nP; j++) if (fscanf(fp, "%lf", &g_P[j]) != 1) goto bad;
    for (i = 0; i < g_nT; i++)
        for (j = 0; j < g_nP; j++)
            if (fscanf(fp, "%lf", &g_V[i * g_nP + j]) != 1) goto bad;
    fclose(fp);
    for (i = 1; i < g_nT; i++)
        if (g_T[i] <= g_T[i - 1])
        {
            if (UDF_IS_WRITER)
                Message("property_from_table: %s, temperatures must increase "
                        "(row %d is %g after %g). Table ignored.\n",
                        FILE_RHO_TP, i, g_T[i], g_T[i - 1]);
            g_nT = g_nP = 0;
            return;
        }
    for (j = 1; j < g_nP; j++)
        if (g_P[j] <= g_P[j - 1])
        {
            if (UDF_IS_WRITER)
                Message("property_from_table: %s, pressures must increase "
                        "(row %d is %g after %g). Table ignored.\n",
                        FILE_RHO_TP, j, g_P[j], g_P[j - 1]);
            g_nT = g_nP = 0;
            return;
        }
    if (UDF_IS_WRITER)
        Message("property_from_table: %s, grid %d x %d, T = %g .. %g, P = %g .. %g\n",
                FILE_RHO_TP, g_nT, g_nP, g_T[0], g_T[g_nT - 1], g_P[0], g_P[g_nP - 1]);
    return;
bad:
    if (UDF_IS_WRITER)
        Message("property_from_table: %s ended early, expected %d values\n",
                FILE_RHO_TP, g_nT * g_nP);
    g_nT = g_nP = 0;
    fclose(fp);
}

static int find_idx(const double *a, int n, double v)
{
    int lo = 0, hi = n - 1, mid;
    if (v <= a[0]) return 0;
    if (v >= a[n - 1]) return n - 2;
    while (hi - lo > 1) { mid = (lo + hi) / 2; if (a[mid] <= v) lo = mid; else hi = mid; }
    return lo;
}

DEFINE_PROPERTY(density_TP, c, t)
{
    double T, P, ft, fp2, dT, dP, v00, v01, v10, v11;
    int i, j;

    t2_load();
    if (g_nT == 0) return 998.2;

    T = C_T(c, t);
    P = C_P(c, t) + RP_Get_Real("operating-pressure");

    i = find_idx(g_T, g_nT, T);
    j = find_idx(g_P, g_nP, P);

    /* A table with a repeated temperature or pressure gives a zero
     * spacing here. Dividing by it yields inf or NaN, which then spreads
     * silently through the density field rather than stopping the run. */
    dT  = g_T[i + 1] - g_T[i];
    dP  = g_P[j + 1] - g_P[j];
    ft  = (dT != 0.0) ? (T - g_T[i]) / dT : 0.0;
    fp2 = (dP != 0.0) ? (P - g_P[j]) / dP : 0.0;
    if (ft  < 0.0) ft  = 0.0; if (ft  > 1.0) ft  = 1.0;
    if (fp2 < 0.0) fp2 = 0.0; if (fp2 > 1.0) fp2 = 1.0;

    v00 = g_V[i * g_nP + j];
    v01 = g_V[i * g_nP + j + 1];
    v10 = g_V[(i + 1) * g_nP + j];
    v11 = g_V[(i + 1) * g_nP + j + 1];

    return (real)((1 - ft) * ((1 - fp2) * v00 + fp2 * v01)
                +      ft  * ((1 - fp2) * v10 + fp2 * v11));
}
