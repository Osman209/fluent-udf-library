/*
 * udf_common.h
 *
 * Small helpers shared by all UDFs in this library.
 *
 *   UDF_IS_WRITER   - true on exactly one process so files are written once
 *   udf_interp1     - linear interpolation in a sorted table
 *   udf_read_table  - read a numeric CSV / whitespace table into memory
 *
 * The file reading and interpolation parts are plain C and are tested
 * outside Fluent by tests/test_table.c.
 *
 * Parallel notes
 *   In a parallel Fluent run there is one host process and N compute
 *   nodes. The host has no mesh data. Loops over cells or faces must
 *   therefore be wrapped in  #if !RP_HOST ... #endif.  Sums over faces
 *   must be reduced with PRF_GRSUM1 and then written by one node only,
 *   which is what UDF_IS_WRITER selects. In serial Fluent RP_HOST and
 *   RP_NODE are both 0 and I_AM_NODE_ZERO_P is true, so the same code
 *   works unchanged.
 */

#ifndef UDF_COMMON_H
#define UDF_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>

#ifdef UDF_INSIDE_FLUENT
  /* inside Fluent: udf.h already defines RP_HOST, RP_NODE, I_AM_NODE_ZERO_P.
     Every UDF in this library does
         #include "udf.h"
         #define UDF_INSIDE_FLUENT
         #include "udf_common.h"
     Add common/udf_common.h and common/wave_theory.h under "Header Files"
     in the Compiled UDFs panel, or copy them next to the .c file.        */
  #if RP_HOST
    #define UDF_IS_WRITER (0)
  #elif RP_NODE
    #define UDF_IS_WRITER (I_AM_NODE_ZERO_P)
  #else
    #define UDF_IS_WRITER (1)
  #endif
#else
  /* standalone compile for tests */
  #define UDF_IS_WRITER (1)
  #define UDF_MESSAGE printf
#endif

#ifndef UDF_MESSAGE
  #define UDF_MESSAGE Message
#endif

/* ------------------------------------------------------------------ */
/* Has the run restarted?                                              */
/*                                                                      */
/* A loaded library keeps its static variables across Initialize. So a  */
/* second run inherits the first run's state, and a "have I done this   */
/* yet" flag never fires again. The symptoms differ by what the flag    */
/* guarded:                                                             */
/*   - a log file is appended to instead of rewritten, and ends up with */
/*     two time histories interleaved                                   */
/*   - a "only log when the clock has advanced" test never passes,      */
/*     because the clock now starts behind where it stopped, so nothing */
/*     is logged at all                                                 */
/*   - worst: a reference position captured on the first call is kept,  */
/*     so a spring datum or a moment centre silently belongs to the     */
/*     previous run                                                     */
/*                                                                      */
/* The flow time going backwards is the reliable signal. Give each file */
/* its own static and call this once per step:                          */
/*                                                                      */
/*     static double last_t = -1.0;                                       */
/*     if (udf_restarted(CURRENT_TIME, &last_t)) { ...reset state... }  */
/*                                                                      */
/* It returns true on the very first call too, which is what you want:  */
/* "start fresh" covers both cases.                                     */
/* ------------------------------------------------------------------ */
static int udf_restarted(double now, double *last)
{
    /* True on the first call of this library load (*last still negative)
     * and whenever the clock has gone backwards. Testing only "now <
     * *last" misses the first call and the file never gets its header. */
    int fresh = (*last < 0.0) || (now < *last);
    *last = now;
    return fresh;
}

/* ------------------------------------------------------------------ */
/* Linear interpolation. xs must be increasing. Clamps at both ends.   */
/* ------------------------------------------------------------------ */
static double udf_interp1(const double *xs, const double *ys, int n, double x)
{
    int lo, hi, mid;

    if (n <= 0) return 0.0;
    if (n == 1 || x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];

    lo = 0; hi = n - 1;
    while (hi - lo > 1)
    {
        mid = (lo + hi) / 2;
        if (xs[mid] <= x) lo = mid; else hi = mid;
    }
    return ys[lo] + (ys[hi] - ys[lo]) * (x - xs[lo]) / (xs[hi] - xs[lo]);
}

/* ------------------------------------------------------------------ */
/* Strict numeric table: # comments and blank lines are accepted.
 * Exactly ncol finite values per row; first column strictly increasing.
 * Returns rows, -1 for I/O/allocation failure, -2 for invalid input.
 * On failure *data is NULL; the caller owns successful allocated data.
 * Text headers must start with #. Malformed rows reject the whole table.
 */
static int udf_read_table(const char *fname, int ncol, double **data)
{
    FILE *fp;
    char line[4096];
    int rows = 0, cap = 256, lineno = 0, c, result = -2;
    double *d = NULL;
    if (!data) return -2;
    *data = NULL;
    if (ncol < 1 || ncol > 64) return -2;
    fp = fopen(fname, "r");
    if (!fp) return -1;
    d = (double *)malloc(sizeof(double) * cap * ncol);
    if (!d) { fclose(fp); return -1; }
    while (fgets(line, sizeof(line), fp))
    {
        char *p = line, *end;
        double vals[64];
        lineno++;
        if (!strchr(line, '\n') && !feof(fp)) goto bad;
        while (isspace((unsigned char)*p)) p++;
        if (*p == '#' || !*p) continue;
        for (c = 0; c < ncol; c++)
        {
            errno = 0;
            vals[c] = strtod(p, &end);
            if (end == p || errno == ERANGE || !isfinite(vals[c])) goto bad;
            p = end;
            if (c + 1 < ncol)
            {
                int had_space = isspace((unsigned char)*p);
                while (isspace((unsigned char)*p)) p++;
                if (*p == ',') { p++; while (isspace((unsigned char)*p)) p++; }
                else if (!had_space) goto bad;
            }
        }
        while (isspace((unsigned char)*p)) p++;
        if (*p && *p != '#') goto bad;
        if (rows && vals[0] <= d[(rows - 1) * ncol]) goto bad;
        if (rows == cap)
        {
            double *grown;
            if (cap > INT_MAX / (2 * ncol)) { result = -1; goto bad; }
            cap *= 2;
            grown = (double *)realloc(d, sizeof(double) * (size_t)cap * ncol);
            if (!grown) { result = -1; goto bad; }
            d = grown;
        }
        for (c = 0; c < ncol; c++) d[rows * ncol + c] = vals[c];
        rows++;
    }
    if (ferror(fp)) { result = -1; goto bad; }
    fclose(fp);
    if (!rows) { free(d); return 0; }
    *data = d;
    return rows;
bad:
    if (UDF_IS_WRITER) UDF_MESSAGE("table %s: rejected at line %d (columns, finite values, or increasing axis)\n", fname, lineno);
    free(d);
    fclose(fp);
    return result;
}

/* Antiderivative of the piecewise-linear property, based at xs[0].
 * prefix[i] is its exact integral at knot i. Outside the table the
 * property is clamped, so the integral continues LINEARLY, not flat.
 */
static double udf_integral1(const double *xs, const double *ys,
                           const double *prefix, int n, double x)
{
    int lo = 0, hi = n - 1, mid;
    double dx, slope;
    if (n <= 0) return 0.0;
    if (n == 1 || x <= xs[0]) return ys[0] * (x - xs[0]);
    if (x >= xs[n - 1]) return prefix[n - 1] + ys[n - 1] * (x - xs[n - 1]);
    while (hi - lo > 1) { mid = (lo + hi) / 2; if (xs[mid] <= x) lo = mid; else hi = mid; }
    dx = x - xs[lo];
    slope = (ys[hi] - ys[lo]) / (xs[hi] - xs[lo]);
    return prefix[lo] + ys[lo] * dx + 0.5 * slope * dx * dx;
}
#endif /* UDF_COMMON_H */
