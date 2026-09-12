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

#ifdef UDF_INSIDE_FLUENT
  /* inside Fluent: udf.h already defines RP_HOST, RP_NODE, I_AM_NODE_ZERO_P.
     Every UDF in this library does
         #include "udf.h"
         #define UDF_INSIDE_FLUENT
         #include "udf_common.h"
     Add common/udf_common.h and common/wave_theory.h under "Header Files"
     in the Compiled UDFs panel, or copy them next to the .c file.        */
  #if RP_NODE
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
/*     static real last_t = -1.0;                                       */
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
/* Read a numeric table. Lines starting with # are comments. Values     */
/* separated by commas, spaces or tabs. Returns number of rows read,    */
/* or -1 if the file cannot be opened. data is allocated with malloc    */
/* and laid out row-major: data[row*ncol + col].                        */
/* ------------------------------------------------------------------ */
static int udf_read_table(const char *fname, int ncol, double **data)
{
    FILE *fp;
    char line[1024];
    int rows = 0, cap = 256, c;
    double *d;

    fp = fopen(fname, "r");
    if (!fp) return -1;

    d = (double *)malloc(sizeof(double) * cap * ncol);
    if (!d) { fclose(fp); return -1; }

    while (fgets(line, sizeof(line), fp))
    {
        char *tok, *p = line;
        double vals[64];
        int got = 0;

        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r' || *p == '\0') continue;

        tok = strtok(p, ", \t\r\n");
        while (tok && got < 64)
        {
            vals[got++] = atof(tok);
            tok = strtok(NULL, ", \t\r\n");
        }
        if (got < ncol) continue; /* skip malformed line */

        if (rows == cap)
        {
            cap *= 2;
            d = (double *)realloc(d, sizeof(double) * cap * ncol);
            if (!d) { fclose(fp); return -1; }
        }
        for (c = 0; c < ncol; c++) d[rows * ncol + c] = vals[c];
        rows++;
    }
    fclose(fp);
    *data = d;
    return rows;
}

#endif /* UDF_COMMON_H */
