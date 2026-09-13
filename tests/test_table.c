/*
 * test_table.c
 * Checks udf_read_table and udf_interp1 from common/udf_common.h with
 * plain gcc: comments, mixed separators, malformed lines, clamping and
 * interior interpolation.
 */
#include <stdio.h>
#include <math.h>
#include "udf_common.h"

static int fails = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); fails++; } else printf("ok:   %s\n", msg); } while (0)

int main(void)
{
    const char *fname = "table_tmp.txt";
    FILE *fp = fopen(fname, "w");
    double *d = NULL, ts[4], vx[4];
    int n, r;

    fprintf(fp, "# t vx vy vz wx wy wz\n");
    fprintf(fp, "0.0, 0, 0, 0, 0, 0, 0\n");
    fprintf(fp, "1.0\t2.0 0 0 0 0 0.5\n");

    fprintf(fp, "2.0 4.0 0 0 0 0 1.0\n");
    fprintf(fp, "\n");
    fprintf(fp, "3.0 4.0 0 0 0 0 1.0\n");
    fclose(fp);

    n = udf_read_table(fname, 7, &d);
    CHECK(n == 4, "4 valid rows read, comments and blank lines skipped");
    for (r = 0; r < n; r++) { ts[r] = d[r * 7]; vx[r] = d[r * 7 + 1]; }

    CHECK(fabs(udf_interp1(ts, vx, n, 0.5) - 1.0) < 1e-12, "interp at 0.5 -> 1.0");
    CHECK(fabs(udf_interp1(ts, vx, n, 1.5) - 3.0) < 1e-12, "interp at 1.5 -> 3.0");
    CHECK(fabs(udf_interp1(ts, vx, n, -5.0) - 0.0) < 1e-12, "clamp below -> 0.0");
    CHECK(fabs(udf_interp1(ts, vx, n, 99.0) - 4.0) < 1e-12, "clamp above -> 4.0");
    CHECK(fabs(udf_interp1(ts, vx, n, 2.0) - 4.0) < 1e-12, "exact node -> 4.0");
    free(d); d = NULL;
    CHECK(udf_read_table("does_not_exist.txt", 7, &d) == -1, "missing file -> -1");

    remove(fname);
    printf("%s\n", fails ? "TESTS FAILED" : "all table checks passed");
    return fails ? 1 : 0;
}
