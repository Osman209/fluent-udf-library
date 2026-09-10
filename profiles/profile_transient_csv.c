/*
 * profile_transient_csv.c
 *
 * Any boundary condition varying with time, read from a text file:
 *
 *     # t [s]   value
 *     0.0       0.0
 *     0.5       101325.0
 *     ...
 *
 * Linear interpolation between rows, clamped outside the table.
 * Use it for a pressure inlet or outlet that ramps, a velocity that
 * follows a measured signal, a temperature schedule, a mass flow that
 * changes with time, anything hooked through DEFINE_PROFILE.
 *
 * Three separate functions with three separate files are provided so
 * you can drive three boundaries at once (for example inlet velocity,
 * inlet temperature and outlet pressure) from one compiled library.
 *
 * Hook: on the boundary panel, set the quantity to "user-defined" and
 * pick transient_1 / transient_2 / transient_3.
 *
 * Note on how Fluent calls this: DEFINE_PROFILE is re-evaluated every
 * iteration, so the value follows the current flow time correctly in a
 * transient run. In a steady run CURRENT_TIME does not advance and the
 * profile stays at the first table value; that is expected.
 *
 * Common mistake this avoids: writing a table with a coarse time step
 * and expecting a smooth ramp. The interpolation is linear between the
 * rows you give, so put in enough rows for the shape you want.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define FILE_1 "transient_1.txt"
#define FILE_2 "transient_2.txt"
#define FILE_3 "transient_3.txt"
/* --------------------------------------------------------------------- */

typedef struct {
    double *t, *v;
    int n, loaded;
    const char *fname;
} TimeTable;

static TimeTable tt1 = {NULL, NULL, 0, 0, FILE_1};
static TimeTable tt2 = {NULL, NULL, 0, 0, FILE_2};
static TimeTable tt3 = {NULL, NULL, 0, 0, FILE_3};

static void tt_load(TimeTable *tt)
{
    double *raw = NULL;
    int r, n;

    if (tt->loaded) return;
    tt->loaded = 1;

    n = udf_read_table(tt->fname, 2, &raw);
    if (n <= 0)
    {
        if (UDF_IS_WRITER)
            Message("profile_transient_csv: cannot read %s, value fixed at 0\n",
                    tt->fname);
        tt->n = 0;
        return;
    }
    tt->t = (double *)malloc(sizeof(double) * n);
    tt->v = (double *)malloc(sizeof(double) * n);
    for (r = 0; r < n; r++) { tt->t[r] = raw[2 * r]; tt->v[r] = raw[2 * r + 1]; }
    tt->n = n;
    free(raw);

    if (UDF_IS_WRITER)
        Message("profile_transient_csv: %s, %d rows, t = %g .. %g s, "
                "value = %g .. %g\n", tt->fname, n, tt->t[0], tt->t[n - 1],
                tt->v[0], tt->v[n - 1]);
}

static double tt_value(TimeTable *tt, double t)
{
    tt_load(tt);
    if (tt->n == 0) return 0.0;
    return udf_interp1(tt->t, tt->v, tt->n, t);
}

DEFINE_PROFILE(transient_1, t, i)
{
    face_t f;
    real val = (real)tt_value(&tt1, CURRENT_TIME);
    begin_f_loop(f, t) { F_PROFILE(f, t, i) = val; } end_f_loop(f, t)
}

DEFINE_PROFILE(transient_2, t, i)
{
    face_t f;
    real val = (real)tt_value(&tt2, CURRENT_TIME);
    begin_f_loop(f, t) { F_PROFILE(f, t, i) = val; } end_f_loop(f, t)
}

DEFINE_PROFILE(transient_3, t, i)
{
    face_t f;
    real val = (real)tt_value(&tt3, CURRENT_TIME);
    begin_f_loop(f, t) { F_PROFILE(f, t, i) = val; } end_f_loop(f, t)
}

/* Edit a table while Fluent is open, then run this from Execute on
 * Demand to pick up the change without restarting. */
DEFINE_ON_DEMAND(reload_transient_tables)
{
    TimeTable *all[3]; int i;
    all[0] = &tt1; all[1] = &tt2; all[2] = &tt3;
    for (i = 0; i < 3; i++)
    {
        if (all[i]->t) free(all[i]->t);
        if (all[i]->v) free(all[i]->v);
        all[i]->t = all[i]->v = NULL;
        all[i]->n = 0;
        all[i]->loaded = 0;
        tt_load(all[i]);
    }
}
