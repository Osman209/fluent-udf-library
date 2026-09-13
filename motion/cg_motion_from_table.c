/*
 * cg_motion_from_table.c
 *
 * Prescribed rigid-body motion read from a text table:
 *
 *     # t [s]   vx  vy  vz [m/s]   wx  wy  wz [rad/s]
 *     0.0       0   0   0          0   0   0
 *     0.5       1   0   0          0   0   0.2
 *     ...
 *
 * Linear interpolation between rows, clamped outside the table.
 * The file is read once, on the first call, by every process (each
 * compute node has access to the working directory in a normal Fluent
 * installation). Put the file in the Fluent working directory or give
 * an absolute path in TABLE_FILE.
 *
 * Hook: Dynamic Mesh Zones > Rigid Body > Motion UDF.
 *
 * Comment lines start with #. Separators may be comma, space or tab.
 * The table reader and interpolator are tested outside Fluent by
 * tests/test_table.c.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define TABLE_FILE "motion_table.txt"
#define NCOL 7   /* t, vx, vy, vz, wx, wy, wz */
/* --------------------------------------------------------------------- */

static double *tab      = NULL;
static double *tab_t    = NULL;   /* time column, contiguous copy */
static double *tab_col  = NULL;   /* one contiguous column buffer   */
static int     tab_rows = 0;
static int     tab_ok   = 0;

static void load_table(void)
{
    int r;

    if (tab_ok) return;

    tab_rows = udf_read_table(TABLE_FILE, NCOL, &tab);
    if (tab_rows <= 0)
    {
        Message("cg_motion_from_table: could not read %s (rows=%d). "
                "Motion set to zero.\n", TABLE_FILE, tab_rows);
        tab_rows = 0;
        tab_ok = 1;
        return;
    }

    tab_t   = (double *)malloc(sizeof(double) * tab_rows);
    tab_col = (double *)malloc(sizeof(double) * tab_rows);
    if (!tab_t || !tab_col) {
        free(tab); free(tab_t); free(tab_col);
        tab = tab_t = tab_col = NULL; tab_rows = 0; tab_ok = 1;
        Message("motion table: allocation failed, motion disabled\n"); return;
    }
    for (r = 0; r < tab_rows; r++)
        tab_t[r] = tab[r * NCOL];

    if (UDF_IS_WRITER)
        Message("cg_motion_from_table: read %d rows from %s, t = %g .. %g s\n",
                tab_rows, TABLE_FILE, tab_t[0], tab_t[tab_rows - 1]);
    tab_ok = 1;
}

static double table_value(int col, double t)
{
    int r;
    if (tab_rows == 0) return 0.0;
    for (r = 0; r < tab_rows; r++)
        tab_col[r] = tab[r * NCOL + col];
    return udf_interp1(tab_t, tab_col, tab_rows, t);
}

DEFINE_CG_MOTION(motion_from_table, dt, vel, omega, time, dtime)
{
    load_table();

    NV_S(vel,   =, 0.0);
    NV_S(omega, =, 0.0);

    vel[0]   = table_value(1, time);
    vel[1]   = table_value(2, time);
    vel[2]   = table_value(3, time);
    omega[0] = table_value(4, time);
    omega[1] = table_value(5, time);
    omega[2] = table_value(6, time);
}

/* Call from Execute on Demand if you edit the table and want to reload
 * without restarting Fluent. */
DEFINE_ON_DEMAND(reload_motion_table)
{
    if (tab)     free(tab);
    if (tab_t)   free(tab_t);
    if (tab_col) free(tab_col);
    tab = tab_t = tab_col = NULL;
    tab_rows = 0;
    tab_ok = 0;
    load_table();
}
