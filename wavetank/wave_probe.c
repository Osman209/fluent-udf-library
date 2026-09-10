/*
 * wave_probe.c
 *
 * Free-surface elevation at N_PROBE stations along the tank, written
 * every time step to wave_probes.csv. This is the numerical equivalent
 * of a wave gauge and is what you compare with theory to validate the
 * tank (validation/wave_probe_analysis.py: wave height, period, and
 * reflection coefficient from two probes).
 *
 * Method: in a vertical column of width PROBE_W centred on x_p, the
 * water volume is  sum(alpha_water * cell volume). Divided by the
 * column footprint (PROBE_W * SPAN_Z) it gives the water height above
 * the bed, and eta = Y_BED + height - Y_SWL. In 2D Fluent C_VOLUME is
 * the cell area times a unit depth, so SPAN_Z = 1.0. In 3D set SPAN_Z
 * to the tank width (the column then spans the whole width).
 *
 * This is an integral measure, so it is smooth and does not depend on
 * finding the alpha = 0.5 cell. Its resolution is the column width:
 * use PROBE_W about 2 to 4 cell widths. Cells are counted if their
 * centroid is inside the column, so pick PROBE_W as a multiple of the
 * local cell size to avoid a small bias.
 *
 * Parallel: begin_c_loop_int visits interior cells only (no ghost
 * cells), sums are reduced with PRF_GRSUM1, one process writes.
 *
 * WATER_PHASE is the index of the water phase in the Phases panel
 * (primary phase = 0). For "air primary, water secondary" it is 1.
 */

#include "udf.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define N_PROBE     4
#define PROBE_X     {5.0, 8.0, 10.0, 12.0}   /* m, station positions   */
#define PROBE_W     0.10                     /* m, column width        */
#define Y_BED      -1.0                      /* m, bed level           */
#define Y_SWL       0.0                      /* m, still water level   */
#define SPAN_Z      1.0                      /* 1.0 in 2D, width in 3D */
#define WATER_PHASE 1
#define LOG_FILE   "wave_probes.csv"
/* --------------------------------------------------------------------- */

static int first = 1;

DEFINE_EXECUTE_AT_END(wave_probes)
{
#if !RP_HOST
    Domain *mix = Get_Domain(1);
    Thread *t, *tw;
    cell_t c;
    real xc[ND_ND];
    real px[N_PROBE] = PROBE_X;
    real vol[N_PROBE];
    real eta[N_PROBE];
    int i;

    for (i = 0; i < N_PROBE; i++) vol[i] = 0.0;

    thread_loop_c(t, mix)
    {
        tw = THREAD_SUB_THREAD(t, WATER_PHASE);
        begin_c_loop_int(c, t)
        {
            C_CENTROID(xc, c, t);
            for (i = 0; i < N_PROBE; i++)
            {
                if (fabs(xc[0] - px[i]) <= 0.5 * PROBE_W)
                    vol[i] += C_VOF(c, tw) * C_VOLUME(c, t);
            }
        }
        end_c_loop_int(c, t)
    }

    for (i = 0; i < N_PROBE; i++)
    {
        vol[i] = PRF_GRSUM1(vol[i]);
        eta[i] = Y_BED + vol[i] / (PROBE_W * SPAN_Z) - Y_SWL;
    }

    if (UDF_IS_WRITER)
    {
        FILE *fp;
        if (first)
        {
            fp = fopen(LOG_FILE, "w");
            if (fp)
            {
                fprintf(fp, "t");
                for (i = 0; i < N_PROBE; i++) fprintf(fp, ",eta_x%g", px[i]);
                fprintf(fp, "\n");
                fclose(fp);
            }
            first = 0;
        }
        fp = fopen(LOG_FILE, "a");
        if (fp)
        {
            fprintf(fp, "%g", CURRENT_TIME);
            for (i = 0; i < N_PROBE; i++) fprintf(fp, ",%g", eta[i]);
            fprintf(fp, "\n");
            fclose(fp);
        }
    }
#endif
}
