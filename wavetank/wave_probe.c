/*
 * wave_probe.c
 *
 * Free-surface elevation at N_PROBE stations along the tank, written
 * every time step to wave_probes.csv.
 *
 * The first two stations are a PAIR: the reflection coefficient is
 * computed from them, and that method needs them between 0.05 and 0.45
 * wavelengths apart. The shipped positions are 1.0 m apart, which is
 * 0.30 of the 3.35 m default wavelength. If you change the wave, move
 * them: at 0.5 or 1.0 wavelengths the two probes see the same phase and
 * the incident and reflected waves cannot be told apart, and the answer
 * comes out wrong without looking wrong. This is the numerical equivalent
 * of a wave gauge and is what you compare with theory to validate the
 * tank (validation/wave_probe_analysis.py: wave height, period, and
 * reflection coefficient from two probes).
 *
 * Method: in a vertical column of width PROBE_W centred on x_p, both
 * the water volume sum(alpha * V) and the total volume sum(V) are
 * accumulated. The water fraction of the column is their ratio, and
 *
 *     eta = Y_BED + (sum alpha V / sum V) * (Y_TOP - Y_BED) - Y_SWL
 *
 * Dividing by the measured column volume rather than by the nominal
 * footprint PROBE_W is what makes this correct. A cell joins the column
 * when its centroid falls inside, so the width actually captured is a
 * whole number of cells and is almost never equal to PROBE_W. Dividing
 * by PROBE_W therefore reports the ratio of captured width to nominal
 * width, not the water level: in a 2D tank with dx = 0.0279 m and
 * PROBE_W = 0.084 m, a column catching three cells read -0.003 m and
 * the one next to it catching two read -0.336 m, on the same flat
 * still-water surface. The ratio above cancels that completely, and
 * the reading no longer depends on PROBE_W at all.
 *
 * This is an integral measure, so it is smooth and does not depend on
 * finding the alpha = 0.5 cell. Its resolution is the column width:
 * use PROBE_W about 2 to 4 cell widths.
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
#define PROBE_X     {5.0, 6.0, 10.0, 14.0}   /* m, station positions   */
#define PROBE_W     0.10                    /* m, column width        */
#define Y_BED      -1.0                      /* m, bed level           */
#define Y_TOP       0.5                      /* m, top of the domain   */
#define Y_SWL       0.0                      /* m, still water level   */
#define WATER_PHASE 1
#define LOG_FILE   "wave_probes.csv"
/* --------------------------------------------------------------------- */

/* A static "have I written the header yet" flag survives re-initialising
 * the case, so a second run appends to the first run's data and the file
 * ends up with two time histories interleaved. udf_restarted watches the
 * clock instead: fresh on the first call, and again whenever the time
 * goes backwards. */
static real last_t = -1.0;

DEFINE_EXECUTE_AT_END(wave_probes)
{
#if !RP_HOST
    Domain *mix = Get_Domain(1);
    Thread *t, *tw;
    cell_t c;
    real xc[ND_ND];
    real px[N_PROBE] = PROBE_X;
    real wvol[N_PROBE];   /* water volume in the column */
    real tvol[N_PROBE];   /* total volume in the column */
    real eta[N_PROBE];
    int i;

    for (i = 0; i < N_PROBE; i++) { wvol[i] = 0.0; tvol[i] = 0.0; }

    thread_loop_c(t, mix)
    {
        tw = THREAD_SUB_THREAD(t, WATER_PHASE);
        begin_c_loop_int(c, t)
        {
            C_CENTROID(xc, c, t);
            for (i = 0; i < N_PROBE; i++)
            {
                if (fabs(xc[0] - px[i]) <= 0.5 * PROBE_W)
                {
                    real v = C_VOLUME(c, t);
                    wvol[i] += C_VOF(c, tw) * v;
                    tvol[i] += v;
                }
            }
        }
        end_c_loop_int(c, t)
    }

    for (i = 0; i < N_PROBE; i++)
    {
        wvol[i] = PRF_GRSUM1(wvol[i]);
        tvol[i] = PRF_GRSUM1(tvol[i]);
        eta[i] = (tvol[i] > 0.0)
                   ? Y_BED + (wvol[i] / tvol[i]) * (Y_TOP - Y_BED) - Y_SWL
                   : 0.0;
    }

    if (UDF_IS_WRITER)
    {
        FILE *fp;
        if (udf_restarted(CURRENT_TIME, &last_t))
        {
            fp = fopen(LOG_FILE, "w");
            if (fp)
            {
                fprintf(fp, "t");
                for (i = 0; i < N_PROBE; i++) fprintf(fp, ",eta_x%g", px[i]);
                fprintf(fp, "\n");
                fclose(fp);
            }
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
