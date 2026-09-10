/*
 * wave_init.c
 *
 * Still-water initialisation for a VOF tank: water fraction 1 below
 * Y_SWL, 0 above, velocities zero. Hook: User-Defined > Function Hooks
 * > Initialization, then run Initialize (standard) first and this UDF
 * fills in the phase field.
 *
 * Why a UDF rather than Patch: patching needs a cell register aligned
 * with the water level, and this is one click less to get wrong. It
 * also sets a smoothed interface over one cell height, which starts
 * the run with a cleaner surface than a hard patch.
 *
 * Pressure: this UDF does NOT set the hydrostatic pressure. With
 * gravity on and the operating density set to the air density (the
 * standard VOF setting), Fluent recovers the hydrostatic field within
 * the first few iterations of the first time step. If you want a
 * hydrostatic start too, use the Open Channel initialisation instead,
 * or add C_P(c,t) = rho_w * g * (Y_SWL - y) below the surface after
 * checking your operating pressure reference.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define Y_SWL       0.0
#define WATER_PHASE 1     /* phase index of water (primary = 0)         */
#define SMOOTH_DY   0.01  /* m, interface smoothing height              */
/* --------------------------------------------------------------------- */

DEFINE_INIT(still_water, mix)
{
#if !RP_HOST
    Thread *t, *tw;
    cell_t c;
    real xc[ND_ND], a;

    thread_loop_c(t, mix)
    {
        tw = THREAD_SUB_THREAD(t, WATER_PHASE);
        begin_c_loop_all(c, t)
        {
            C_CENTROID(xc, c, t);
            a = 0.5 + (Y_SWL - xc[1]) / (SMOOTH_DY > 0.0 ? SMOOTH_DY : 1.0e-9);
            if (a > 1.0) a = 1.0;
            if (a < 0.0) a = 0.0;
            C_VOF(c, tw) = a;
            C_U(c, t) = 0.0;
            C_V(c, t) = 0.0;
#if ND_ND == 3
            C_W(c, t) = 0.0;
#endif
        }
        end_c_loop_all(c, t)
    }
#endif
}
