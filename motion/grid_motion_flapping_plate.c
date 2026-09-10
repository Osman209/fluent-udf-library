/*
 * grid_motion_flapping_plate.c
 *
 * Deforming boundary: a cantilever plate clamped at x = X_ROOT whose
 * transverse displacement follows a prescribed mode shape
 *
 *     y(x, t) = AMP * s(x)^2 * sin(2 pi f t),   s = (x - X_ROOT) / LENGTH
 *
 * (quadratic mode shape, zero at the root, maximum at the tip).
 * Hook: Dynamic Mesh Zones > User-Defined > Mesh Motion UDF, on the
 * plate wall zone. Use smoothing (and remeshing if needed) on the
 * fluid zone around it.
 *
 * Parallel-safe pattern (the one that fails most often on the forums):
 *   - the whole node loop is inside #if !RP_HOST
 *   - SET_DEFORMING_THREAD_FLAG is called on the adjacent cell thread
 *   - each node is moved only once per step via NODE_POS_NEED_UPDATE /
 *     NODE_POS_UPDATED, which matters because a node is shared by
 *     several faces and, in parallel, may be seen on several partitions
 *
 * The UDF moves nodes by the INCREMENT between the previous and the
 * current time, so the undeformed node positions do not need to be
 * stored. This also makes it restart-safe as long as the restart time
 * is the time the case was saved at.
 */

#include "udf.h"

/* ---- user parameters ------------------------------------------------ */
#define X_ROOT   0.0     /* m, clamped end                              */
#define LENGTH   1.0     /* m, plate length along x                     */
#define AMP      0.05    /* m, tip displacement amplitude               */
#define FREQ     2.0     /* Hz                                          */
#define DISP_DIR 1       /* 0 = x, 1 = y, 2 = z; direction of displacement */
#define SPAN_DIR 0       /* axis along the plate length                 */
/* --------------------------------------------------------------------- */

static real shape(real x)
{
    real s = (x - X_ROOT) / LENGTH;
    if (s < 0.0) s = 0.0;
    if (s > 1.0) s = 1.0;
    return s * s;
}

DEFINE_GRID_MOTION(flapping_plate, domain, dt, time, dtime)
{
#if !RP_HOST
    Thread *tf = DT_THREAD(dt);
    face_t f;
    Node *v;
    int n;
    real w  = 2.0 * M_PI * FREQ;
    real dq = sin(w * time) - sin(w * (time - dtime));

    SET_DEFORMING_THREAD_FLAG(THREAD_T0(tf));

    begin_f_loop(f, tf)
    {
        f_node_loop(f, tf, n)
        {
            v = F_NODE(f, tf, n);
            if (NODE_POS_NEED_UPDATE(v))
            {
                NODE_POS_UPDATED(v);
                NODE_COORD(v)[DISP_DIR] += AMP * shape(NODE_COORD(v)[SPAN_DIR]) * dq;
            }
        }
    }
    end_f_loop(f, tf)
#endif
}
