/*
 * sdof_speed_heading_controller.c
 *
 * Six-DOF body driven to a target forward speed and a target heading by
 * a simple controller: proportional thrust on surge, PD moment on yaw.
 * This is the "self-propelled hull turning at speed" case. The thrust
 * and yaw moment are applied in the BODY frame (SDOF_LOAD_LOCAL = TRUE).
 *
 * Axis convention used below (change the three macros to match yours):
 *   FWD  = body forward axis index at zero heading (0 = x)
 *   VERT = vertical axis index, yaw is about it            (1 = y)
 *   LAT  = the remaining horizontal axis                    (2 = z)
 * With Y up and X forward, a positive yaw about +Y turns the bow toward
 * -Z (right-hand rule). If your case is Z up, set FWD 0, VERT 2, LAT 1
 * and the forward vector becomes (cos psi, sin psi, 0).
 *
 * Staging: the thrust can be held back for T_HOLD seconds so the hull
 * settles to its hydrostatic draft first (float phase). Starting the
 * thrust and the turn on an unsettled hull is the usual reason a 6DOF
 * VOF hull run diverges in the first second. See docs/pitfalls_6dof.md.
 *
 * Output: controller_log.csv (t, speed, heading_deg, thrust, yaw_moment).
 */

#include "udf.h"
#include "dynamesh_tools.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define MASS   8400.0        /* kg, must be the real displacement mass  */
#define IXX    5000.0        /* kg m^2 */
#define IYY    40000.0
#define IZZ    40000.0

#define FWD  0
#define VERT 1
#define LAT  2

#define V_TARGET      5.0                 /* m/s forward speed target   */
#define PSI_TARGET   (15.0*M_PI/180.0)    /* rad, heading target        */
#define T_HOLD        2.0                 /* s, float phase, no thrust  */

#define K_SURGE     3000.0   /* N per (m/s) of speed error              */
#define THRUST_MAX 20000.0   /* N, saturation                            */
#define KP_YAW       125.0   /* N m per rad of heading error             */
#define KD_YAW       500.0   /* N m per (rad/s) of yaw rate              */
#define MOMENT_MAX  5000.0   /* N m, saturation                          */

#define LOCK_HEAVE 0   /* usually free for a floating hull               */
#define LOCK_ROLL  0
#define LOCK_PITCH 0
#define LOCK_SWAY  0

#define LOG_FILE "controller_log.csv"
/* --------------------------------------------------------------------- */

static real clamp(real v, real lo, real hi)
{
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

/* Both of these used to survive Initialize. The consequence was worse
 * than a duplicated log: last_logged_time kept the end time of the
 * previous run, so "log only when the clock has advanced" never passed
 * again and the second run recorded nothing at all. */
static real last_logged_time = -1.0;
static double log_last_t = -1.0;

DEFINE_SDOF_PROPERTIES(speed_heading, prop, dt, time, dtime)
{
    real *v   = DT_VEL_CG(dt);
    real *th  = DT_THETA(dt);
    real *w   = DT_OMEGA_CG(dt);
    real psi, r, fwd_dir[3], u_fwd, thrust, moment, err;

    prop[SDOF_MASS] = MASS;
    prop[SDOF_IXX]  = IXX;
    prop[SDOF_IYY]  = IYY;
    prop[SDOF_IZZ]  = IZZ;
    prop[SDOF_IXY]  = 0.0;
    prop[SDOF_IXZ]  = 0.0;
    prop[SDOF_IYZ]  = 0.0;

    prop[SDOF_LOAD_LOCAL] = TRUE;

    /* heading and yaw rate about the vertical axis */
    psi = th[VERT];
    r   = w[VERT];

    /* body forward direction in global coordinates (Y-up convention) */
    NV_S(fwd_dir, =, 0.0);
    fwd_dir[FWD] =  cos(psi);
    fwd_dir[LAT] = (VERT == 1) ? -sin(psi) : sin(psi);
    u_fwd = v[0]*fwd_dir[0] + v[1]*fwd_dir[1] + v[2]*fwd_dir[2];

    if (time < T_HOLD)
    {
        thrust = 0.0;
        moment = 0.0;
    }
    else
    {
        err    = V_TARGET - u_fwd;
        thrust = clamp(K_SURGE * err, -THRUST_MAX, THRUST_MAX);

        err    = PSI_TARGET - psi;
        moment = clamp(KP_YAW * err - KD_YAW * r, -MOMENT_MAX, MOMENT_MAX);
    }

    /* loads in the body frame: thrust along body forward axis,
       moment about the body vertical axis */
    prop[SDOF_LOAD_F_X] = 0.0;
    prop[SDOF_LOAD_F_Y] = 0.0;
    prop[SDOF_LOAD_F_Z] = 0.0;
    prop[SDOF_LOAD_M_X] = 0.0;
    prop[SDOF_LOAD_M_Y] = 0.0;
    prop[SDOF_LOAD_M_Z] = 0.0;

    if (FWD == 0) prop[SDOF_LOAD_F_X] = thrust;
    if (FWD == 1) prop[SDOF_LOAD_F_Y] = thrust;
    if (FWD == 2) prop[SDOF_LOAD_F_Z] = thrust;

    if (VERT == 0) prop[SDOF_LOAD_M_X] = moment;
    if (VERT == 1) prop[SDOF_LOAD_M_Y] = moment;
    if (VERT == 2) prop[SDOF_LOAD_M_Z] = moment;

    prop[SDOF_ZERO_TRANS_X] = FALSE;
    prop[SDOF_ZERO_TRANS_Y] = FALSE;
    prop[SDOF_ZERO_TRANS_Z] = FALSE;
    prop[SDOF_ZERO_ROT_X]   = FALSE;
    prop[SDOF_ZERO_ROT_Y]   = FALSE;
    prop[SDOF_ZERO_ROT_Z]   = FALSE;

    if (VERT == 1)
    {
        if (LOCK_HEAVE) prop[SDOF_ZERO_TRANS_Y] = TRUE;
        if (LOCK_SWAY)  prop[SDOF_ZERO_TRANS_Z] = TRUE;
        if (LOCK_ROLL)  prop[SDOF_ZERO_ROT_X]   = TRUE;
        if (LOCK_PITCH) prop[SDOF_ZERO_ROT_Z]   = TRUE;
    }
    else if (VERT == 2)
    {
        if (LOCK_HEAVE) prop[SDOF_ZERO_TRANS_Z] = TRUE;
        if (LOCK_SWAY)  prop[SDOF_ZERO_TRANS_Y] = TRUE;
        if (LOCK_ROLL)  prop[SDOF_ZERO_ROT_X]   = TRUE;
        if (LOCK_PITCH) prop[SDOF_ZERO_ROT_Y]   = TRUE;
    }

    if (UDF_IS_WRITER)
    {
        FILE *fp;
        if (udf_restarted(time, &log_last_t))
        {
            fp = fopen(LOG_FILE, "w");
            if (fp) { fprintf(fp, "t,u_fwd,heading_deg,thrust,yaw_moment\n"); fclose(fp); }
            last_logged_time = -1.0;   /* so the new run logs from its first step */
        }
        if (time > last_logged_time)
        {
            fp = fopen(LOG_FILE, "a");
            if (fp)
            {
                fprintf(fp, "%g,%g,%g,%g,%g\n",
                        time, u_fwd, psi * 180.0 / M_PI, thrust, moment);
                fclose(fp);
            }
            last_logged_time = time;
        }
    }
}
