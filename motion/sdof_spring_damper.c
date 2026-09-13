/*
 * sdof_spring_damper.c
 *
 * Six-DOF rigid body with linear springs and dampers on any DOF, plus
 * DOF constraints. This is the vortex-induced-vibration / elastically
 * mounted cylinder / bridge deck / mooring-spring setup.
 *
 * Hook: Dynamic Mesh Zones > Rigid Body > Six DOF > Six DOF UDF/Properties.
 * Enable "Six DOF" in the dynamic mesh model, and set the same zone
 * (and only that zone as an active body, see docs/pitfalls_6dof.md)
 * to Six DOF. Gravity: if gravity is switched on in Operating
 * Conditions, the Fluent 6DOF solver includes it. Do NOT add it here
 * as well.
 *
 * Forces:   F_i = -K_i (x_i - x0_i) - C_i v_i      (global frame)
 * Moments:  M_i = -KR_i (theta_i - theta0_i) - CR_i omega_i
 * x0, theta0 are explicit fixed reference parameters, preserved on restart.
 * The rotational spring is written for small angles about each axis.
 *
 * Output: sdof_motion.csv in the working directory, one line per time
 * step: t, x, y, z, vx, vy, vz, thx, thy, thz, wx, wy, wz.
 * Compare the free-decay response with validation/sdof_reference.py.
 *
 * Macro compatibility note (measured, Fluent 2025 R1): DT_OMEGA(dt) is
 * a scalar there and DT_OMEGA(dt)[i] does not compile. The array is
 * DT_OMEGA_CG(dt), which is what this file uses. DT_CG, DT_VEL_CG and
 * DT_THETA are arrays.
 */

#include "udf.h"
#include "dynamesh_tools.h"
#define UDF_INSIDE_FLUENT
#include "udf_common.h"

/* ---- user parameters ------------------------------------------------ */
#define MASS   10.0      /* kg */
#define IXX     1.0      /* kg m^2 */
#define IYY     1.0
#define IZZ     1.0

/* translational spring stiffness [N/m] and damping [N s/m], per axis  */
#define KX  0.0
#define KY  500.0
#define KZ  0.0
#define CX  0.0
#define CY  2.0
#define CZ  0.0

/* rotational spring [N m/rad] and damping [N m s/rad], per axis       */
#define KRX 0.0
#define KRY 0.0
#define KRZ 0.0
#define CRX 0.0
#define CRY 0.0
#define CRZ 0.0

/* constraints: 1 = lock that DOF                                        */
#define LOCK_TX 1
#define LOCK_TY 0
#define LOCK_TZ 1
#define LOCK_RX 1
#define LOCK_RY 1
#define LOCK_RZ 1

/* Global spring equilibrium position [m] and small-angle orientation [rad].
 * Set these to the physical equilibrium datum, NOT a displaced restart CG. */
#define SPRING_REF_X {0.0, 0.0, 0.0}
#define SPRING_REF_TH {0.0, 0.0, 0.0}
#define LOG_FILE "sdof_motion.csv"
/* --------------------------------------------------------------------- */

/* Fixed physical datum does not change on Initialize or library reload. */
static const real x0[3] = SPRING_REF_X;
static const real th0[3] = SPRING_REF_TH;
static real last_logged_time = -1.0;
static double run_last_t = -1.0;

DEFINE_SDOF_PROPERTIES(spring_damper, prop, dt, time, dtime)
{
    real *x  = DT_CG(dt);
    real *v  = DT_VEL_CG(dt);
    real *th = DT_THETA(dt);
    real *w  = DT_OMEGA_CG(dt);
    int i;

    if (udf_restarted(time, &run_last_t))
    {
        last_logged_time = -1.0;
        if (UDF_IS_WRITER)
        {
            FILE *fp = fopen(LOG_FILE, time > 0.0 ? "a+" : "w");
            if (fp)
            {
                fseek(fp, 0, SEEK_END);
                if (ftell(fp) == 0) fprintf(fp, "t,x,y,z,vx,vy,vz,thx,thy,thz,wx,wy,wz\n");
                fclose(fp);
            }
            Message("sdof_spring_damper: reference CG = (%g, %g, %g)\n",
                    x0[0], x0[1], x0[2]);
        }
    }

    prop[SDOF_MASS] = MASS;
    prop[SDOF_IXX]  = IXX;
    prop[SDOF_IYY]  = IYY;
    prop[SDOF_IZZ]  = IZZ;
    prop[SDOF_IXY]  = 0.0;
    prop[SDOF_IXZ]  = 0.0;
    prop[SDOF_IYZ]  = 0.0;

    /* loads are given in the global frame */
    prop[SDOF_LOAD_LOCAL] = FALSE;

    prop[SDOF_LOAD_F_X] = -KX * (x[0] - x0[0]) - CX * v[0];
    prop[SDOF_LOAD_F_Y] = -KY * (x[1] - x0[1]) - CY * v[1];
    prop[SDOF_LOAD_F_Z] = -KZ * (x[2] - x0[2]) - CZ * v[2];

    prop[SDOF_LOAD_M_X] = -KRX * (th[0] - th0[0]) - CRX * w[0];
    prop[SDOF_LOAD_M_Y] = -KRY * (th[1] - th0[1]) - CRY * w[1];
    prop[SDOF_LOAD_M_Z] = -KRZ * (th[2] - th0[2]) - CRZ * w[2];

    prop[SDOF_ZERO_TRANS_X] = LOCK_TX ? TRUE : FALSE;
    prop[SDOF_ZERO_TRANS_Y] = LOCK_TY ? TRUE : FALSE;
    prop[SDOF_ZERO_TRANS_Z] = LOCK_TZ ? TRUE : FALSE;
    prop[SDOF_ZERO_ROT_X]   = LOCK_RX ? TRUE : FALSE;
    prop[SDOF_ZERO_ROT_Y]   = LOCK_RY ? TRUE : FALSE;
    prop[SDOF_ZERO_ROT_Z]   = LOCK_RZ ? TRUE : FALSE;

    /* This UDF can be called more than once per time step (6DOF
     * sub-iterations). Log only when the flow time has advanced.      */
    if (UDF_IS_WRITER && time > last_logged_time)
    {
        FILE *fp = fopen(LOG_FILE, "a");
        if (fp)
        {
            fprintf(fp, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                    time, x[0], x[1], x[2], v[0], v[1], v[2],
                    th[0], th[1], th[2], w[0], w[1], w[2]);
            fclose(fp);
        }
        last_logged_time = time;
    }
}
