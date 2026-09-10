/*
 * MOCK udf.h  --  NOT the Ansys header.
 *
 * Minimal stand-in so every UDF in this library can be syntax-checked
 * with plain gcc (tests/run_tests.sh). It defines the macro NAMES the
 * UDFs use with plausible types. It does not reproduce Fluent
 * behaviour, and a UDF that passes here can still fail to compile in
 * Fluent if a macro name or signature differs in your Fluent version.
 * The real test is the Fluent build; this mock only catches typos and
 * C errors before you get there.
 */
#ifndef MOCK_UDF_H
#define MOCK_UDF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef ND_ND
#define ND_ND 3
#endif
#define TRUE  1
#define FALSE 0
#define RP_HOST 0
#define RP_NODE 0
#define I_AM_NODE_ZERO_P 1
#define PRF_GRSUM1(x) (x)

typedef double real;
typedef int face_t;
typedef int cell_t;
typedef struct Thread_s { int id; struct Thread_s *t0; } Thread;
typedef struct Domain_s { int id; Thread *threads; } Domain;
typedef struct Node_s { real x[3]; int flag; } Node;
typedef struct DynamicThread_s {
    real cg[3], vel[3], theta[3], omega[3]; Thread *thread;
} Dynamic_Thread;

#define Message printf
#define CURRENT_TIME (0.0)

#define NV_S(a, op, s)  do { (a)[0] op (s); (a)[1] op (s); (a)[2] op (s); } while (0)
#define NV_VEC(a) real a[ND_ND]

#define DEFINE_CG_MOTION(name, dt, vel, omega, time, dtime) \
    void name(Dynamic_Thread *dt, real vel[], real omega[], real time, real dtime)
#define DEFINE_SDOF_PROPERTIES(name, prop, dt, time, dtime) \
    void name(real *prop, Dynamic_Thread *dt, real time, real dtime)
#define DEFINE_PROFILE(name, t, i) void name(Thread *t, int i)
#define DEFINE_SOURCE(name, c, t, dS, eqn) real name(cell_t c, Thread *t, real dS[], int eqn)
#define DEFINE_INIT(name, d) void name(Domain *d)
#define DEFINE_EXECUTE_AT_END(name) void name(void)
#define DEFINE_ON_DEMAND(name) void name(void)
#define DEFINE_PROPERTY(name, c, t) real name(cell_t c, Thread *t)
#define DEFINE_SPECIFIC_HEAT(name, T, Tref, h, yi) \
    real name(real T, real Tref, real *h, real *yi)
#define C_T(c, t) (300.0)
#define RP_Get_Real(s) (101325.0)
#define DEFINE_GRID_MOTION(name, d, dt, time, dtime) \
    void name(Domain *d, Dynamic_Thread *dt, real time, real dtime)

/* face / cell loops */
#define begin_f_loop(f, t)  { int f##_n = 0; for (f = 0; f < f##_n; f++) {
#define end_f_loop(f, t)    } }
#define begin_c_loop_int(c, t) { int c##_n = 0; for (c = 0; c < c##_n; c++) {
#define end_c_loop_int(c, t)   } }
#define begin_c_loop_all(c, t) { int c##_n = 0; for (c = 0; c < c##_n; c++) {
#define end_c_loop_all(c, t)   } }
#define thread_loop_c(t, d) for (t = (d)->threads; t != NULL; t = NULL)
#define f_node_loop(f, t, n) for (n = 0; n < 0; n++)

static real mock_vec3[3];
static real mock_scalar;
static Node mock_node;

#define F_CENTROID(x, f, t) do { (x)[0] = 0; (x)[1] = 0; if (ND_ND == 3) (x)[2] = 0; } while (0)
#define F_AREA(A, f, t)     F_CENTROID(A, f, t)
#define C_CENTROID(x, c, t) F_CENTROID(x, c, t)
#define F_PROFILE(f, t, i)  mock_scalar
#define F_P(f, t)           (0.0)
#define F_STORAGE_R_N3V(f, t, sv) mock_vec3
#define SV_WALL_SHEAR 0
#define PRINCIPAL_FACE_P(f, t) 1
#define Get_Domain(i)       ((Domain *)NULL)
#define Lookup_Thread(d, id) ((Thread *)NULL)
#define THREAD_SUB_THREAD(t, i) (t)
#define THREAD_T0(t) (t)
#define C_VOF(c, t)     mock_scalar
#define C_VOLUME(c, t)  (1.0)
#define C_R(c, t)       (1000.0)
#define C_U(c, t)       mock_scalar
#define C_V(c, t)       mock_scalar
#define C_W(c, t)       mock_scalar
#define C_P(c, t)       mock_scalar
#define C_MU_EFF(c, t)  (1.0e-3)
#define C_U_G(c, t)     mock_vec3
#define C_V_G(c, t)     mock_vec3
#define C_W_G(c, t)     mock_vec3
#define F_C0(f, t)      (0)
#define NV_MAG(a)       sqrt((a)[0]*(a)[0] + (a)[1]*(a)[1] + (a)[2]*(a)[2])

/* dynamic mesh */
#define DT_CG(dt)       ((dt)->cg)
#define DT_VEL_CG(dt)   ((dt)->vel)
#define DT_THETA(dt)    ((dt)->theta)
#define DT_OMEGA_CG(dt) ((dt)->omega)
#define DT_THREAD(dt)   ((dt)->thread)
#define SET_DEFORMING_THREAD_FLAG(t) ((void)0)
#define F_NODE(f, t, n) (&mock_node)
#define NODE_POS_NEED_UPDATE(v) ((v)->flag == 0)
#define NODE_POS_UPDATED(v) ((v)->flag = 1)
#define NODE_COORD(v) ((v)->x)

/* SDOF property indices (values arbitrary in the mock) */
enum {
    SDOF_MASS, SDOF_IXX, SDOF_IYY, SDOF_IZZ, SDOF_IXY, SDOF_IXZ, SDOF_IYZ,
    SDOF_LOAD_LOCAL, SDOF_LOAD_F_X, SDOF_LOAD_F_Y, SDOF_LOAD_F_Z,
    SDOF_LOAD_M_X, SDOF_LOAD_M_Y, SDOF_LOAD_M_Z,
    SDOF_ZERO_TRANS_X, SDOF_ZERO_TRANS_Y, SDOF_ZERO_TRANS_Z,
    SDOF_ZERO_ROT_X, SDOF_ZERO_ROT_Y, SDOF_ZERO_ROT_Z, SDOF_N_PROPS
};

#endif
