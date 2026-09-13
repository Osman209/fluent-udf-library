#include <assert.h>
#include "mock/udf.h"
#undef F_AREA
#if ND_ND == 3
#define F_AREA(A,f,t) do{(A)[0]=0;(A)[1]=-1;(A)[2]=0;}while(0)
#else
#define F_AREA(A,f,t) do{(A)[0]=0;(A)[1]=-1;}while(0)
#endif
static real gu[3]={0,2,0},gz[3]={0,0,0};
#undef C_U_G
#undef C_V_G
#undef C_W_G
#define C_U_G(c,t) gu
#define C_V_G(c,t) gz
#define C_W_G(c,t) gz
#include "../motion/force_moment_logger.c"
int main(void){Thread t={0};real fv[3];visc_force_gradient(0,&t,fv);assert(fabs(fv[0]-.002)<1e-8);assert(fabs(fv[1])<1e-8);gu[1]=-2;visc_force_gradient(0,&t,fv);assert(fabs(fv[0]+.002)<1e-8);puts("PASS actual viscous-force direction and reversed gradient");return 0;}
