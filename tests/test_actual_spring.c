#include <assert.h>
#include "../motion/sdof_spring_damper.c"
int main(void){Dynamic_Thread dt={0};real prop[SDOF_N_PROPS]={0};
 dt.cg[1]=.2;dt.vel[1]=.3;
 spring_damper(prop,&dt,10,.01);assert(fabs(prop[SDOF_LOAD_F_Y]+100.6)<1e-4);
 dt.cg[1]=.1;spring_damper(prop,&dt,11,.01);assert(fabs(prop[SDOF_LOAD_F_Y]+50.6)<1e-4);
 run_last_t=-1; /* simulate a library reload at a displaced restart */
 spring_damper(prop,&dt,11,.01);assert(fabs(prop[SDOF_LOAD_F_Y]+50.6)<1e-4);
 spring_damper(prop,&dt,0,.01);assert(fabs(prop[SDOF_LOAD_F_Y]+50.6)<1e-4);
 remove(LOG_FILE);puts("PASS fixed spring datum through time reset and reload");return 0;}
