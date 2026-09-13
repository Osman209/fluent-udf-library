/* Execute the shipped DEFINE_SPECIFIC_HEAT, not a copied formula. */
#include <assert.h>
#include "../profiles/property_from_table.c"
static void near(double a,double b){assert(fabs(a-b)<0.05);}
int main(void){
    real h,hp,hm,cp; FILE *fp=fopen("cp_T.txt","w");
    assert(fp);fprintf(fp,"300 1000\n400 2000\n");fclose(fp);
    cp=cp_T(325,300,&h,NULL);near(cp,1250);near(h,28125);
    cp_T(325.125,300,&hp,NULL);cp_T(324.875,300,&hm,NULL);near((hp-hm)/.25,cp);
    cp=cp_T(450,300,&h,NULL);near(cp,2000);near(h,250000);
    cp=cp_T(250,300,&h,NULL);near(cp,1000);near(h,-50000);
    cp_T(325,350,&h,NULL);near(h,-34375);
    cp_T(350,325,&h,NULL);near(h,34375);
    cp_T(325,325,&h,NULL);near(h,0);
    free(tb_cp.x);free(tb_cp.y);free(tb_cp.h);
    tb_cp.x=tb_cp.y=tb_cp.h=NULL;tb_cp.n=tb_cp.loaded=0;
    fp=fopen("cp_T.txt","w");assert(fp);fprintf(fp,"300 1500\n");fclose(fp);
    cp=cp_T(450,250,&h,NULL);near(cp,1500);near(h,300000);
    free(tb_cp.x);free(tb_cp.y);free(tb_cp.h);remove("cp_T.txt");
    puts("PASS actual cp: off-knot integral, derivative, endpoints, reverse reference, one row");return 0;
}
