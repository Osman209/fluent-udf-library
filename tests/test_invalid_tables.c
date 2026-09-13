#include <assert.h>
#include "udf_common.h"
int main(void){const char *bad[]={"time,value\n","0,10\n1,broken\n","0,10\n1,nan\n","0,10\n1,inf\n","0,10\n0,20\n","1,10\n0,20\n","0,,10\n","0 10 20\n","0 1e999\n","0 10\n1\n","0,10,\n"};unsigned i;double *p=NULL;
 for(i=0;i<sizeof(bad)/sizeof(bad[0]);i++){FILE *f=fopen("invalid.txt","w");assert(f);fputs(bad[i],f);fclose(f);assert(udf_read_table("invalid.txt",2,&p)==-2);assert(p==NULL);}
 assert(udf_read_table("invalid.txt",65,&p)==-2);
 {FILE *f=fopen("invalid.txt","w");assert(f);fputs("# header\n0, 10 # value\n1\t20\n",f);fclose(f);assert(udf_read_table("invalid.txt",2,&p)==2);assert(p[3]==20);free(p);}
 remove("invalid.txt");puts("PASS strict table input validation");return 0;}
