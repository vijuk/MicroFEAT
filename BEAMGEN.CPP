#include <stdio.h>
FILE *outf;
main()
{
int i,N;
N=100;
outf=fopen("DAT.DAT","w");
for (i=1;i<=N+1;i++)
	fprintf(outf,"%3d %6.3lf %6.3lf %6.3lf \n",i,(i-1)*(100.0/N),0,0);
for (i=1;i<=N;i++)
	fprintf(outf,"%3d 1 %3d %3d 0.24 0.0072 0.0072 0.0072 0 0\n",i,i,i+1);
fcloseall();
}
