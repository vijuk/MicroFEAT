#include "stdafx.h"
#include <stdio.h>

double KK[10][10],M0[10][10],U[10],U_t[10],U_tt[10],u[10],u_t[10],u_tt[10],_u[10];
//P[10][10]

main()
{
int Simple_Solve();

double D /*>=*/ =0.50;   /* delta */
double a /*>=*/ =0.25*sqr(0.5+D);  /* alpha */
double dt=.28;
double th=1.40;

double a0,a1,a2,a3,a4,a5,a6,a7,a8;

int DirInt = WILSONT;
int i,j,k;

/* Compute effective stiffness matrix  = K + a0 * M */

for (i=0;i<no_dof;i++)
	{
	u[i]=u_t[i]=u_tt[i]=0.0;
	U[i]=U_t[i]=U_tt[i]=0.0;
	_u[i]=0.0;
	}

U_tt[1]=u_tt[1]=10;

switch (DirInt) {
case NEWMARK :
	{
		 a0=1.0/(a*sqr(dt));
		 a1=D/(a*dt);
		 a2=1.0/(a*dt);
		 a3=1.0/(2*a)-1.0;
		 a4=D/a-1.0;
		 a5=dt/2.0*(D/a-2.0);
		 a6=dt*(1-D);
		 a7=dt*D;
		 break;
	}

case WILSONT :
	{
		 a0=6.0/sqr(th*dt);
		 a1=3.0/(th*dt);
		 a2=2.0*a1;
		 a3=th*dt/2.0;
		 a4=a0/th;
		 a5=-a2/th;
		 a6=1.0-3.0/th;
		 a7=dt/2.0;
		 a8=sqr(dt)/6.0;
		 break;
	}

case CENTRAL :
	{
		a0=1.0/sqr(dt);
		a1=1.0/(2.0*dt);
		a2=2*a0;
		a3=1.0/a2;
		a4=0;
		a5=0;
		a6=0;
		a7=0;
		for (i=0;i<no_dof;i++)
			_u[i]=U[i]-dt*U_t[i]+a3*U_tt[i];
		break;
	}

}

for (i=0;i<no_dof;i++)
for (j=0;j<no_dof;j++)
		{
		if ((DirInt==NEWMARK) || (DirInt==WILSONT))	K[i][j]+=a0*M[i][j];
		if (DirInt==CENTRAL)	{
					M0[i][j]=M[i][j];
					M[i][j]=a0*M[i][j];
					}
		}

printf("\nK\n");
for (i=0;i<no_dof;i++)
	{
	printf("\n");
	for (j=0;j<no_dof;j++)
		printf("%10.6lf ",K[i][j]);
	}

printf("\nM\n");
for (i=0;i<no_dof;i++)
	{
	printf("\n");
	for (j=0;j<no_dof;j++)
		printf("%10.6lf ",M[i][j]);
	}

/****************************** LOOOOOP */

for (k=0;k<10;k++) {

for (i=0;i<no_dof;i++)
		{
		U[i]=u[i];
		U_t[i]=u_t[i];
		U_tt[i]=u_tt[i];
		}

for (i=0;i<no_dof;i++)
	{
	AE[i]=P[i][k];      /* Input */
	for (j=0;j<no_dof;j++)
		{
		if (DirInt==NEWMARK)
			AE[i]=AE[i]+M[i][j]*(a0*U[j]+a2*U_t[j]+a3*U_tt[j]);
		if (DirInt==WILSONT)
			AE[i]=AE[i]+th*(P[i][k+1]-P[i][k])+M[i][j]*(a0*U[j]+a2*U_t[j]+2.0*U_tt[j]);
		if (DirInt==CENTRAL)
			AE[i]=AE[i]-(K[i][j]-a2*M0[i][j])*U[j]-(M[i][j])*_u[j];
		}
	}

printf("\n AE \n");
for (i=0;i<no_dof;i++)  printf("%10.6lf ",AE[i]);

for (i=0;i<no_dof;i++)
	for (j=0;j<no_dof;j++)
	   {
	   if ((DirInt==NEWMARK) || (DirInt==WILSONT)) KK[i][j]=K[i][j];
	   if (DirInt==CENTRAL)
			{
			KK[i][j]=K[i][j];
			K[i][j]=M[i][j];
			M[i][j]=KK[i][j];
 			}
	   }

if (DirInt==CENTRAL)
	for (i=0;i<no_dof;i++)
	for (j=0;j<no_dof;j++)   KK[i][j]=K[i][j];

Simple_Solve();

if (DirInt==CENTRAL)
	for (i=0;i<no_dof;i++)
	for (j=0;j<no_dof;j++)    K[i][j]=KK[i][j];

for (i=0;i<no_dof;i++)
	for (j=0;j<no_dof;j++)
	   {
	   if ((DirInt==NEWMARK) || (DirInt==WILSONT)) K[i][j]=KK[i][j];
	   if (DirInt==CENTRAL)
			{
			KK[i][j]=K[i][j];
			K[i][j]=M[i][j];
			M[i][j]=KK[i][j];
 			}
	   }

printf("\n Time  = %ddt   \n",k);
/*
for (i=0;i<no_dof;i++)
	   printf("%10.6lf ",U[i]);
*/

/* Calculate accelerations and displacements at (t+dt) */

for (i=0;i<no_dof;i++)
      {
	u[i]=AE[i];

if (DirInt==NEWMARK)
	{
	u_tt[i]=a0*(u[i]-U[i])-a2*U_t[i]-a3*U_tt[i];
	u_t[i]=U_t[i]+a6*U_tt[i]+a7*u_tt[i];
	}

if (DirInt==WILSONT)
	{
	u_tt[i]=a4*(u[i]-U[i])+a5*U_t[i]+a6*U_tt[i];
	u_t[i]=U_t[i]+a7*(u_tt[i]+U_tt[i]);
	u[i]=U[i]+dt*U_t[i]+a8*(u_tt[i]+2.0*U_tt[i]);
	}

if (DirInt==CENTRAL)
	{
	u_tt[i]=a0*(_u[i]-2*U[i]+u[i]);
	u_t[i]=a1*(-_u[i]+u[i]);
	_u[i]=U[i];
	}
printf("\n U[%d] = %10.6lf  U_t[%d] = %10.6lf  U_tt[%d] = %10.6lf \n",i,u[i],i,u_t[i],i,u_tt[i]);
      }
getchar();
}
}

