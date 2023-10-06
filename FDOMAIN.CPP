#include "stdafx.h"
#include "space.h"

void Fill_Upper(double **A,int n)
{
for (int i=0;i<n-1;i++)
	for (int j=i+1;j<n;j++)
		A[i][j]=A[j][i];   
}


void matout(double **a,int m,int n)
{
	int i,j;
	for (i=0;i<m;i++)
	{
		fprintf(outfile,"\n");
		for (j=0;j<n;j++)
		{
			fprintf(outfile,"%10.4lf",a[i][j]);
		}
	}
	fprintf(outfile,"\n");
}

/*
matinp(a,m,n)
float a[10][10];
int m,n;
{
	int i,j;
	if ((m>10) || (n>10)) return;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			scanf("%f",&a[i][j]);
}

store(b,a,m,n)
float a[10][10],b[10][10];
int m,n;
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			b[i][j]=a[i][j];
}

*/



/****   Functions used for dynamics starts here   ****/

/*
main()
{
	int i,j,n;
	float a[10][10],b[10][10],c[10][10];
	printf("\nSize of the matrix  ? ");
	scanf("%d",&n);
	printf("%d\n",n);
	printf("Input matrix : \n",n,n);
	matinp(a,n,n);
	printf("\nInput > \n");
	matout(a,n,n);

	for (j=0;j<no_dof;j++)
	for (k=0;k<no_dof;k++)
		U[j][k]=SS[j][k];
	
	IDCMP=0;

	decomp(U,no_dof,3) 
	Decompose stiffness matrix into the factored form U'U *

	for (i=0;i<no_dof;i++)
	for (j=0;j<no_dof;j++)
		{
		sum=A[i][j];
		K1=I-1;
		for (k=0;k<K1;k++)
			sum-=A[k][i]*A[k][j];
		if (i==j)
			if (sum<=0) return 1
		else
			A[i][i]=sqrt(sum);
		}
		else
		A[i][j]=sum/A[i][j];
	}

printf(" Stiffness matrix decomposed.\n"); 

	for (j=0;j<no_dof;j++)
	for (k=0;k<no_dof;k++)
		CMU[j][k]=CMS[j][k];

	/* decompose mass matrix and copy stiffness matrix to CMU

	for (j=0;j<no_dof;j++)
	for (k=0;k<no_dof;k++)
		CMU[j][k]=SS[j][k];

	for (j=0;j<no_dof;j++)
	for (k=0;k<no_dof;k++)
		U[j][k]=CMS[j][k];
	
	IDCMP=0;
	decomp(U,no_dof,6);

	jacobi(a,b,n);

	printf("\nOutput > EIGEN SOLUTION :\n\n");

	printf("EIGEN VALUES >\n");
	matout(a,n,n);
	printf("EIGEN VECTOR >\n");
	matout(b,n,n);

	for (i=0;i<n;i++)
	{
		printf(" Eigen vector for %8.4f :\n\n",a[i][i]);
		for (j=0;j<n;j++)
			printf("                           %10.4f\n",b[i][j]);
		printf("\n");
	}
}
*/

/*
Cstructure::standard_jacobi()
{
a=K
b=M
n=no_dof;

	int i,j,k,m;
	double th,cf;
	int flag,ix,jx;
	int fg=0;
	double temp[10][10],bt[10][10],temp2[10][10],temp3[10][10];
	int sweep=1;
	flag=1;
	ix=0;
	jx=0;
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
		{
			if (i==j) temp2[i][j]=1; 
			else temp2[i][j]=0.0;
		}

	while (flag)
	{
		for (i=0;i<n;i++)
			for (j=0;j<n;j++)
				if (i==j) b[i][j]=1; 
				else b[i][j]=0.0;

		jx++;
		if (jx>n-1)
		{
			ix++;
			jx=ix+1;
			if (ix>n-2) {
				ix=0;
				jx=1;
				sweep++;
			}
		}

		cf=sqrt((a[ix][jx]*a[ix][jx]/(a[ix][ix]*a[jx][jx])));
		if (cf>pow(10.0,-2.0*sweep))
		{
			store(temp3,a,n,n);
			if (a[ix][ix]==a[jx][jx]) th=PI/4;
			else 
		             th=atan(2*a[ix][jx]/(a[ix][ix]-a[jx][jx]))/ 2.0;

			b[ix][ix]=cos(th);
			b[ix][jx]=-sin(th);
			b[jx][ix]=sin(th);
			b[jx][jx]=cos(th);

			mattran(b,bt,n,n);
			matmul(bt,a,temp,n,n,n);
			matmul(temp,b,a,n,n,n);
			matmul(temp2,b,temp,n,n,n);

			store(temp2,temp,n,n);
			flag=0;
			fg++;

			for (i=0;i<n;i++)
			{
			if (fabs((double)(temp3[i][i]-a[i][i])/temp3[i][i])
			        > 1e-12)
				{
					flag=1;
					break;
				}
				for (j=i+1;j<n;j++)
					if ((fabs(a[i][j])>1e-12))
					{
						flag=1;
						if (fg>100)
						{
							store(b,temp,n,n);
printf("Avoiding infinite loop > No. of sweeps           = %d \n", sweep);
printf("No. of transformations  = %d \n\n", fg);
							return;
						}
						break;
					}
			}

			if (!flag)
			{
				store(b,temp,n,n);
				printf("No. of sweeps          = %d \n", sweep);
				printf("No. of transformations = %d \n\n", fg);
				return;
			}
		}
	}
}

/*
mattran(a,b,m,n)
float a[10][10],b[10][10];
int m,n;
{
	int i,j,t;
	char c;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			b[j][i]=a[i][j];
}

matmul(a,b,c,m,p,n)
float a[10][10],b[10][10],c[10][10];
int m,p,n;
{
	int i,j,k;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
		{
			c[i][j]=0.0;
			for (k=0;k<p;k++)
				c[i][j]=c[i][j]+a[i][k]*b[k][j];
		}
}
*/

/* jacobi */

// 4. Page 525

/*
void Cstructure::Pre_Eigen()
{
// Copy Stiffness matrix K to U
for (j=0;j<no_dof;j++)
	for (k=0;k<no_dof;k++)
		U[j][k]=K[j][k];

idcmp="K";
// decompose stiffness matrix into factored form U'U
if (decompose(U,no_dof,3))
	{
	// Copy mass matrix to CMU
	for (j=0;j<no_dof;j++)
		for (k=0;k<no_dof;k++)
			CMU[j][k]=M[j][k];
	}
else
	{
	// Copy Stiffness matrix K to U
	for (j=0;j<no_dof;j++)
		for (k=0;k<no_dof;k++)	
					U[j][k]=M[j][k];

	idcmp="M";
	decompose(U,no_dof,6); // decompose mass matrix into factored form U'U

	// Copy mass matrix to CMU
	for (j=0;j<no_dof;j++)
		for (k=0;k<no_dof;k++)
			CMU[j][k]=K[j][k];
	}
}


// Obtain inverse transpose of U
Invert_U(U,no_dof);
standardise();

}


// Inverse transpose of U
void Cstructure::invert_U()      
{
for (i=0;i<no_dof;i++)
	{
	a[i][i]=1.0/a[i][i];
	n1=no_dof-1;
	for (i=0;i<n1;i++)
			{
			ip1=i+1;
			for (j=ip1-1;j<no_dof;j++)
					{
					sum=0;
					jm1=j-1;
					for (k=i-1;k<jm1;k++)
						sum-=a[k][i]*a[k][j];
					a[j][i]=sum*a[j][j];
					}
			}
	}
}


//d. Transform to standard symmetric form
void Cstructure::standardise()
{

//  Premultiply by (1/U)'
for (int ii=0;ii<no_dof;ii++)
	{
	int i=no_dof-ii+1;
	for (int ji=0;ji<no_dof;ji++)
		{
		int j=no_dof-ji+1;
		temp=0;
		for (int k=0;k<j;k++) temp+=U[j][k]*CMU[i][k];
		CMU[i][j]=temp;
		}
	}

// Post-multiply by (1/U)'
for (ii=0;ii<no_dof;ii++)
	{
	int i=no_dof-ii+1;
	for (int j=0;j<no_dof;j++)
		{
		double temp=0;
		for (int k=0;k<i;k++) temp+=CMU[k][j]*U[i][k];
		CMU[i][j]=temp;
		}
	}

}


 5.
{
Eigen2();
}
*/

/*
{
/
// Subprogram  TRAVEC for subprogram VIB
tolerance=1e-8;   				// Set tolerance for zero eigenvalues
if (idcmp=1)  goto 2; 			// if mass matrix was decomposed
for (i=0;i<no_dof;i++)                            
	cv[i]=1.0/fsqrt(cv[i]));  // Calculate angular freqencies ( 3.6-12)
2: 
for (i=0;i<no_dof;i++)
	if (cv[i]<tolerance) cv[i]=0;  // Set small eigen values = 0 
			else cv[i]=fsqrt(cv[i]);  // Calculate angular freqencies ( 3.6-12)

// Back-Transform Eigen  vectors
for (j=0;j<no_dof;j++)
for (i=0;i<no_dof;i++)
{
sum=0;
for (k=j;k<no_dof;k++)
	sum+=U[k][j]*phi[k][i]; // Back-transform eigenvectors ( 3.6-12)
phi[j][i]=sum;
}

// c. Normalize Eigenvectors with  respect to Largest values
for (i=0;i<no_dof;i++)
BIG=0;
for (j=0;j<no_dof;j++) //Search for largest vector in I
{
  c1=fabs(phi[j][i]); // Set up constants in c1,c2
  c2=fabs(big);
  if (c1>c2) big=c1;
}

for (j=0;j<no_dof;j++)
	phi[j][i]/=big; // Normalize vector I with respect BIG

}

// Subprogram Resipt for subprogram vibpt
// a. Reorder angular frequencies and Eigenvectors
if (idcmp!=1)      
	for (i=0;i<ndf/2;i++)  // if mass matrix was not decomposed
		{
		j=ndf-I+1;
		temp=cv[i];  //Interchange symmetrically placed angular frequencies
		cv[i]=cv[j];
		cv[j]=temp;
		for (k=0;k<no_dof;k++)
		   {
		   temp=phi[k][i]; //interchange symmetrically placed eigen vectors
		   phi[k][i]=phi[k][j];
		   phi[k][j]=temp;
		   }

// b. Write angular frequencies and Expanded eigen vectors
for (i=0;i<no_dof;i++)
for (j=0;j<no_dof;j++)
		if (k>ndf) EE[j]=0; // for restrained dof
			else EE[j]=phi[ID[J]][i]; // for free dof
// MODE(I)
// ANGULAR FREQUECY (CV[i])
// NODE DJ1 DJ2

// for (j=0;j<no_nodes;j++) cout << j,EE[2*j-1],EE[2*j];

// c.Normalize Eigenvectors with respect to Mass Matrix
for (i=0;i<no_dof;i++)
for (j=0;j<no_dof;j++)
{
sum=0;
for (k=0;k<no_dof;k++)
	sum+=phi[k][i]*CMS[k][j];  // Premultiply M by EV'
EE[j]=sum; // temp storage
sum=0;
for (j=0;j<no_dof;j++)
	sum+=EE[j]*phi[j][i]; // Post multiply by the eigen vector
sum=fsqrt(sum); // Calculate normalization constant
for (j=0;j<no_dof;j++)
	phi[j][i]/=sum; // Normalize eigen vector
}

// 8. Subprogram dylopt for program dynapt
// a. Dynamic Parameters
scanf(" ",isolve,nts,dampr);//  Read dynamic parameters
// b. Initail conditions
for (j=0;j<no_dof;j++)
	{d0[j]=0; v0[j]=0;}
scanf(" ",nind,ninv) // Read initial condition parameters

// Initial displacements
if (nnid!=0)
{
for (i=0;i<nnid;i++)
	scanf(" " ,j,D0[2*j-],D0[2*j]); // REad initial displacements
for (i=0;i<nnid;i++)
	{
	id[i];
	if (k<=no_dof) D0[k]=D0[i];
	}
}
// 2.Initail velocities
if (nniv!=0)
{
for (i=0;i<nniv;i++)
	scanf(" " ,j,V0[2*j-],V0[2*j]); // REad initial displacements
for (i=0;i<nniv;i++)
	{
	id[i];
	if (k<=no_dof) V0[k]=V0[i];
	}
}
// c. Applied Actions
{
for (j=0;j<nnd;j++)
	AS[j]=0; // Clear Action vector
scanf(" ",nln,nel); // read load parameters
//1. Nodal loads
if (nln!=0)
for (i=0;i<nln;i++)
	scanf(" " ,j,AS[2*j-],AS[2*j]); // Read scale factors for nodal loads
// 2) element loads
if (nel!=0)
for (j=0;j<nel;j++)
	{
	scanf(" " ,i,b[1],b[2],b[3],b[4]); // Read scale factors for line loads
	pbs[1]=(2*b[1]+b[3])*el[i]/6; // Calculate equivalent nodal loads (4.10-4)
	pbs[2]=(2*b[2]+b[3])*el[i]/6;
	pbs[3]=(2*b[3]+b[1])*el[i]/6;
	pbs[4]=(2*b[4]+b[2])*el[i]/6;
	j1=2*jn[i]-1;
	j2=j1+1;
	k1=2*kn[i]-1;
	k2=k1+1;
	AS[j1]+=PBS[1];
	AS[j2]+=PBS[2];
	AS[k1]+=PBS[3];
	AS[k2]+=PBS[4];
	}

// d. Ground Accelerations
scanf(" ",iga); // read acceleration parameter
if (iga!=0) 
{
scanf(" ",gax,gay);
// Fill vector GA with ground accelrations
for (i=0;i<no_nodes;i++)
	{
	j1=id[2*i-1];
	j2=id[2*i];
	if (j1<=no_dof) GA[j1]=GAX;
	if (j2<=no_dof) GA[j2]=GAY;
	}
// Multiply accelration with mass matrix
for (i=0;i<nnd;i++)
	{
	j=id[j];
	if j<=no_dof
		{
		sum=0;
		for (k=0;k<no_dof;k++)
			sum+=CMS[j][k]*GA[k];
		AS[i]-=sum;
		}
	}
}

// Read forcing function 
n=nln+nel+iga; // Add the parameters requiring a forcing function
nfo=0; 
if (n!=0) 
	{
	scanf(" " ,nfo); // Read function ordinates
	if( j=0;j<nfo;j++)
		scanf(" ",k,t[k],fo[k]);
	}
// f . Calculate step and ramp for each increment
	if( j=0;j<nts;j++)
		{
		step[j]=0;
		ramp[j]=0;
		}
	nfom1=nfo-1;
j=0;
for(k=0;k<nfom1;k++)
	{
ntss=(t[k+1]-t[k])/dt;
if (ntss!=0) 
	{
dstep=(fo[k+1]-fo[k])/ntss;
	for (i=0;i<ntss;i++)
		{
		j++;
step[j]=fo[k]+(i-1)*dstep;
ramp[j]=dstep;
}

// 9. Subprogram tranor for subprogram normod
// a. Read number of modes
scanf(" ",no_modes);
// Normal mode solution
// b. Calculate transformation operator 1/phiN
for( i=0;i<no_modes;i++)
for( j=0;j<no_dof;j++)
	{
sum=0;
for( k=0;k<no_dof;k++)
	{
sum+=phi[k][j]*CMS[k][j];
trop[i][j]=sum;
}
// Transform initial conditions
for( i=0;i<no_modes;i++)
{
dsum=0;
vsum=0;
for( j=0;j<no_dof;j++)
	{
	dsum+=trop[i][j]*D0[j];
	vsum+=trop[i][j]*V0[j];
	}
dn0[i]=dsum;
vn0[i]=vsum;
}

// d. Transform applied forces
for( i=0;i<no_modes;i++)
{
sum=0;
for( j=0;j<no_dof;j++)
	sum+=phi[j][i]*AF[j];
AN[i]=sum;
}
// 10 . Subprogram tihist for subprogram normod
for( i=0;i<no_modes;i++)
{
omega=cv[i];
	// response constants
if (!eq(omega,0) // for non rigid body mode 
	{
omega2=omega*omega;
omegad=omega*fsqrt(1-dampr*dampr);
cn=dampr*omega;
ex=exp(-cn*dt);
c=cos(omegad*dt);
s=sin(omegad*dt);
c1=ex*c+cn*s/omegad;
c2=es*s/omegad;
c3=(dt-c2-2.0*cn*c3)/(omega2*dt);
c5=-omega*c2;
c6=c1-2.0*cn*c2;
c7=c2;
c8=c3/dt;
}
else
{
c1=\0;
c2=dt;
c3=dt*dt/2.0;
c4=dt*dt/6.0;
c5=0;
c6=1;
c7=c2;
c8=c3/dt;
}
dnj=dn0[i];
vnj=vn0[i];
for(j=0;j<nts;j++)
	{
dn[i][j]=c1*dnj+c2*vnj+c3*an[i]*step[j]+c4*an[i]*ramp[j];
vnj=c5*dnj+c6*vnj+c7*an[i]*step[j]+c8*an[i]*ramp[j];
dnj=dn[i][j];
	}
}

// 11. Subprogram for trabac for subprogram normod
for (j=0;j<nts;j++)
for (k=0;k<no_dof;k++)
	{
	sum=0;
	for (i=0;i<no_modes;i++)
		sum+=phi[k][i]*dn[i][j];
	df[k][j]=sum;
	}
//12.subprogram numint for subprogram dyna
scanf(" ",alpha,beta,gamma);
// Direct numerical intergration
// b. Calculate strucure damping matrix
// determine factor SA
for (i=0;i<no_dof;i++)
	{
c1=fsqrt(2*dampr*cv[i]);
for (j=0;j<no_dof;j++)
	{
sum=0;
	for (k=0;k<no_dof;k++)
		sum+=CMS[j][k]*phi[k][j];
	SA[j][i]=sum*c1;
}
// multiplty SA with its transpose
for (k=0;k<no_dof;k++)
for (j=0;j<no_dof;j++)
	{
sum=0;
	for (i=0;i<no_dof;i++)
		sum+=SA[k][i]*SA[j][i];
	CS[k][j]=sum;
	}
// c. Calculate initial accelerations

for (i=0;i<no_dof;i++)
for (j=0;j<no_dof;j++)
	SA[i][j]=CMS[i][j];
FACTOR(SA,NO_DOF,&8); //Factor mass matrix, check for singularity 

for (i=0;i<no_dof;i++)
	{
	sum=0;
for (j=0;j<no_dof;j++)
	sum-=ss[i][j]*d0[j]-cs[i][j]*v0[j]; // Set up right hand side vector
aa[i]=sum+AF[i]*step[i];
}
SOLVER(SA,aa,a0,no_dof); // solve (5.1-4) for initial accelerations
for (i=0;i<no_dof;i++) 
	a0[i]=0; // if mass matrix is singular use zero initial displacements

// d . Calculate and factor equivalent stiffness matrix

c1=1.0+alpha;  // determine constants for eqns. (5.5-15) and (5.5-33)
c2=1.0/(beta*dt*dt);
c3=gamma(beta*dt);
for (i=0;i<no_dof;i++)
for (j=0;j<no_dof;j++)
	SA[i][j]=c1*ss[i][j]+c2*CMS[i][j]+c3*cs[i][j];
FACTOR(sa,no_dof,&101);
// Calculate response for first time step
c1=1.0/beta*dt);
c2=1.0/(2*beta);
c3=gamma/beta;
c4=(gamma*c2-1.0)*dt;
c5=c1/dt;
c6=gamma*c1;
for (i=0;i<no_dof;i++)
	{
	Q[i]=c1*vo[i]+c2*a0[i];
	R[i]=c3*v0[i]+c4*a0[i];
	}
for (i=0;i<no_dof;i++)
	{
	sum=0;
	for (k=0;k<no_dof;k++)
		sum+=CMS[i][k]*Q[k]+CS[i][k]*R[k]; // add mass and damping terms in 5.5-16
	aa[i]=sum+AF[i]*ramp[i];
SOLVER(sa,aa,dd,no_dof);
for (i=0;i<no_dof;i++)
	{
	df[i][1]=d0[i]+dd[i];
	a0[i]+=c5*dd[i]-Q[i];
	v0[i]+=c6*dd[i]-R[i];
	}
// Calculate response for subsequent time steps
for (j=1;j<nts;j++)
for (i=0;i<no_dof;i++)
	{
	Q[i]=c1*v0[i]+c2*a0[i];
	R[i]=c3*v0[i]+c4*a0[i];
	}
for (i=0;i<no_dof;i++)
	{
	sum=0;
	for (k=0;k<no_dof;k++)
		sum+=CMS[i][k]*Q[k]+CS[i][k]*R[k]+alpha*ss[i][k]*dd[k]; // add mass and damping terms in 5.5-16

	aa[i]=sum+AF[i]*(step[j]+*ramp[j]-step[j-1]-ramp[j-1]));
	}

SOLVER(sa,aa,dd,no_dof);

for (i=0;i<no_dof;i++)
	{
	df[i][j]=df[i][j-1]+dd[i];
	a0[i]+=c5*dd[i]-Q[i];
	v0[i]+=c6*dd[i]-R[i];
	}
}
// 13. subprogram res2pt for subprogram dynapt
// a.read output selections
scanf(" ",iwr,ipl,nno,neo);
for (j=0;j<nno;j++)
	scanf(" ",jno[j]);
for (j=0;j<nno;j++)
	scanf(" ",ieo[j]);
if (ipl==1) // write to plotter

// b. Transfer displacements from DF to DS
for (i=0;i<nnd;i++)
	{
	k=id[i]; 
	if (k>no_dof)
		{
		ds0[i]=0;
		for (j=0;j<nts;j++) ds[i][j]=0;
		}
	else
		{
		ds0[i]=d0[k];
		for (j=0;j<nts;j++)
			ds[i][j]=df[k][j];
		}
	}

// write displacements
for (kk=0;kk<nno;kk++)
	{
	k=jno[kk];
	k0=2*(k-1);
	j=nrl[k0+1]+nrl[k0+2];
	if (j==2) continue;
	j=0;
	time=0.0;
	if (iwr==1) 
		for (k=0;k<2;k++) printf(" " ,j,time,(ds0(k0+k));
	// Some output (?) code left here.

// Calculate and write axial forces
for (ii=0;ii<neo;kk++)
	{
	i=ieo[ii];
	j1=2*jn[i]-1;
	j2=j1+1;
	k1=2*kn[i]+1;
	k2=k1+1;
	scm=e*ax[i]/el[i];
	j=0;
	time=0;
	amd=scm*((ds0[k1])-ds0[j1])*cx[i]+(ds0[k2]-ds0[j2])*cy[i]);
	amdmax=amd;
	amdin=amd;
	tmdmax=0;
	tmdmin=0;
	for (j=0;j<nts;j++)
		{
		time=dt*j;
		amd=scm*((ds[k1][j])-DS[j1][j])*cx[i]+(ds[k2][j]-ds[j2][j])*cy[i]);
		if (amd >= amdmax) 
			{
			amdmax=amd; 
			tmdmax=time;
			}
		
		if (amd <= amdmin) 
			{
			amdmin=amd; 
			tmdmin=time;
			}
		}


void Cstructure::decompose(double &&a)      // Decompose M or K
{
for (int i=0;i<no_dof;i++)
	for (int j=0;j<no_dof;j++)
		{
			double sum=a[i][j];
			for (int k=0;k<i-1;k++)	sum-=a[k][i]*a[k][j];
			if (i==j)
				{
					if (sum<=0)  return;
					a[i][i]=sqrt(sum);
				}
				else
					a[i][j]=sum/a[i][i];
		}
}


// Subprogram FACTOR for subprogram numint
void Cstructure::factor()
{
if a[1][1]<=0 return;
for (j=1;j<n;j++)
	{
	//j1=j-1;
	for (i=1;i<j-1;i++)
		{
		sum=a[i][j];
		ii=i-1;
		for (k=0;k<ii;k++)
			sum-=a[k][i]*a[k][j];
		a[i][j]=sum;
		}
	sum=a[j][j];
	for (k=0;k<j1;k++)
		{
		temp=a[k][j]/a[k][k];
		sum-=temp*a[k][j];
		a[k][j]=temp;
		}
	if (sum<=0) return;
	a[i][j]=sum;
	}
}

// 17. subprogram solver for subprogram numint

void Cstructure::solver()
{
for (int i=0;i<no_dof;i++)
	{
	double sum=b[i];
	//k1=i-1;
	for (int k=0;k<i-1;k++)
		sum-=U[k][i]*x[k];
	x[i]=sum;
	for (i=0;i<n;i++)
		x[i]/=U[i][i];
	for (int ii=0;ii<n;ii++)
		{
		i=n-ii+1;
		sum=x[i];
		//k2=i+1;
		for (k=i+1;k<n;k++)
			{
			sum-=U[i][k]*x[k];
			x[i]=sum;
			}
		}
	}
}
	
}

*/

void Make_Unity(double **A,int n)
{
for (int i=0;i<n;i++)  
for (int j=0;j<n;j++)  
	if (i==j) A[i][j] = 1.0; else  A[i][j] = 0.0;
}

BOOL Cstructure::compute_EigenVectors()
//BOOL jacobi()
{		
   int nsmax = 15; 						// maximum number of iterations
/*   double d[40],eptol_K,eptol_M,eps,akk,ajj,ab,
   check,sqch,d1,d2,den,cg,ca,ak,bk,epsa,epsb,aj,bj,tol,rtol,dif,bb;
   int    i,j,k,n,nsweep,nr,jj,jm1,kp1,km1,jp1;
*/

no_dof=4;

K[0][0]=5;
K[1][0]=-4;
K[2][0]=1;
K[3][0]=0;

K[1][1]=6;
K[2][1]=-4;
K[3][1]=1;

K[2][2]=6;
K[3][2]=-4;

K[3][3]=5;

Make_Unity(M,4);

int n = no_dof ;
double rtol = 1e-20;

Fill_Upper(K,n);
Fill_Upper(M,n);

matout(K,n,n);
matout(M,n,n);

for(int i = 0;i < n;++i)
    	{
	  	if ((K[i][i] <= 0.0) || (M[i][i] <= 0.0))
			{
			fprintf(outfile,"MATRICES NOT POSITIVE DEFINITE\n");
			return 1;		  
			}
 		eig_val[i] = d[i] = K[i][i] / M[i][i];
  		}

Make_Unity(eig_vect,n);

if (n == 1) return(0);

/* Initialize Sweep Counter and Begin Iteration */

   int nsweep = 0;
   int nr = n - 1;

do {
   nsweep++;

   fprintf(outfile,"\nSWEEP NUMBER IN *** JACOBI *** = %d\n",nsweep);

/* Check if present off-diagonal element is large enough to require zeroing */
 
   double eps = sqr(pow(0.01,nsweep));
   for(int j = 0;j < n-1;++j)
     {
	// int jj = j+1;
	 for(int k = j+1;k < n;++k)
	    {
	      double eptol_K = K[j][k]*K[j][k] / (K[j][j]*K[k][k]);
	      double eptol_M = M[j][k]*M[j][k] / (M[j][j]*M[k][k]);

	      if (eptol_K < eps && eptol_M < eps) continue;

		/* If zeroing is required calculate the rotation matrix elements ca and cg */

	      double akk = K[k][k] * M[j][k] - M[k][k] * K[j][k];
	      double ajj = K[j][j] * M[j][k] - M[j][j] * K[j][k];
	      double ab  = K[j][j] * M[k][k] - K[k][k] * M[j][j];

	      double check = (ab * ab + 4.0 * akk * ajj) / 4.0;
	      
	      if(check <0.0)
				{
				fprintf(outfile,"\n**** ERROR : MATRICES NOT POSITIVE DEFINITE\n");
				exit(0);
				}

	      double sqch = sqrt(check);
	      
	      double d1 = ab / 2.0 + sqch;
	      double d2 = ab / 2.0 - sqch;
	      double den=d1;
	    
	      if (fabs(d2) > fabs(d1)) den = d2;

	      double ca,cg;
	      if (den==0)
					{
					ca = 0.0;
					cg = -K[j][k] / K[k][k];
					}
				else
					{
					ca = akk / den;
					cg = -ajj / den;
					}

/* Perform the Generalized Rotation to Zero the Present Off-diagonal Element */

            double kjj = K[j][j] + cg * cg * K[k][k] + 2 * cg * K[j][k];
            double kjk = K[j][k] + ca * K[j][j]  + cg * K[k][k] +  ca * cg * K[j][k];
            double kkk = K[k][k] + ca * ca * K[j][j] + 2 * ca * K[j][k] ;
            K[j][j] = kjj;
            K[k][k] = kkk;
            K[j][k] = kjk;
            K[k][j] = K[j][k];
            
            kjj = M[j][j] + cg * cg * M[k][k] + 2 * cg * M[j][k];
            kjk = M[j][k] + ca * M[j][j]  + cg * M[k][k] +  ca * cg * M[j][k];
            kkk = M[k][k] + ca * ca * M[j][j] + 2 * ca * M[j][k] ;

            M[j][j] = kjj;
            M[k][k] = kkk;
            M[j][k] = kjk;
            M[k][j] = M[j][k];

            for(i = 0;  i<n; i++)
	  	       if(i !=j && i != k)
			       {
		             double  kji = K[j][i] + cg*K[k][i];
		  	         double  kki = K[k][i] + ca*K[j][i];
		       
		             K[j][i] = kji;
		 	         K[k][i] = kki;
		 	   
		 	         double  kij = K[i][j] + cg*K[i][k];
			         double  kik = K[i][k] + ca*K[i][j];
		       
		             K[i][j] = kij;
		 	         K[i][k] = kik;
			   
			         kji = M[j][i] + cg*M[k][i];
			         kki = M[k][i] + ca*M[j][i];
		       
		             M[j][i] = kji;
		 	         M[k][i] = kki;
			   
			         kij = M[i][j] + cg*M[i][k];
			         kik = M[i][k] + ca*M[i][j];
		       
		             M[i][j] = kij;
			         M[i][k] = kik;
			      }

/*      if (n-2 != 0)
	      {
	      int jp1 = j + 1;
	      int jm1 = j - 1;
	      int kp1 = k + 1;
	      int km1 = k - 1;
	      if ( jm1-1 >= 0)   // 1 >0
	      for(i = 0;i < jm1;++i)
				 {
				    double aj = K[i][j];
				    double bj = M[i][j];
				    double ak = K[i][k];
				    double bk = M[i][k];
				    K[i][j] = aj + cg * ak;
				    M[i][j] = bj + cg * bk;
				    K[i][k] = ak + ca * aj;
				    M[i][k] = bk + ca * bj;
				 }

	      if ( (kp1-n) <= -1)
	      for(i = kp1;i < n;++i)
				 {
				    double aj = K[j][i];
				    double bj = M[j][i];
				    double ak = K[k][i];
				    double bk = M[k][i];
				    K[j][i] = aj + cg * ak;
				    M[j][i] = bj + cg * bk;
				    K[k][i] = ak + ca * aj;
				    M[k][i] = bk + ca * bj;
				 }
	      
	      if ( (jp1-km1) <= 0)
	      for(int i = jp1;i < km1;++i)
				 {
				    double aj = K[j][i];
				    double bj = M[j][i];
				    double ak = K[i][k];
				    double bk = M[i][k];
				    K[j][i] = aj + cg * ak;
				    M[j][i] = bj + cg * bk;
				    K[i][k] = ak + ca * aj;
				    M[i][k] = bk + ca * bj;
				 }
	      }         // end if ((n-2)!=0)
	      
	      double ak = K[k][k];
	      double bk = M[k][k];
	      K[k][k] = ak + 2.0 * ca * K[j][k] + ca * ca * K[j][j];
	      M[k][k] = bk + 2.0 * ca * M[j][k] + ca * ca * M[j][j];
	      K[j][j] = K[j][j] + 2.0 * cg * K[j][k] + cg * cg * ak;
	      M[j][j] = M[j][j] + 2.0 * cg * M[j][k] + cg * cg * bk;
	      K[j][k] = 0.0;
	      M[j][k] = 0.0;

*/

/* Update the Eigen Values After Each Rotation */
	      for(int i = 0;i < n;++i)
				 {
				 	double tj=eig_vect[i][j];
					double tk=eig_vect[i][k];
				    eig_vect[i][j] = tj + cg * tk;
				    eig_vect[i][k] = tk + ca * tj;
				 }

      }       // k
   } 	// j

/* Update the Eigen Values After Each Sweep */

for(i = 0;i < n;++i)
		    if (K[i][i] <= 0.0 || M[i][i] <= 0.0) 
		    		{
		    		fprintf(outfile,"\n**** ERROR : MATRICES NOT POSITIVE DEFINITE\n");
		    		exit(0);
					}
		    else 
		    eig_val[i] = K[i][i] / M[i][i];
	    
fprintf(outfile,"\nCURRENT EIGEN VALUES IN ** JACOBI ** ARE\n\n");
for(i = 0;i < n; ++i) fprintf(outfile,"%10.4lf      ",eig_val[i]);

fprintf(outfile,"\nSTIFFNESS MATRIX\n");
matout(K,n,n);
fprintf(outfile,"\nMASS MATRIX\n");
matout(M,n,n);

/* Check for Convergence */
double dif,tol;

for(i = 0;i < n;++i)
		 {
		    tol = fabs(rtol * d[i]);
		    dif = fabs( eig_val[i] - d[i]);
		    if ( dif > tol) break;
		 }

double epsa,epsb;

if (dif<=tol)
	 {
/* Check all Off-diagonal Elements to see if Another Sweep is Required */
      eps = pow(rtol,2.0);
      for(j = 0;j < nr;++j)
		 {
		    int jj = j + 1;
		    for(int k = jj;k < n;++k)
		       {
				  epsa = fabs((K[j][k] * K[j][k]) / (K[j][j] * K[k][k]));
				  epsb = fabs((M[j][k] * M[j][k]) / (M[j][j] * M[k][k]));
				  if (epsa >= eps || epsb >= eps) break;
		       }
		 }

OVER:  if (epsa < eps && epsb < eps)
	{
	/* Fill the Bottom Triangle of Resultant Matrices and Scale Eigen Vectors */
      for(i = 0;i < n;++i)
		    for(j = 0;j < n;++j)
		       {
				  K[j][i] = K[i][j];
				  M[j][i] = M[i][j];
		       }
		 
	      for(j = 0;j < n;++j)
			 {
			    double temp = sqrt(M[j][j]);
			    for(int k = 0;k < n;++k) 
			    		eig_vect[k][j] /= temp;
			 }

// Sort eigen values in ascending order

   for(j=1; j<n; j++)
   for(i=0; i<n-j; i++)
	 if( eig_val[i] > eig_val[i+1] )
         {
			double kkj = eig_val[i];
			eig_val[i] = eig_val[i+1];
			eig_val[i+1] = kkj;
			for(int k=0; k<n; k++)
				{
				kkj = eig_vect[k][i];
				eig_vect[k][i] = eig_vect[k][i+1];
				eig_vect[k][i+1] = kkj;
				}
         }
 
    /*
    for(j=0; j<n; j++)
	   { 
	     if	(eig_val[j]<0) eig_val[j] = 0;
		 				else     
		 		eig_val[j] = sqrt(eig_val[j])/(2.0*_PI);
	   }
	*/

     return(0);
     }
 }

if (nsweep>=nsmax) goto OVER;

/* Update d Matrix and Start new Sweep */
for (i = 0;i < n;++i)	d[i] = eig_val[i];

} while (1);
}





