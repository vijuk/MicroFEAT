// plate8()
switch (ind) {
	case 1:	scanf(nel);
	for(i=1;i<=no_nodes;i++) ELC[i]=0;
	kk=0;
	for(i=1;i<=nel;i++) 
		{
		scanf(lnum,mg,(lnc[j],j=1,4),th,udl,ki)
		kk++;
		}

	// Generate intermediate nodes ( not coded )
	// Calculate element stiffness matrix EK
	R={-1,1,1,-1,0,1,0,-1}
	S={-1,-1,1,1,-1,0,1,0}
	U={-1,1,1,-1}
	V={-1,-1,1,1}
	WG={1,1,1,1}

	for (i=0;i<5;i++)
	for (j=0;j<5;j++)
		C[i][j]=0;

	// Compute constitutive relationship matrix

	C[1][1]=E[MG]*th*th*th/(12.0*(1-pr[mg]*pr[mg]));
	C[1][2]=C[1][1]*pr[mg];
	C[2][1]=C[1][2];
	C[2][2]=C[1][1];
	C[3][3]=C[1][1]*(1-pr[mg])/2.0;
	C[4][4]=E[mg]*th/(2.4*(1+pr[mg]);
	C[5][5]=C[4][4];

	// Initialize element stiffness matrix and load vector
	for (i=0;i<24;i++)
	{
	Q[i]=0;
	for (j=0;j<24;j++)
		EK[i][j]=0;
	}
	// Shear energy terms EKS (1x1 integration is used)
	for (i=0;i<12;i++)
	for (j=0;j<12;j++)
		EKS[i][j]=0;
		
	// Enter loops for numerical integration (2x2 used)
	for (ix=1;ix<=2;ix++)
	for (iy=1;iy<=2;iy++)
	// Calculate shape functions, derivatives,, Jacobian,
	// Determinant Jacobian, Jacobian Inverse

	// Evaluate shape functions
	for(i=1;i<=8;i++)
		{
		aa=(1.0+r[i]*gp[ix]);
		bb=(1.0+s[i]*gp[iy]);

		if (i<4)
			{
		sf[i]=0.25*aa*bb*(aa+bb-3.0);
		sf[1][i]=0.25*(2.0*aa+bb-3.0)*bb*r[i];
		sf[2][i]=0.25*(2.0*aa+bb-3.0)*aa*s[i];
			}
			else
			{
		aa+=(r[i]*2.0-1.0)*gp[ix]*gp[ix];
		bb+=(s[i]*2.0-1.0)*gp[iy]*gp[iy];
		sf[i]=0.5*aa*bb;
		sf[1][i]=0.5*(r[i]+2*(r[i]*r[i]-1)*gp[ix])*bb;
		sf[2][i]=0.5*(s[i]+2*(s[i]*s[i]-1)*gp[iy])*aa;
			}
		}

	//Compute Jacobian Matrix
	for (i=0;i<2;i++)
	for (j=0;j<2;j++)
		JM[i][j]=0;

	for (i=0;i<8;i++)
	for (k=0;k<2;k++)
		{
		JM[k][1]+=xL[i]*sfd[k][i];
		JM[k][2]+=yL[i]*sfd[k][i];
		}

	// Calcualte determinant of Jacobian Matrix
	DET_J=JM[1][1]*JM[2][2]-JM[1][2]*JM[2][1];
	// Compute Inverse of Jacobian
	Inv_JM[1][1]=JM[2][2]/DET_J;
	Inv_JM[1][2]=-JM[1][2]/DET_J;
	Inv_JM[2][1]=-JM[2][1]/DET_J;
	Inv_JM[2][2]=JM[1][1]/DET_J;
	da=DET_J*wg[ix]*wg[iy];

	// Form global derivatives SFG

	for (i=0;i<2;i++)
	for (j=0;j<8;j++)
		SFG[i][j]=0;
		
	for (i=0;i<2;i++)
	for (k=0;k<8;k++)
	for (j=0;j<2;j++)
		SFG[i][k]+=Inv_JM[i][j]*N[j][k];

	// Compute CBS - the shear contribution to stress matrix

	for (i=0;i<5;i++)
	for (j=0;j<24;j++)
		B[i][j]=0;

	for (i=0;i<8;i++)
		{
		k1=3*(i-1)+1;
		k2=k1+1;
		k3=k2+1;
		B[1][k3]= sfg[1][i];
		B[2][k2]=-sfg[2][i];
		B[3][k3]= sfg[2][i];
		B[3][k2]=-sfg[1][i];

		B[4][k1]=-sfg[1][i];
		B[4][k3]=-sf[i];
		B[5][k1]=-sfg[2][i];
		B[5][k2]= sfg[i];
		}


	for (i=1;i<=5;i++)
	for (j=0;j<24;j++)
		{
		CB[i][j]=C[i][i]*B[i][j];
		if(i<=2)
			{
			ii=3-i;
			CB[i][j]+=C[i][ii]*B[ii][j];
			}
		}

	// Calculate the element stiffness matrix

	for(i=1;i<=24;i++)
	for(j=1;j<=24;j++)
	for(k=1;k<=5;k++)
		EK[i][j]+=B[k][i]*CB[k][j]*da;
	
	// Calculation of load vector due to UDL


	for(i=1;i<=8;i++)
		{
		nn=3*i-2;
		Q[nn]+=SF[i]*UDL[*DA;
		}
	for(i=2;i<=24;i++)
		{
		k=i-1;
	for(j=1;j<=k;j++)
		EK[i][j]=EK[j][i];
		}


for (j=1;j<=no_elem_dof;j++)
	{
	jj=nd[j];
	if (jj<=0) PQ[jj]+=Q[j];
	}
	nued=no_elem_dof;
	passem(sk,ek,nds,nd,ned,neq1,nsky,nued);
}

// Calculation of bending stress resultants at nodes

val=3.0;
term=sqrt(val)/2.0;
tms[1][1]=1.0+term;
tms[1][2]=-0.5;
tms[1][3]=-0.5;
tms[1][4]=1.0-term;

tms[2][1]=-0.5;
tms[2][2]=tms[1][4];
tms[2][3]=tms[1][1];
tms[2][4]=-0.5;

tms[3][1]=tms[2][2];
tms[3][2]=-0.5;
tms[3][3]=-0.5;
tms[3][4]=tms[1][1];

tms[4][1]=-0.5;
tms[4][2]=tms[1][1];
tms[4][3]=tms[2][2];
tms[4][4]=-0.5;

for (i=1;i<=nsn;i++)
	for(j=1;j<=ns;j++)
		str[i][j]=0.0;

for (elem_no=1;elem_no<=nel;elem_no++)
{

// retrive variable from istress

for (i=1;i<=ns;i++)
		sig[i]=0;
for(j=1;j<=no_elem_dof;j++)
	{
	jj=nd[j];
	if (jj!=0)
		for (i=1;i<=ns;i++)
			sig[i]+=CB[i][j]*P[jj];
	}

// Shear stress resultants within an element at R=0 , S=0
if (k==1) printf("");
for (j=1;j<=ns;j++)
	for (i=1;i<=4;i++)
		{
		nn=lnc[i];
		str[nn][j]+=tms[i][k]*sig[j]/elc[nn];
		}
	}


	



