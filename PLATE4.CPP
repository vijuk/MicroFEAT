// plate4()

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
	R={-1,1,1,-1}
	S={-1,-1,1,1,}
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
	for (i=0;i<12;i++)
	{
	Q[i]=0;
	for (j=0;j<12;j++)
		EK[i][j]=0;
	}
	// Shear energy terms EKS (1x1 integration is used)
	for (i=0;i<12;i++)
	for (j=0;j<12;j++)
		EKS[i][j]=0;
		
	gpsx=0;
	gpsy=0;
	// Evaluate shape functions
	for(i=1;i<=4;i++)
		{
		aa=(1.0+r[i]*gpsx);
		bb=(1.0+s[i]*gpsy);
		sf[1][i]=0.25*r[i]*bb;
		sf[2][i]=0.25*s[i]*aa;
		sf[3][i]=0.25*aa*bb;
		}
	//Compute Jacobian Matrix
	for (i=0;i<2;i++)
	for (j=0;j<2;j++)
		JM[i][j]=0;

	for (i=0;i<4;i++)
	for (k=0;k<2;k++)
		{
		JM[k][1]+=xL[i]*sf[k][i];
		JM[k][2]+=yL[i]*sf[k][i];
		}
	// Calcualte determinant of Jacobian Matrix
	DET_J=JM[1][1]*JM[2][2]-JM[1][2]*JM[2][1];
	// Compute Inverse of Jacobian
	Inv_JM[1][1]=JM[2][2]/DET_J;
	Inv_JM[1][2]=-JM[1][2]/DET_J;
	Inv_JM[2][1]=-JM[2][1]/DET_J;
	Inv_JM[2][2]=JM[1][1]/DET_J;
	db=4*DET_J;
	// Form global derivatives SFG
	for (i=0;i<2;i++)
	for (j=0;j<4;j++)
		SFG[i][j]=0;
		
	for (i=0;i<2;i++)
	for (k=0;k<4;k++)
	for (j=0;j<4;j++)
		SFG[i][k]+=Inv_JM[i][j]*N[j][k];

	for (i=0;i<4;i++)
		SFG[3][i]=SF[3][i];
	// Compute CBS - the shear contribution to stress matrix
	for (i=0;i<NS;i++)
	for (j=0;j<12;j++)
		B[i][j]=0;
	
	for (i=0;i<4;i++)
		{
		k1=3*(i-1)+1;
		k2=k1+1;
		k3=k2+1;
		B[4][k1]=sfg[1][i];
		B[4][k3]=sfg[3][i];
		B[5][k1]=sfg[2][i];
		B[5][k2]=sfg[3][i];
		}

	for (i=1;i<=3;i++)
	for (j=0;j<12;j++)
		CBS[i][j]=0;
	for (i=4;i<=5;i++)
	for (j=0;j<12;j++)
		CBS[i][j]=C[4][4]*B[i][j];

	ii=0;
	for(ia=1;ia<=4;ia++)
	{
	jj=0;
	for(ib=1;ib<=4;ib++)
	{
	for(i=1;i<=2;i++)
		{
		EKS[ii+1][jj+1]+=C[4][4]*SFG[i][ia]*sfg[i][ib]*db;
		EKS[ii+i][jj+3-i]-=C[4][4]*SFG[i+1][ia]*sfg[4-i][ib]*db;
		EKS[ii+i+1][jj+i+1]+=C[4][4]*SFG[3][ia]*sfg[3][ib]*db;
		}
	
	for(i=1;i<=3;i+=2)
		EKS[ii+i][jj+4-i]+=C[4][4]*SFG[i][ia]*sfg[4-i][ib]*db;
	jj+=3;
	}
	ii+=3;
	}

// Initialize bending stiffness matrix 'EKB' ( 2*2) integration is used

	for (i=1;i<=12;i++)
	for (j=0;j<12;j++)
		EBK[i][j]=0;

	// Enter loops for numerical integration
	nsp=4;
	for (ix=1;ix<=nsp;ix++)
	// Evaluate shape functions, its derivatives, Jacobian and 
	// Determinant of Jacobian
	GP=0.5773502691896;

	// Shape functions
	for(i=1;i<=4;i++)
		{
		aa=(1.0+r[i]*u[ix]*gp);
		bb=(1.0+s[i]*v[ix]*gp);
		sf[1][i]=0.25*r[i]*bb;
		sf[2][i]=0.25*s[i]*aa;
		sf[3][i]=0.25*aa*bb;
		}
	//Compute Jacobian Matrix
	for (i=0;i<2;i++)
	for (j=0;j<2;j++)
		JM[i][j]=0;

	for (i=0;i<4;i++)
	for (k=0;k<2;k++)
		{
		JM[k][1]+=xL[i]*sf[k][i];
		JM[k][2]+=yL[i]*sf[k][i];
		}
	// Calcualte determinant of Jacobian Matrix
	DET_J=JM[1][1]*JM[2][2]-JM[1][2]*JM[2][1];
	// Compute Inverse of Jacobian
	Inv_JM[1][1]=JM[2][2]/DET_J;
	Inv_JM[1][2]=-JM[1][2]/DET_J;
	Inv_JM[2][1]=-JM[2][1]/DET_J;
	Inv_JM[2][2]=JM[1][1]/DET_J;
	da=WG[ix]*DET_J;
	// Form global derivatives SFG
	for (i=0;i<2;i++)
	for (j=0;j<4;j++)
		SFG[i][j]=0;
		
	for (i=0;i<2;i++)
	for (k=0;k<4;k++)
	for (j=0;j<4;j++)
		SFG[i][k]+=Inv_JM[i][j]*N[j][k];

	for (i=0;i<4;i++)
		SFG[3][i]=SF[3][i];
	// Compute CBS - the shear contribution to stress matrix
	for (i=0;i<NS;i++)
	for (j=0;j<12;j++)
		B[i][j]=0;
	
	for (i=0;i<4;i++)
		{
		k1=3*(i-1)+1;
		k2=k1+1;
		k3=k2+1;
		B[1][k3]=sfg[1][i];
		B[2][k2]=-sfg[2][i];
		B[3][k3]=sfg[2][i];
		B[3][k2]=-sfg[3][i];
		}


	for (i=1;i<=3;i++)
	for (j=0;j<12;j++)
		CBB[i][j]=C[i][i]*B[i][j];

	if (i<=2)
		{
		ii=3-i;
		CBB[i][j]=C[i][ii]*B[ii][j]+CBB[i][j];
		}

	for (i=4;i<=5;i++)
	for (j=0;j<12;j++)
		CBS[i][j]=0;

	// Compute the stress matrix CB

	for (i=1;i<=NS;i++)
	for (j=0;j<12;j++)
		CB[i][j]=0;

	for (i=1;i<=NS;i++)
	for (j=0;j<12;j++)
		{
		if (ix!=1) CB[i][j]+=CBB[i][j];
			else
			CB[i][j]+=CBB[i][j]+CBS[i][j];

		}

	ii=0;
	for(ia=1;ia<=4;ia++)
	{
	jj=0;
	for(ib=1;ib<=4;ib++)
	{
	for(i=2;i<=3;i++)
		{
		EKB[ii+i][jj+i]+=C[1][1]*SFG[4-i][ia]*sfg[4-i][ib]*da
			+C[3][3]*SFG[i-1][ia]*SFG[i-1][ib]*da;
		EKB[ii+i][jj+5-i]-=C[1][2]*SFG[4-i][ia]*sfg[i-1][ib]
			+C[3][3]*SFG[i-1][ia]*SFG[4-i][ib]*da;
		}
	
	jj+=3;
	}
	ii+=3;
	}


// Calculate element load vector due to UDL

for (i=1;i<=4;i++)
	{
	nn=3*i-2;
	Q[nn]+=SFG[3][i]*UDL*da;
	}

// Compute element stiffness matrix EK
for (i=1;i<=12;i++)
for (j=1;j<=12;j++)
	EK[i][j]+=EKB[i][j]+EKS[i][j];

for (i=2;i<=12;i++)
for (j=1;j<=i-1;j++)
	EK[i][j]=EK[j][i];

for (j=1;j<=no_elem_dof;j++)
	{
	jj=nd[j];
	if (jj>0) PQ[jj]+=Q[j];
	}
nued=no_elem_dof;
passem(sk,ek,nds,nd,ned,neq1,nsky,nued);
return;
}

// Calculation of bending stress resultants at nodes
val=3.0;
term=sqrt(val)/2.0;
tms[1][1]=1.0+term;
tms[1][2]=-0.5;
tms[1][3]=1.0+term;
tms[1][4]=-0.5;

tms[2][1]=-0.5;
tms[2][2]=1.0+term;
tms[2][3]=-0.5;
tms[2][4]=1.0+term;

tms[3][1]=1.0+term;
tms[3][2]=-0.5;
tms[3][3]=1.0+term;
tms[3][4]=-0.5;

tms[4][1]=-0.5;
tms[4][2]=1.0+term;
tms[4][3]=-0.5;
tms[4][4]=1.0+term;

for ( i=1;i<=nsn;i++)
	for(j=1;j<=3;j++)
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
for (j=1;j<=3;j++)
	for (i=1;i<=3;i++)
		{
		nn=lnc[i];
		str[nn][j]+=tms[i][k]*sig[j]/elc[nn];
		}
	}



	



