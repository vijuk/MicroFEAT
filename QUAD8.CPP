#include "stdafx.h"
#include "space.h"
/*
INPUT FOR Cquad8 element

{
read NPUT ; // if NPUT <> 0  stress averaging is done
read no_Mat

// Material ID, Youngs modulus,Poisson's ratio, unit weight

LOOP i=1,no_Mat 
	read Mat_ID,Y,meu,uwt 
READ node_no,RAD,ANG,COX,COY  // 

LOOP node_no,disx,disy // specified displacements ( only when dof is constraint )

LOOP node_no,rex,rey,edgep // nodal loads (already incorporated)
							// and edge load at node

wirl,gravx,gravy // angular rotation of structure about global element,
				// components of gravity in x & y dirn.s

*/



void Cquad8::Read_Properties()	// Added on 08/10/96.
{
Celem::Read_Properties();
//fscanf(inpfile," %

fscanf(inpfile,"%lf",&Thickness);
fprintf(outfile," %10.4lf\n",Thickness);
}

void Cquad8::Compute_C(double &c1,double &c2, double &c3)
{
	// sets up material property matrix
	double E,MEU,WD,G;
	material->get(E,MEU,WD,G);

	C[3][0]=C[0][0]=E/(1.0-MEU*MEU);
	C[2][0]=C[1][0]=C[0][0]*MEU;
	C[4][0]=G;
}


void Cquad8::Compute_SMs_Qbody(double cc1,double cc2,double cc3)
{
double x1,y1,z1,x2,y2,z2;

int Gauss_Order = 3;

for (int i=0;i<no_elem_dof;i++) AML[i]=0.0;         // element load vector  = 0
for (i=0;i<136;i++) S[i]=0.0;		      // element stiffness    = 0

for (i=0;i<no_elem_nodes;i++)  
	{
	nodeman->getnode(enode[0]->no,x1,y1,z1);
	nodeman->getnode(enode[1]->no,x2,y2,z2);

	// intcrd();   // midside nodes are generated, if not defined.
	//pload();      // checks for any nodal loads for element

	for (int j=0;j<Gauss_Order;j++)  
		{
		rst[0]=Gauss_Point[j][Gauss_Order-1];			  	// psi: local coordinate 1
		double wx=Gauss_Weight[j][Gauss_Order-1];
		for (int k=0;k<Gauss_Order;k++)  
			{
			rst[1]=Gauss_Point[k][Gauss_Order-1];		  	// zeta: local coordinate 2
			double wy=Gauss_Point[k][Gauss_Order-1];         

			Compute_N();
			Compute_DN_rst();
			Compute_JM();
			Compute_DET_J();
			Compute_Inv_JM();
			Compute_DN_xyz();

			Compute_Stiffness_PartI(); // products of Ni'Nj' and W is assembled ino SMm[]

			//if ((gravx!=0) || (gravy!=0)) grav();   // gravity loads
			//if (wirl!=0) rotat();                   // rotational loads
			}
		}

	Compute_Stiffness_PartII();
	//if (maxprs!=0) pres1();
	//if (maxtns!=0) transf();
}
}


void Cquad8::Compute_Stiffness_PartI()
{
double wx,wy,DET_J;
int n=-1;
double W=wx*wy*DET_J;
double DN_xL[100],DN_yL[100];
for (int i=0;i<no_elem_nodes;i++)
	for (int j=0;j<no_elem_nodes;j++)
		{
		S[n+1]+=W*DN_xL[i]*DN_xL[j];
		S[n+2]+=W*DN_xL[i]*DN_yL[j];
		S[n+3]+=W*DN_yL[i]*DN_yL[j];
		n+=3;
		if (i!=j) {
					S[n+4]+=W*DN_yL[i]*DN_xL[j];
					n++;
				  }
		}
}


void Cquad8::Compute_Stiffness_PartII()
{
int n=-1;

for (int i=0;i<no_elem_nodes;i++)
	for (int j=0;j<no_elem_nodes;j++)
	{
// Products of derivative of N and Guass_Weigths are retrieved from elem stiffness
		double s1=S[n+1];
		double s2=S[n+2];
		double s3=S[n+3];
		if (i!=j)  // 2 x 2 sub matrix is computed
		{
			double s4=S[n+4];
			S[n+1]=s1*C[0][0]+s3*C[4][0];
			S[n+2]=s2*C[1][0]+s4*C[4][0];
			S[n+3]=s3*C[2][0]+s2*C[4][0];
			S[n+4]=s4*C[3][0]+s1*C[4][0];
			n+=4;
		}
		else 	 // Upper triangle of 2 x 2 sub-matrix is computed
		{
			S[n+1]=s1*C[0][0]+s3*C[4][0];
			S[n+2]=s2*C[2][0]+s2*C[4][0];
			S[n+3]=s3*C[3][0]+s1*C[4][0];
			n+=3;
		}
	}	
}


void Cquad8::Compute_N()
{
double xx[10],yy[10];

for (int i=0;i<no_elem_nodes;i++) 
			{
			// N and Derivatives of N wrt. to natural cordinates for each nodes
			switch (i) {
				case 0 : case 2 : case 4 : case 6 :	  	// Corner nodes
					{
					double a1=rst[0]*xx[i];
					double a2=rst[1]*yy[i];
					N[i]=0.25*(1.0+a1)*(1.0+a2)*(a1+a2-1.0);
					break;
					}

				case 1 : case 5 :						// Midside nodes rst[1]=0
					{
					N[i]=0.5*(1.0+rst[0]*xx[i])*(1.0-rst[1]*rst[1]);
					break;
					}

				case 3 : case 7 :					   	// Midside nodes rst[0]=0
					{
					N[i]=0.5*(1.0+rst[1]*yy[i])*(1.0-rst[0]*rst[0]);
					break;
					}
				}
}
}


void Cquad8::Compute_DN_rst()
{
double xx[10],yy[10];
for (int i=0;i<no_elem_nodes;i++) 
			{
			// N and Derivatives of N wrt. to natural cordinates for each nodes
			switch (i) {
				case 0 : case 2 : case 4 : case 6 :	  	// Corner nodes
					{
					double a1=rst[0]*xx[i];
					double a2=rst[1]*yy[i];
					DN_rst[i][0]=0.25*xx[i]*(1.0+a2)*(2.0*a1+a2);
					DN_rst[i][1]=0.25*yy[i]*(1.0+a1)*(2.0*a2+a1);
					break;
					}

				case 1 : case 5 :						// Midside nodes rst[1]=0
					{
 					DN_rst[i][0]=0.5*xx[i]*(1.0-rst[1]*rst[1]);
					DN_rst[i][1]=-rst[1]+(1.0+rst[0]*xx[i]);
					break;
					}
				case 3 : case 7 :					   	// Midside nodes rst[0]=0
					{
					DN_rst[i][0]=-rst[0]+(1.0+rst[1]*yy[i]);
 					DN_rst[i][1]=0.5*yy[i]*(1.0-rst[0]*rst[0]);
					break;
					}
				}
}
}


void Cquad8::Compute_JM()
{
for (int i=0;i<4;i++) JM[i][0]=0.0; // Initialize JM = [0]
for ( i=0;i<no_elem_nodes;i++) 
		{
			// Compute Jacobian Matrix			
			JM[0][0]+=DN_rst[i][0]*enode[i]->x[0];
			JM[1][0]+=DN_rst[i][0]*enode[i]->x[1];
			JM[2][0]+=DN_rst[i][1]*enode[i]->x[0];
			JM[3][0]+=DN_rst[i][1]*enode[i]->x[1];
		}
}


void Cquad8::Compute_DET_J()
{
	DET_J=JM[0][0]*JM[3][0]-JM[1][0]*JM[2][0];
	if (DET_J<=0) fprintf(outfile,"Jacobian is negetive or zero");
}

void Cquad8::Compute_Inv_JM()
{
	// Inverse of Jacobian is computed
	Inv_JM[0][0]=JM[3][0];
	Inv_JM[1][0]=-JM[1][0];
	Inv_JM[2][0]=-JM[2][0];
	Inv_JM[3][0]=JM[0][0];
	for (int i=0;i<4;i++) Inv_JM[i][0]/=DET_J;
}

void Cquad8::Compute_DN_xyz()
{
	// Derivative of N wrt x,y axes
	for (int i=0;i<no_elem_nodes;i++) 
		{
		DN_xyz[i][0]=Inv_JM[0][0]*DN_rst[0][i]+Inv_JM[1][0]*DN_rst[1][i];
		DN_xyz[i][1]=Inv_JM[2][0]*DN_rst[0][i]+Inv_JM[3][0]*DN_rst[1][i];
		}
}



