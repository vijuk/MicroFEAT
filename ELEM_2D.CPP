#include "stdafx.h"
#include "space.h"

void Ctriangle_Constant::Read_Properties()	// Added on 08/10/96.
{
Celem::Read_Properties();
fscanf(inpfile,"%lf",&Thickness);
//fprintf(outfile," %10.4lf\n",Thickness);
}

void Ctriangle_Constant::compute_R()
{
for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=0;
R[0][0]=R[1][1]=R[2][2]=1.0;
} 

void Ctriangle_Constant::Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
		double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
for (int i=0;i<elemman->first->no_elem_dof;i++)
for (int j=0;j<elemman->first->no_elem_dof;j++) SMs[i][j]=SMm[i][j];
}

void Ctriangle_Constant::Compute_SMm(double Axial_Load)
{
/*double pr=E/(2.0*G)-1.0;           // Compute poisson's ratio*/

double e1=material->E/(1.0-material->MEU*material->MEU);			// elastic coeffs
double e2=material->MEU;
double e12=e1*(1-e2)/2.0;

double x1,y1,z1,x2,y2,z2,x3,y3,z3;

nodeman->getnode(enode[0]->no,x1,y1,z1);	  // get coords
nodeman->getnode(enode[1]->no,x2,y2,z2);
nodeman->getnode(enode[2]->no,x3,y3,z3);

double bi=y2-y3;
double ci=x3-x2;

double bj=y3-y1;
double cj=x1-x3;

double bm=y1-y2;
double cm=x2-x1;

Init_Lower(SMm);

SMm[0][0]=e1*bi*bi+e12*ci*ci;

SMm[1][0]=(e1*e2+e12)*bi*ci;
SMm[1][1]=e1*ci*ci+e12*bi*bi;

SMm[2][0]=e1*bi*bj   +e12*ci*cj;
SMm[2][1]=e1*e2*bj*ci+e12*bi*cj;
SMm[2][2]=e1*bj*bj   +e12*cj*cj;

SMm[3][0]=e1*e2*bi*cj+e12*bj*ci;
SMm[3][1]=e1*ci*cj   +e12*bi*bj;
SMm[3][2]=(e1*e2+e12)*bj*cj;
SMm[3][3]=e1*cj*cj   +e12*bj*bj;

SMm[4][0]=e1*bi*bm+e12*ci*cm;
SMm[4][1]=e1*e2*bm*ci+e12*bi*cm;
SMm[4][2]=e1*bj*bm+e12*cj*cm;
SMm[4][3]=e1*e2*bm*cj+e12*bj*cm;
SMm[4][4]=e1*bm*bm+e12*cm*cm;

SMm[5][0]=e1*e2*bi*cm+e12*bm*ci;
SMm[5][1]=e1*ci*cm+e12*bi*bm;
SMm[5][2]=e1*e2*bj*cm+e12*bm*cj;
SMm[5][3]=e1*cj*cm+e12*bj*bm;
SMm[5][4]=(e1*e2+e12)*bm*cm;
SMm[5][5]=e1*cm*cm+e12*bm*bm;
double Area=(x2*y3-x3*y2-x1*(y3-y2)+y1*(x3-x2))/2.0;		 // compute area
for (int i=0;i<no_elem_dof;i++) 
	for(int j=i;j<no_elem_dof;j++)
			SMm[j][i]*=Thickness/(4.0*Area);

Fill_Upper(SMm);

Print_M(SMm);
}

void Ctriangle_Constant::compute_P()
/* function get_R() should be called before this code */
{
int MDI[MAX_NO_ELEM_DOF];

/* Compute nodal load vector */
if (is_loaded())
	{
	compute_gdir_vector(MDI);  /* on 9/9/96 */

for (int i=0;i<no_elem_dof;i++)
{
AE[MDI[i]]+=AML[i];
//fprintf(outfile,"**** m_i=%d g_i=%d  AE = %10.4lf\n",i,MDI[i],AE[MDI[i]]);
}
	}
/* End compute nodal vector */
}


