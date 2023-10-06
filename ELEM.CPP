#include "stdafx.h"
#include "space.h"

//********* Celem definition (base element)
Celem::Celem(int eID,int ned,int nen,unsigned EDMA[])
{		  
elem_ID=eID; // Element Identity Number
next=0;							   
no=0;
no_elem_dof=ned;         
no_elem_nodes=nen;
enode=new Cnode* [no_elem_nodes];
no_enode_dof=new int[no_elem_nodes];

AML=new double [no_elem_dof]; // ToDo:Allocate only if AML is present
for (int i=0;i<no_elem_dof;i++) AML[i]=0;

// Find what are the elem dof's from Element Dof Map Array for the present class.
enode_dofmap=EDMA;

for (i=0;i<no_elem_nodes;i++) {
		unsigned int tempmap=enode_dofmap[i];
		no_enode_dof[i]=0;
		do {
		    no_enode_dof[i]+=(tempmap) & 1;
		   } while ((tempmap>>=1)!=0);
		if (no_enode_dof[i]<=0) DE_MSG(1);
		}
}					 

Celem::~Celem()
{
delete enode;
delete no_enode_dof;
delete AML;
if (next) delete next;
}

void Celem::Assemble_K()
{
int gdir_for[MAX_NO_ELEM_DOF];
compute_gdir_vector(gdir_for);
	/* Generation of semibanded stiffness matrix */
	for(int m_i=0;m_i<no_elem_dof;++m_i)
		{
		int g_i=gdir_for[m_i];
		// fprintf(outfile,"mem dof = %d glob dof = %d\n",m_i,g_i);
		
		if ((g_i>-1) && (g_i<no_dof))
		for(int m_j=0;m_j<no_elem_dof;++m_j)
			{
			int g_j=gdir_for[m_j];
			if (g_j<no_dof)
				{
					if (Sol==BANDED) g_j -= g_i;    // Assemble as banded matrix
					if ((g_j<0) || (g_j>sbw-1)) continue;
					K[g_i][g_j] += SMs[m_i][m_j];
					if (dyna) M[g_i][g_j] += MMs[m_i][m_j];
				}
			}
		}
}

/*
void Cbeam3D::Assemble_K()
{
Celem::Assemble_K();
}
*/

Cbeam_3D::Cbeam_3D(int eID,int ned,int nen,unsigned EDMA[]):Celem(eID,ned,nen,EDMA)
{		  
knode=0;
PreTension=0;
}

// Check whether the element is loaded 
BOOL Celem::is_loaded()
{
for (int i=0;i<no_elem_dof;i++) 
		if (AML[i]!=0) return 1;
return 0;
}	

void Celem::Read_Properties()	// Added on 25/09/96.
{
int node_no;	
					
fscanf(inpfile,"%d",&no);
fprintf(outfile,"Element No.%3d ",no);

int MatID;	
fscanf(inpfile,"%d",&MatID);
material=materialman->getmaterial(MatID);
fprintf(outfile,"Material No. %3d : ",MatID);

for (int i=0;i<no_elem_nodes;i++) 
	{
	fscanf(inpfile,"%d",&node_no);
	if (nodeman->assign_elem_node(this,enode[i],i,node_no)==0)  // anchor enode[] pointers  
		{
    	fprintf(outfile,"Error in assigning nodes\n");
		fclose(outfile);
    	exit(1);
		}
	fprintf(outfile,"%d ",enode[i]->no); 
    }
fprintf(outfile,"\n"); 
}


void Ctruss_3D::Read_Properties()	// Modified on 25/09/96.
{
Celem::Read_Properties();
fscanf(inpfile,"%lf",&A);
fprintf(outfile," %10.4lf ",A);
}


void Cgap_3D::Read_Properties()	// Added on 21/01/97.
{
Celem::Read_Properties();
//fscanf(inpfile,"%lf",&A);
//fprintf(outfile," %10.4lf ",A);
}

void Cbeam_3D::Read_Properties()	// Modified on 25/09/96.
{
int knode_no;
Celem::Read_Properties();
fscanf(inpfile,"%lf %lf %lf %lf %d %lf",&A,&Ix,&Iy,&Iz,&knode_no,&PreTension); // PreTension added on 24/1/97
fprintf(outfile," %10.4lf %10.4lf %10.4lf %10.4lf %d %10.4lf\n",A,Ix,Iy,Iz,knode_no,PreTension);
if (knode_no) knode=nodeman->getnode(knode_no);
}

// Read element loads
void Celem::Read_Loads()
{
for (int i=0;i<no_elem_dof;i++)
			{
			fscanf(inpfile,"%lf",&AML[i]);
			fprintf(outfile,"Elem(%d) : dir=%d  AML[i]=%6.3lf\n",no,i,AML[i]);
			}
fprintf(outfile,"\n");
}

/* Computes rotation stiffness matrix */
void Cbeam_3D::compute_R()
{
double x1,x2,x3,y1,y2,y3,z1,z2,z3,cx,cy,cz,cxz,xps,yps,zps,ypv,zpv,sina,cosa;
double sq;

nodeman->getnode(enode[0]->no,x1,y1,z1);
nodeman->getnode(enode[1]->no,x2,y2,z2); 

cx=(x2-x1)/L;
cy=(y2-y1)/L;
cz=(z2-z1)/L;
cxz=sqrt(sqr(cx)+sqr(cz));

if (knode!=0) 
{
nodeman->getnode(knode->no,x3,y3,z3);
xps= x3-x1;
yps= y3-y1;
zps= z3-z1; 
}

if (eq(cxz,0)) 
	{
	R[0][1]=cy;
	R[1][0]=-cy;
	R[2][2]=1.0;
	R[0][0]=R[0][2]=R[1][1]=R[1][2]=R[2][0]=R[2][1]=0;

	if (knode!=0) 	 
		{
		sq=sqrt(sqr(xps)+sqr(zps));
		cosa=-xps*cy/sq;
		sina=zps/sq;
		R[1][0]=-cy*cosa;
		R[1][2]=sina;
		R[2][0]=cy*sina;
		R[2][2]=cosa;
		}
	} 
else
	{
	R[0][0]=cx;
	R[0][1]=cy;
	R[0][2]=cz;
	R[1][0]=-cx*cy/cxz;
	R[1][1]=cxz;
	R[1][2]=-cy*cz/cxz;
	R[2][0]=-cz/cxz;
	R[2][1]=0;
	R[2][2]=cx/cxz;





	if (knode!=0) 	
		{
		ypv=R[1][0]*xps+R[1][1]*yps+R[1][2]*zps;
		zpv=R[2][0]*xps+R[2][1]*yps+R[2][2]*zps;
		
		sq=sqrt((sqr(ypv)+sqr(zpv)));

		sina=zpv/sq;
		cosa=ypv/sq;

		R[1][0]=(-cx*cy*cosa-cz*sina)/cxz;
		R[1][1]=cxz*cosa;
		R[1][2]=(-cy*cz*cosa+cx*sina)/cxz;

		R[2][0]=(cx*cy*sina-cz*cosa)/cxz;
		R[2][1]=-cxz*sina;
		R[2][2]=(cy*cz*sina+cx*cosa)/cxz;
		}
	}

/*
fprintf(outfile,"Rotation matrix\n");
for (i=0;i<3;i++)
	{
	for (j=0;j<3;j++)
		fprintf(outfile,"%10.3lf ",R[i][j]);
	fprintf(outfile,"\n");
	}
*/
}

// Initialize lower traiangle of SMm[] to zero
void Celem::Init_Lower(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
// For each col. all elements below diag = 0			//  x  x  x
for (int i=0;i<no_elem_dof-1;i++)                		//	0  x  x
	for (int j=i;j<no_elem_dof;j++) A[j][i]=0.0;	//	0  0  x
}
													
// Generation of Upper symmetrical half of SMm[]
void Celem::Fill_Upper(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
for (int i=0;i<no_elem_dof-1;i++)
	for (int j=i+1;j<no_elem_dof;j++)
		A[i][j]=A[j][i];   //  modified on 9/9/96
}

//Print element stiffness matrix
void Celem::Print_M(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
for (int i=0;i<no_elem_dof;i++)
	{
	for (int j=0;j<no_elem_dof;j++) 
		fprintf(outfile,"%12.6lf ",A[i][j]);
	fprintf(outfile,"\n");
	}
}

void Ctruss_3D::compute_MMm()	  // Computes SMs directly
{
double E,MEU,WD,G;
material->get(E,MEU,WD,G);
double t=WD*A*L/6.0;

for (int i=0;i<no_elem_dof;i++) 
	for (int j=0;j<no_elem_dof;j++) 
			if (i==j) MMs[i][i]=2*t; else MMs[i][j]=0.0;

MMs[3][0]=MMs[4][1]=MMs[5][2]=t;
}


void Ctruss_3D::Compute_SMm(double Axial_Load)	  // Computes SMs directly
{
double E,MEU,WD,G;
material->get(E,MEU,WD,G);

double t=A*E/L;

// Compute R
double x1,y1,z1,x2,y2,z2;

nodeman->getnode(enode[0]->no,x1,y1,z1);
nodeman->getnode(enode[1]->no,x2,y2,z2); 

double cx=(x2-x1)/L;
double cy=(y2-y1)/L;
double cz=(z2-z1)/L;

double xx=cx*cx;
double xy=cx*cy;
double xz=cx*cz;

double yy=cy*cy;
double yz=cy*cz;
double zz=cz*cz;

// Compute SMs

SMs[3][3]=SMs[0][0]= t*xx;
SMs[4][4]=SMs[1][1]= t*yy;
SMs[5][5]=SMs[2][2]= t*zz;

SMs[4][3]=SMs[1][0]= t*xy;
SMs[5][3]=SMs[2][0]= t*xz;

SMs[3][0]=-t*xx;
SMs[4][0]=-t*xy;
SMs[5][0]=-t*xz;

SMs[2][1]= t*yz;
SMs[3][1]=-t*xy;
SMs[4][1]=-t*yy;
SMs[5][1]=-t*yz;

SMs[3][2]=-t*xz;
SMs[4][2]=-t*yz;
SMs[5][2]=-t*zz;

SMs[5][4]= t*yz;

for (int i=0;i<no_elem_dof-1;i++)
	for (int j=i+1;j<no_elem_dof;j++)
					SMs[i][j]=SMs[j][i];   
Print_M(SMs);
}


void Cgap_3D::Compute_SMm(double Axial_Load)	  // Computes SMs directly
{
double E,MEU,WD,G;
material->get(E,MEU,WD,G);

double t=1e-50; //A*E/L;  Initial stiffness = 0

// Compute R
double x1,y1,z1,x2,y2,z2;

nodeman->getnode(enode[0]->no,x1,y1,z1);
nodeman->getnode(enode[1]->no,x2,y2,z2); 

if (fabs(L)<1e-5) t=1e50;

double cx=(x2-x1)/L;
double cy=(y2-y1)/L;
double cz=(z2-z1)/L;

double xx=cx*cx;
double xy=cx*cy;
double xz=cx*cz;

double yy=cy*cy;
double yz=cy*cz;
double zz=cz*cz;

// Compute SMs

SMs[3][3]=SMs[0][0]= t*xx;
SMs[4][4]=SMs[1][1]= t*yy;
SMs[5][5]=SMs[2][2]= t*zz;

SMs[4][3]=SMs[1][0]= t*xy;
SMs[5][3]=SMs[2][0]= t*xz;

SMs[3][0]=-t*xx;
SMs[4][0]=-t*xy;
SMs[5][0]=-t*xz;

SMs[2][1]= t*yz;
SMs[3][1]=-t*xy;
SMs[4][1]=-t*yy;
SMs[5][1]=-t*yz;

SMs[3][2]=-t*xz;
SMs[4][2]=-t*yz;
SMs[5][2]=-t*zz;

SMs[5][4]= t*yz;

for (int i=0;i<no_elem_dof-1;i++)
	for (int j=i+1;j<no_elem_dof;j++)	SMs[i][j]=SMs[j][i];   
Print_M(SMs);
}


void Ctruss_3D::compute_P() //duplicated code
{
int MDI[MAX_NO_ELEM_DOF];

/* Compute nodal load vector */
if (is_loaded())
	{
	compute_gdir_vector(MDI);  
	for (int i=0;i<no_elem_dof;i++)
	{
	AE[MDI[i]]+=AML[i];
	//fprintf(outfile,"**** m_i=%d g_i=%d  AE = %10.4lf\n",i,MDI[i],AE[MDI[i]]);
	}
	}
/* End compute nodal vector */
}


// Transformation applied to element stiffness and mass matrices
// **** Element dof must be 12 for using this routine.
void Cbeam_3D::Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
		double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
int j,k,k1,k2,k3,j1,j2,j3;
double SMmR[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];

for (k=1;k<=4;k++)                 /* index for blocks of columns in SMm */
	{
	k1=3*k-3;                  /* determine column indices */
	k2=3*k-2;
	k3=3*k-1;
	for (j=0;j<NEDOF_Cbeam_3D;j++)        /* index for rows of SMm */
		{
SMmR[j][k1]=SMm[j][k1]*R[0][0]+SMm[j][k2]*R[1][0]+SMm[j][k3]*R[2][0];
SMmR[j][k2]=SMm[j][k1]*R[0][1]+SMm[j][k2]*R[1][1]+SMm[j][k3]*R[2][1];
SMmR[j][k3]=SMm[j][k1]*R[0][2]+SMm[j][k2]*R[1][2]+SMm[j][k3]*R[2][2];
		}
	}                         /* SMm * R */

for (j=1;j<=4;j++)
	{
	j1=3*j-3;
	j2=3*j-2;
	j3=3*j-1;
	for (k=j1;k<NEDOF_Cbeam_3D;k++)
		{
SMs[j1][k]=R[0][0]*SMmR[j1][k]+R[1][0]*SMmR[j2][k]+R[2][0]*SMmR[j3][k];
SMs[j2][k]=R[0][1]*SMmR[j1][k]+R[1][1]*SMmR[j2][k]+R[2][1]*SMmR[j3][k];
SMs[j3][k]=R[0][2]*SMmR[j1][k]+R[1][2]*SMmR[j2][k]+R[2][2]*SMmR[j3][k];
		}
	}                        /* R' * SMm * R */

/*
for(i=1;i<=12;i++)
	for(k=1;k<=12;k++)
		FEAs[i]+=R[k][i]*FEAm[mno][k];
*/
for (int i=0;i<no_elem_dof-1;i++)
	for (int j=i+1;j<no_elem_dof;j++)
					SMs[j][i]=SMs[i][j];   // Not required  when assembling. 
											// Used in NLinear.cpp when finding internal forces
											// avoid later
}


// Transformation applied to element forces
void Cbeam_3D::Change_2_Global(double *X)
{
int j,j1,j2,j3;
double Xs[MAX_NO_ELEM_DOF];
for (j=1;j<=4;j++)
	{
	j1=3*j-3;
	j2=3*j-2;
	j3=3*j-1;
Xs[j1]=R[0][0]*X[j1]+R[1][0]*X[j2]+R[2][0]*X[j3];
Xs[j2]=R[0][1]*X[j1]+R[1][1]*X[j2]+R[2][1]*X[j3];
Xs[j3]=R[0][2]*X[j1]+R[1][2]*X[j2]+R[2][2]*X[j3];
	}                        /* R' * SMm * R */
for(int i=0;i<12;i++) X[i]=Xs[i];
}

// Transformation applied to element forces
void Cbeam_3D::Change_2_Local(double *X)
{
int j,j1,j2,j3;
double Xm[MAX_NO_ELEM_DOF];
for (j=1;j<=4;j++)
	{
	j1=3*j-3;
	j2=3*j-2;
	j3=3*j-1;
Xm[j1]=R[0][0]*X[j1]+R[0][1]*X[j2]+R[0][2]*X[j3];
Xm[j2]=R[1][0]*X[j1]+R[1][1]*X[j2]+R[1][2]*X[j3];
Xm[j3]=R[2][0]*X[j1]+R[2][1]*X[j2]+R[2][2]*X[j3];
	}                        /* R' * SMm * R */
for(int i=0;i<12;i++) X[i]=Xm[i];
}

// map local directions to global directions & verify no_elem_dof.
void Celem::compute_gdir_vector(int gdir_for[])
{
// modified on 9.9.96						 
int Tno_elem_dof=0;
int dir;
for(int i=0;i<no_elem_nodes;i++)
	{
	for(dir=0; dir<enode[i]->no_node_dof; dir++)
				gdir_for[Tno_elem_dof+dir]=enode[i]->ID[dir];
	Tno_elem_dof+=enode[i]->no_node_dof;
	}
if (no_elem_dof!=Tno_elem_dof) disp_error(1);
}

int Celem::get_elem_sbw()
{
// modified on 9.9.96
int max=0;
int min=no_dof;
for (int j=0;j<no_elem_nodes;j++)
for (int i=0;i<enode[j]->no_node_dof;i++)
		{
		if (enode[j]->ID[i]>max) 
								max=enode[j]->ID[i];
		if (enode[j]->ID[i]>-1 && enode[j]->ID[i]<min) 
							min=enode[j]->ID[i];
		}
return (max-min+1);
}

void Cbeam_3D::compute_P()
/* function get_R() should be called before this code */
{
int MDI[MAX_NO_ELEM_DOF];
int j,i1,i2,i3,j1,j2,j3;
/* Compute nodal load vector */
if (is_loaded())
	{
	compute_gdir_vector(MDI);  /* on 9/9/96 */
	for (j=1;j<=4;j++) /* index for blocks of terms in AML */
		{				   
		j1=3*j-3;   /* indexes for individual terms in AML & AE */
		j2=3*j-2;
		j3=3*j-1;
		i1=MDI[j1];
		i2=MDI[j2];
		i3=MDI[j3];

/* AE = -R'*AML */
if (i1>=0) AE[i1]+=-R[0][0]*AML[j1]-R[1][0]*AML[j2]-R[2][0]*AML[j3];
if (i2>=0) AE[i2]+=-R[0][1]*AML[j1]-R[1][1]*AML[j2]-R[2][1]*AML[j3];
if (i3>=0) AE[i3]+=-R[0][2]*AML[j1]-R[1][2]*AML[j2]-R[2][2]*AML[j3];
		}
	}
/* End compute nodal vector */
}

void Cbeam_2D_VArch::compute_P()
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

/********** stiffness for straight beam element */

void Cbeam_3D::Compute_SMm(double Axial_Load)
{
double EAbyL,EIzbyL2,EIybyL2,EIybyL,EIzbyL,GJbyL;

double E,MEU,WD,G;
material->get(E,MEU,WD,G);

EAbyL  =E/L*A;
EIzbyL2=E*Iz/sqr(L);
EIybyL2=E*Iy/sqr(L);
EIybyL =E*Iy/L;
EIzbyL =E*Iz/L;
GJbyL  =G*Ix/L;

Init_Lower(SMm);

// Populate lower triangle in SMm

	SMm[0][0]   =EAbyL;
	SMm[6][0]   =-EAbyL;
	SMm[1][1]   =12*EIzbyL2/L;
	SMm[5][1]   =6*EIzbyL2;
	SMm[7][1]   =-SMm[1][1];
	SMm[11][1]  =SMm[5][1];

	SMm[2][2]   =12*EIybyL2/L;
	SMm[4][2]   =-6*EIybyL2;
	SMm[8][2]   =-SMm[2][2];
	SMm[10][2]  =SMm[4][2];

	SMm[3][3]   =GJbyL;
	SMm[9][3]  =-SMm[3][3];

	SMm[4][4]   =4*EIybyL;
	SMm[8][4]   =-SMm[10][2];
	SMm[10][4]  =SMm[4][4]/2;

	SMm[5][5]   =4*EIzbyL;
	SMm[7][5]   =-SMm[5][1];
	SMm[11][5]  =SMm[5][5]/2;

	SMm[6][6]   =SMm[0][0];
	SMm[7][7]   =SMm[1][1];
	SMm[11][7]  =SMm[7][5];
	SMm[8][8]   =SMm[2][2];
	SMm[10][8]  =SMm[8][4];
	SMm[9][9]   =SMm[3][3];
	SMm[10][10] =SMm[4][4];
	SMm[11][11] =SMm[5][5];

if (!eq(Axial_Load+PreTension,0)) Include_SMmP(Axial_Load+PreTension);  // Include only on user request.
Fill_Upper(SMm);
// Print_M(SMm);
}

void Cbeam_3D::Include_SMmP(double Axial_Load)  //ToDo(vks): Verify all the coefficients in the matrix
{
Init_Lower(SMmP);


 double  a36=36.0/(L*L);
 double  a3 = 3.0/L;

  SMmP[1][1]=  a36;
  SMmP[5][1]=  a3;
  SMmP[7][1]= -a36;
 SMmP[11][1]=  a3;

  SMmP[2][2]=  a36;
  SMmP[4][2]= -a3; //
  SMmP[8][2]= -a36;
 SMmP[10][2]= -a3; //

  //
 
  SMmP[4][4]=   4;
  SMmP[8][4]=  a3; //
 SMmP[10][4]=  -1;

  SMmP[5][5]=   4;
  SMmP[7][5]=   -a3;
  SMmP[11][5]=  -1;
 
  //

  SMmP[7][7]=   a36;
 SMmP[11][7]=  -a3;

  SMmP[8][8]=   a36;
 SMmP[10][8]=   a3; //
 
 //
 SMmP[10][10]=   4;
 SMmP[11][11]=   4;
 
 //
//  SMmP by J.L.MEEK ends


/*

SMmP[1][4]=SMmP[1][5]=SMmP[1][10]=SMmP[1][11]=3L;
SMmP[1][7]=SMmP[1][8]=-(SMmP[1][1]=SMmP[1][2]=36);

SMmP[4][4]=SMmP[4][5]=4*L*L;
SMmP[4][10]=SMmP[4][11]=-L*L;
SMmP[4][7]=SMmP[4][8]=(SMmP[4][1]=SMmP[4][2]=3*L);

SMmP[7][1]=SMmP[7][2]=-36*L;
SMmP[7][4]=SMmP[7][5]=4*L*L;						
SMmP[7][7]=SMmP[7][8]=SMmP[7][10]=SMmP[7][11]=-3*L;

SMmP[10][1]=SMmP[10][2]=3*L;
SMmP[10][4]=SMmP[10][5]=SMmP[10][10]=SMmP[10][11]=-3*L;
SMmP[10][7]=SMmP[10][8]=36.0;

for(int i=0;i<no_elem_dof;i++) {
					SMmP[2][i]=SMmP[1][i];
					SMmP[5][i]=SMmP[5][4];
					SMmP[8][i]=SMmP[7][i];
					SMmP[11][i]=SMmP[10][i];
				  }
*/
/*  u1  v1	w1	  tu1	tv1  tw1  u2    v2   w2    tu2   tv2    tw2
    1  2     3     4     5    6   7     8     9     10    11    12
u1  0  0	 0	   0	 0	  0   0     0     0     0     0      0
v1  0  36    36	   0    3L	 3L	  0     -36  -36    0     3L    3L
w1  0  36    36	   0    3L	 3L	  0     -36  -36    0     3L    3L
tu1 0  0      0    0     0    0   0     0     0     0     0      0 
tv1 0  3L    3L    0   4LL	4LL   0		-3L	 -3L	0	 -LL   -LL
tw1 0  3L    3L    0   4LL	4LL   0		-3L	 -3L	0	 -LL   -LL
u2  0  0	 0	   0	 0	  0   0      0     0    0     0      0
v2  0 -36L  -36L   0   4LL	4LL   0     -3L  -3L    0    -3L   -3L				  					
w2  0 -36L  -36L   0   4LL	4LL   0     -3L  -3L    0    -3L   -3L				  					
tu2 0  0      0    0     0    0   0     0     0     0     0      0 
tv2 0  3L    3L    0   -3L	-3L   0		36	  36  	0	 -3L   -3L
tw2 0  3L    3L    0   -3L	-3L   0		36	  36  	0	 -3L   -3L
*/

// Used  only while debugging
/*
for (int i=0;i<no_elem_dof-1;i++)                	
	for (int j=i+1;j<no_elem_dof;j++) if (!eq(SMmP[j][i],SMmP[i][j])) RE_MSG(4);
*/



/*     SMmP[1][1]=36;
       SMmP[5][1]=3*L;
       SMmP[7][1]=-36;
       SMmP[11][1]=3*L;
       SMmP[5][5]=4*(L*L);
       SMmP[7][5]=-3*L;
       SMmP[11][5]=-1*(L*L);
       SMmP[7][7]=36;
       SMmP[11][7]=-3*L;
       SMmP[11][11]=4*(L*L);
*/

       for (int i=0;i<12;i++)
       for (int j=0;j<12;j++)
       		SMmP[i][j]=(Axial_Load/(30*L))*SMmP[i][j]        *L*L;//Meek only

Fill_Upper(SMmP); 

/*
fprintf(outfile," GEOMETRIC STIFFNESS MATRIX\n");
       for (i=0;i<12;i++)
	   {
       for (int j=0;j<12;j++)
	       fprintf(outfile," %15.0lf",SMmP[i][j]);
	   fprintf(outfile,"\n");
	   }

*/
/*
fprintf(outfile,"SMmP\n");
Print_M(SMmP);
*/

for ( i=0;i<no_elem_dof;i++)                	
	for (int j=i;j<no_elem_dof;j++) SMm[j][i]+=SMmP[j][i];
}

void Cbeam_3D::compute_MMm()
{
double	r=0; /* #change */
MMm[0][0]=MMm[6][6]                          = 140;
MMm[1][1]=MMm[2][2]=MMm[7][7]=MMm[8][8]      = 156;
MMm[3][3]=MMm[9][9]                          = 140*r*r;
MMm[4][4]=MMm[5][5]=MMm[10][10]=MMm[11][11]  = 4*L*L;
MMm[4][2]=MMm[11][7]=-(MMm[5][1]=MMm[10][8]  = 22.0*L);
MMm[8][4]=MMm[11][1]=-(MMm[10][2]=MMm[7][5]  = 13*L);
int  	i2,j2;
double	temp=material->WD*A*L/420.0;

for(i2=0;i2<4;++i2)
for(j2=i2+1;j2<4;++j2)
	MMm[j2][i2]=(MMm[i2][j2]=temp*MMm[i2][j2]);
}                      

void Cbeam_3D_curved::Read_Properties()
{
int knode_no;
Celem::Read_Properties();
fscanf(inpfile,"%lf %lf %lf %lf %d %lf %lf",&A,&Ix,&Iy,&Iz,&knode_no,&r,&IncAng); // PreTension added on 24/1/97
fprintf(outfile," %10.4lf %10.4lf %10.4lf %10.4lf %d r = %lf IA = %lf \n",A,Ix,Iy,Iz,knode_no,r,IncAng);
if (knode_no) knode=nodeman->getnode(knode_no);
}

void Cbeam_3D_curved::Compute_SMm(double Axial_Load)
{
double E,MEU,WD,G;
material->get(E,MEU,WD,G);

double d=PI*IncAng/180.0; 
double J=Ix;
double sind=sin(d);
double cosd=1-cos(d);
double sin2d=sin(2*d);

double t0=r*r*r/(E*Iz);
double t1=r/(A*E)+t0;
double t2=d/2-sin2d/4;
double t3=d/2+sin2d/4;
double t4=r/(A*G);

double p1=t1*t2+t3*t4;
double p2=-(t1-t4)*sqr(sind)/2+t0*cosd;
double p3=t0/r*cosd;
double p4=t1*t3+t0*(d-2*sind)+t4*t2;
double p5=t0/r*(d-sind);
double p6=t0/(r*r)*d;

double D=p1*(p4*p6-p5*p5)-p2*(p2*p6-p3*p5)+p3*(p2*p5-p3*p4);

//Populate Lower triangle in SMm()

Init_Lower(SMm);

	SMm[0][0]   =(p4*p6-p5*p5)/D;
	SMm[1][0]   =(p3*p5-p2*p6)/D;
	SMm[5][0]   =(p2*p5-p3*p4)/D;

	SMm[1][1]   =(p1*p6-p3*p3)/D;
 	SMm[5][1]   =(p2*p3-p1*p5)/D;
	SMm[5][5]   =(p1*p4-p2*p2)/D;

	SMm[6][0]   =-SMm[0][0];
	SMm[6][1]   =-SMm[1][0];
	SMm[6][5]   =-SMm[5][0];
	SMm[7][6]   =-SMm[6][1];
	SMm[6][6]   =-SMm[6][0];
	SMm[7][0]   =-SMm[1][0];
	SMm[7][1]   =-SMm[1][1];
	SMm[7][7]   =-SMm[7][1];
	SMm[7][5]   =-SMm[5][1];  

	SMm[11][7]  =-(SMm[11][1]= -SMm[1][0]*r*sind- SMm[1][1]*r*cosd- SMm[5][1]); 
	SMm[11][5]  =              -SMm[5][0]*r*sind- SMm[5][1]*r*cosd- SMm[5][5];
	SMm[11][6]  =-(SMm[11][0]= -SMm[0][0]*r*sind- SMm[1][0]*r*cosd- SMm[5][0]);
	SMm[11][11] =             -SMm[11][0]*r*sind-SMm[11][1]*r*cosd-SMm[11][5];

	t1=r*r/(G*J);
	t2=r*r/(E*Iy);
	t3=t1*r;
	t4=t2*r;
	t0=t1/r;
	double t5=t2/r;

	p1= (t5+t0)*d/2+(t5-t0)*sin2d/4;
	p2= (t5-t0)*sqr(sind)/2;
	p3= (t1-t2)*sqr(sind)/2+t1*-cosd;
	p4= (t0+t5)*d/2+(t2-t1)*sin2d/4+t1*sind;
	p5=-(t1+t2)*d/2+(t2-t1)*sin2d/4+t1*sind;     
	p6= (t3+t4)*d/2+(t3-t4)*sin2d/4+t3*(d-2*sind)+r*d/(A*G);

	D=p1*(p4*p6-p5*p5)-p2*(p2*p6-p3*p5)+p3*(p2*p5-p3*p4);

	SMm[2][2]=(p1*p4-p2*p2)/D;

	SMm[3][2]=(p2*p5-p3*p4)/D;
	SMm[3][3]=(p4*p6-p5*p5)/D;

	SMm[4][2]=(p2*p3-p1*p5)/D;
	SMm[4][3]=(p3*p5-p2*p6)/D;
	SMm[4][4]=(p1*p6-p3*p3)/D;

	SMm[8][2]=-SMm[3][2]*r*sind-SMm[4][2]*r*cosd-SMm[2][2];
	SMm[8][3]=-SMm[3][3]*r*sind-SMm[4][3]*r*cosd-SMm[3][2];
	SMm[8][4]=-SMm[4][3]*r*sind-SMm[4][4]*r*cosd-SMm[4][2];
	SMm[8][8]=-SMm[8][3]*r*sind+SMm[4][3]*r*cosd-SMm[8][2];

	SMm[9][2]=-SMm[3][2];
	SMm[9][3]=-SMm[3][3];
	SMm[9][4]=-SMm[4][3];
	SMm[9][9]= SMm[3][3];
	SMm[9][8]=-SMm[9][3]*r*sind-SMm[9][4]*r*cosd-SMm[9][2];

	SMm[10][2]=-SMm[4][2];
	SMm[10][3]=-SMm[4][3];
	SMm[10][4]=-SMm[4][4];
	SMm[10][8]=-SMm[10][3]*r*sind-SMm[10][4]*r*cosd-SMm[10][2];
	SMm[10][9]= SMm[4][3];
	SMm[10][10]=SMm[4][4];

Fill_Upper(SMm);
/*
fprintf(outfile,"Cbeam_3D_Curved No.%d : SMm\n",no);
Print_M(SMm);
*/
}

void Cbeam_3D_curved::Local_2_Global(
			double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
		 	double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
Cbeam_3D::Local_2_Global(SMm,SMs);
}

void Cbeam_2D_VArch::Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
		double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF])
{
for (int i=0;i<no_elem_dof;i++)
for (int j=0;j<no_elem_dof;j++) SMs[i][j]=SMm[i][j];
}




/* Computes rotation stiffness matrix */
void Cbeam_3D_curved::compute_R()
{
double x1,x2,x3,y1,y2,y3,z1,z2,z3,cx,cy,cz,cxz,xps,yps,zps,ypv,zpv,sina,cosa;
double sq;

//nodeman.getnode(enode[0]->no,x1,y1,z1);

x1=0;y1=0;z1=0;   // origin
nodeman->getnode(enode[0]->no,x2,y2,z2);  // enode[1] == enode[0]

knode=enode[1];
L=sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2));

cx=(x2-x1)/L;
cy=(y2-y1)/L;
cz=(z2-z1)/L;
cxz=sqrt(sqr(cx)+sqr(cz));

if (knode!=0) 
{
nodeman->getnode(knode->no,x3,y3,z3);  // enode[1] == knode 
xps= x3-x1;
yps= y3-y1;
zps= z3-z1; 
}

if (eq(cxz,0)) 
	{
	R[0][1]=cy;
	R[1][0]=-cy;
	R[2][2]=1.0;
	R[0][0]=R[0][2]=R[1][1]=R[1][2]=R[2][0]=R[2][1]=0;

	if (knode!=0) 	 
		{
		sq=sqrt(sqr(xps)+sqr(zps));
		cosa=-xps*cy/sq;
		sina=zps/sq;
		R[1][0]=-cy*cosa;
		R[1][2]=sina;
		R[2][0]=cy*sina;
		R[2][2]=cosa;
		}
	} 
else
	{
	R[0][0]=cx;
	R[0][1]=cy;
	R[0][2]=cz;
	R[1][0]=-cx*cy/cxz;
	R[1][1]=cxz;
	R[1][2]=-cy*cz/cxz;
	R[2][0]=-cz/cxz;
	R[2][1]=0;
	R[2][2]=cx/cxz;

	if (knode!=0) 	
		{
		ypv=R[1][0]*xps+R[1][1]*yps+R[1][2]*zps;
		zpv=R[2][0]*xps+R[2][1]*yps+R[2][2]*zps;
		
		sq=sqrt((sqr(ypv)+sqr(zpv)));

		sina=zpv/sq;
		cosa=ypv/sq;

		R[1][0]=(-cx*cy*cosa-cz*sina)/cxz;
		R[1][1]=cxz*cosa;
		R[1][2]=(-cy*cz*cosa+cx*sina)/cxz;

		R[2][0]=(cx*cy*sina-cz*cosa)/cxz;
		R[2][1]=-cxz*sina;
		R[2][2]=(cy*cz*sina+cx*cosa)/cxz;
		}
	}

/*
fprintf(outfile,"Rotation matrix\n");
for (i=0;i<3;i++)
	{
	for (j=0;j<3;j++)
		fprintf(outfile,"%10.3lf ",R[i][j]);
	fprintf(outfile,"\n");
	}
*/
}

void Cbeam_2D_VArch::compute_R()
{
for (int i=0;i<3;i++) for (int j=0;j<3;j++) R[i][j]=0;
R[0][0]=R[1][1]=R[2][2]=1.0;
}


void Cbeam_2D_VArch::Compute_SMm(double Axial_Load)
{
// Ref: FE Struct. Analysis by T.Y.Yang Page 190.
// Cross checked on 17/9/96. 

double E,MEU,WD,G;
material->get(E,MEU,WD,G);
 
double b=IncAng*PI/180.0;
double a=Ix/A;// ********VERIFY !!! : 22/9/1996.
L=b*r;

Init_Lower(SMm);

//Kuu
double p=E*A/5*(1+a/(r*r));
double q=E*A/r;
double s=E*A/L;


	SMm[0][0]   = 6/L  *p;
	SMm[1][0]   = 1/2.0*p;
    SMm[2][0]   =-1/2.0*q;
    SMm[3][0]   =(L/10+a/L)*q;
	SMm[4][0]   =-6/L   *p;
	SMm[5][0]   = 1/2.0 *p;
	SMm[6][0]   = 1/2.0*q;
	SMm[7][0]   =-(L/10+a/L)*q;


	SMm[1][1]   = 2*L/3 *p;
    SMm[2][1]   =-(L/10+a/L)*q;
    SMm[3][1]   = a/2*q;
	SMm[4][1]   =-1/2.0 *p;
	SMm[5][1]   =-L/6   *p;
	SMm[6][1]   = (L/10+a/L)*q;
	SMm[7][1]   =-(L*L/60+a/2)*q;


  	SMm[2][2]   =(13*b*b/35+12*a/(L*L))*s;
	SMm[3][2]   =(11*L*b*b/210+6*a/L)*s;
	SMm[4][2]   = -1/2.0*q;
	SMm[5][2]   =(L/10+a/L)*q;
	SMm[6][2]   =(9*b*b/70-12*a/(L*L))*s;
	SMm[7][2]   =(-13*L*b*b/420+6*a/L)*s;
	

	SMm[3][3]   =(L*L*b*b/105+4*a)*s;
	SMm[4][3]   = -(L/10+a/L)*q;
	SMm[5][3]   =(L*L/60+a/2)*q;
	SMm[6][3]   =(13*L*b*b/420-6*a/L)*s;		        // corrected on 18/9/96
	SMm[7][3]   =(-L*L*b*b/140+2*a)*s;
	

	SMm[4][4]   = 6/L   *p;
	SMm[5][4]   =-1/2.0   *p;
	SMm[6][4]   = 1/2.0*q;
	SMm[7][4]   = (L/10+a/L)*q;

	SMm[5][5]   = 2*L/3.0 *p;
	SMm[6][5]   =-(L/10+a/L)*q;
	SMm[7][5]   =-a/2*q;

	SMm[6][6]   =(13*b*b/35+12*a/(L*L))*s;
	SMm[7][6]   =-(11*L*b*b/210+6*a/L)*s;

	SMm[7][7]   =(L*L*b*b/105+4*a)*s;


/*
	  0                1               4             5                         
0   6/L*p 		     1/2*p		      -6/L*p		1/2*p	   -1/2*q  (L/10+a/L)*q  1/2*q   -(L/10+a/L)*q
1   1/2*p            2*L/3*p	     -1/2*p         -L/6*p     -(L/10+a/L)*q   a/2*q   (L/10+a/L)*q  -(L*L/60+a/2)*q
4  -6/L*p           -1/2*p           6/L*p         -1/2*p	   -1/2*q -(L/10+a/L)*q  1/2*q  (L/10+a/L)*q
5   1/2*p           -L/6*p          -1/2*p          2*L/3*p     (L/10+a/L)*q  (L*L/60+a/2)*q  -(L/10+a/L)*q  -a/2*q  

2  -1/2*q	          -(L/10+a/L)*q   -1/2*q          (L/10+a/L)*q   (13*b*b/35+12*a/(L*L))*s (11*L*b*b/210+6*a/L)*s (9*b*b/70-12*a/(L*L))*s  (-13*L*b*b/420+6*a/L)*s
3   (L/10+a/L)*q        a/2*q         -(L/10+a/L)*q   (L*L/60+a/2)*q (11*L*b*b/210+6*a/L)*s    (L*L*b*b/105+4*a)*s  13*L*b*b/420-6*a/L)*s   (-L*L*b*b/140+2*a)*s
6   1/2*q;		       (L/10+a/L)*q    1/2*q         -(L/10+a/L)*q  (9*b*b/70-12*a/(L*L))*s   (13*L*b*b/420-6*a/L)*s (13*b*b/35+12*a/(L*L))*s -(11*L*b*b/210+6*a/L)*s
7  -(L/10+a/L)*q	  -(L*L/60+a/2)*q  (L/10+a/L)*q   -a/2*q          (-13*L*b*b/420+6*a/L)*s   (-L*L*b*b/140+2*a)*s    -(11*L*b*b/210+6*a/L)*s  (L*L*b*b/105+4*a)*s
	 2				   3               6                 7
0  -1/2*q	          -(L/10+a/L)*q   -1/2*q          (L/10+a/L)*q   (13*b*b/35+12*a/(L*L))*s (11*L*b*b/210+6*a/L)*s (9*b*b/70-12*a/(L*L))*s  (-13*L*b*b/420+6*a/L)*s
1  (L/10+a/L)*q        a/2*q         -(L/10+a/L)*q   (L*L/60+a/2)*q (11*L*b*b/210+6*a/L)*s    (L*L*b*b/105+4*a)*s  13*L*b*b/420-6*a/L)*s   (-L*L*b*b/140+2*a)*s
4   1/2*q;		       (L/10+a/L)*q    1/2*q         -(L/10+a/L)*q  (9*b*b/70-12*a/(L*L))*s   (13*L*b*b/420-6*a/L)*s (13*b*b/35+12*a/(L*L))*s -(11*L*b*b/210+6*a/L)*s
5  -(L/10+a/L)*q	  -(L*L/60+a/2)*q  (L/10+a/L)*q   -a/2*q          (-13*L*b*b/420+6*a/L)*s   (-L*L*b*b/140+2*a)*s    -(11*L*b*b/210+6*a/L)*s  (L*L*b*b/105+4*a)*s
*/


Fill_Upper(SMm);
//fprintf(outfile,"Element stiffness matrix in local coordinates \n");
//Print_M(SMm);
}

void Cbeam_2D_HArch::Compute_SMm(double Axial_Load)
{
double SMmb[12][12],SMmt[12][12],SMmw[12][12];
int i,j;

for (i=0;i<no_elem_dof;i++)  // Operation done on diagonal elements also.
	for (j=i;j<no_elem_dof;j++) 
		SMm[j][i]=SMmb[j][i]=SMmt[j][i]=SMmw[j][i]=0.0;	

 double L=r*IncAng/180*PI;
 double L2=L*L;
 double L4=L2*L2;
 double r2=r*r;
 double rL=r*L;
 double r2L2=r2*L2;
 double rL2=r*L2;
 double r2L=r2*L;
 double L3=L2*L;

// Stiffness matrix Kt from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402. SMmb[0][0]=  12./L4;

 SMmb[0][0]=  12./L4; 
 SMmb[1][0]=  6./L3; 
 SMmb[2][0]=  6./(5*rL2);
 SMmb[3][0]=  1./(10*rL);
 SMmb[4][0]= -12/L4;
 SMmb[5][0]=  6/L3;
 SMmb[6][0]=  -6./(5*rL2);
 SMmb[7][0]=  1./(10*rL);

 SMmb[1][1]=  4/L2;
 SMmb[2][1]=  11./(10*rL);
 SMmb[3][1]= 2./(15*r);
 SMmb[4][1]= -6/L3;
 SMmb[5][1]= 2/L2;
 SMmb[6][1]= -1./(10*rL);
 SMmb[7][1]= -1./(30*r);

 SMmb[2][2]= 13./(35*r2);
 SMmb[3][2]= 11.*L/(210*r2);
 SMmb[4][2]= -6./(5*rL2);
 SMmb[5][2]= 1./(10*rL);
 SMmb[6][2]= 9./(70*r2); 			
 SMmb[7][2]= -13.*L/(420*r2); 

 SMmb[3][3]= L2/(105*r2);
 SMmb[4][3]= -1./(10*rL);
 SMmb[5][3]= -1./(30*r);
 SMmb[6][3]= 13*L/(420*r2);
 SMmb[7][3]= -L2/(140*r2);

 SMmb[4][4]=  12/L4;   // Value given in journal is negetive 12/L4 !!!!!!!!!!!!!!
 SMmb[5][4]= -6./L3;
 SMmb[6][4]=  6./(5*rL2);		
 SMmb[7][4]=  -1./(10*rL);

 SMmb[5][5]=  4/L2;
 SMmb[6][5]= -11./(10*rL);
 SMmb[7][5]=  2./(15*r);

 SMmb[6][6]=   13./(35*r2);
 SMmb[7][6]=   -11./(210*r2);

 SMmb[7][7]=   L2/(105*r2);

// Stiffness matrix Kt from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402.

  SMmt[0][0]=   6./(5*r2L2); 
  SMmt[1][0]=   1./(10*r2L); 
  SMmt[2][0]=   6./(5*rL2);
  SMmt[3][0]=   1./(10*rL); 
  SMmt[4][0]=   6./(5*r2L2);
  SMmt[5][0]=   1./(10*r2L);
  SMmt[6][0]=  -6./(5*rL2);
  SMmt[7][0]=   1./(10*rL);

  SMmt[1][1]=    2./(15*r2);
  SMmt[2][1]=    1./(10*rL);
  SMmt[3][1]=    2./(15*r);
  SMmt[4][1]=   -1./(10*r2L);
  SMmt[5][1]=   -1./(20*r2);
  SMmt[6][1]=   -1./(10*rL);
  SMmt[7][1]=   -1./(30*r);

  SMmt[2][2]=    6./(5*L2);
  SMmt[3][2]=    1./(10*L);
  SMmt[4][2]=   -6./(5*rL2);
  SMmt[5][2]=    1./(10*rL);
  SMmt[6][2]=   -6./(5*L2);
  SMmt[7][2]=   1./(10*L);

  SMmt[3][3]=  2./15;
  SMmt[4][3]=  1./(10*rL);
  SMmt[5][3]=  -1./(10*rL);
  SMmt[6][3]=  -1./(10*L);
  SMmt[7][3]=  -1./30; 

  SMmt[4][4]=   6./(5*r2L2); 
  SMmt[5][4]=  -1./(10*r2L); 
  SMmt[6][4]=   6./(5*rL2) ;
  SMmt[7][4]=   -1./(10*rL);
 
  SMmt[5][5]=   2./(15*r2);
  SMmt[6][5]=   -1./(10*rL);
  SMmt[7][5]=   2./(15*r);

  SMmt[6][6]=   6./(5*L2);
  SMmt[7][6]=   1./(10*L);

  SMmt[7][7]=  2./15;

// Stiffness matrix Kw from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402.

double r2L3=r2L2*L;
double r2L4=r2L2*L2;
double rL3=r*L3;
double rL4=r*L4;


  SMmw[0][0]=   12./(r2L4);
  SMmw[1][0]=  	6./(r2L3);
  SMmw[2][0]=  	12./(rL4);
  SMmw[3][0]=  	6./(rL3);
  SMmw[4][0]=  	-12/(r2L4);
  SMmw[5][0]=  	6/(r2L3);
  SMmw[6][0]=  	-12/(rL4);
  SMmw[7][0]=  	6/(rL3);

  SMmw[1][1]=  -4/r2L2;
  SMmw[2][1]=  6/(rL3);
  SMmw[3][1]=  4/(rL2);
  SMmw[4][1]=  -6/r2L3;  
  SMmw[5][1]=   2/r2L2; 
  SMmw[6][1]=   -6/(rL3); 
  SMmw[7][1]=   2./(rL2); 
		  
  SMmw[2][2]= 12/L4;
  SMmw[3][2]=  6/L3;
  SMmw[4][2]=  -12/rL4;
  SMmw[5][2]=  6/rL3;
  SMmw[6][2]=  -12/L4;
  SMmw[7][2]=  6/L3;
  
  SMmw[3][3]=  4/L2;
  SMmw[4][3]=  -6/rL3;
  SMmw[5][3]=  2./rL2;
  SMmw[6][3]=  -6/L3;
  SMmw[7][3]=  2/L2;

  SMmw[4][4]=   12/r2L4;
  SMmw[5][4]=   -6/r2L3;
  SMmw[6][4]=   12/rL4;
  SMmw[7][4]=  -6/rL3;
  	   
  SMmw[5][5]=  4/r2L2; 
  SMmw[6][5]=   -6/rL3;
  SMmw[7][5]=  4/rL2; 

  SMmw[6][6]=  12/L4; 
  SMmw[7][6]=  -6/L3;
 
  SMmw[7][7]=  4/L2;


// Conventions violated ; Write a note on this in Notes.FEM

double J=Ix;
double I=Iy;
double Iw=Iz;

double E,MEU,WD,G;
material->get(E,MEU,WD,G);

double EIL=E*I*L;
double GJL=G*J*L;
double EIwL=E*Iw*L;

for (i=0;i<no_elem_dof;i++)                	
	for (j=i;j<no_elem_dof;j++) 
// Operation done on diagonal elements also.
		SMm[j][i]=(SMmb[j][i]=(EIL)*SMmb[j][i])+
				  (SMmt[j][i]=(GJL)*SMmt[j][i])+
				  (SMmw[j][i]=(EIwL)*SMmw[j][i]);	

Fill_Upper(SMm);
//fprintf(outfile,"Element stiffness matrix in local coordinates \n");
//Print_M(SMm);
}

void Cbeam_2D_HArch_12::Compute_SMm(double Axial_Load)
{
double C[12][12],SMmb[12][12],SMmt[12][12],SMmw[12][12];
int i,j,k;

double L=r*IncAng/180*PI;
 double L2=L*L;
 double L4=L2*L2;
 double r2=r*r;
 double rL=r*L;
 double r2L2=r2*L2;
 double rL2=r*L2;
 double r2L=r2*L;
 double L3=L2*L;

for (i=0;i<no_elem_dof;i++)                	
	for (j=0;j<no_elem_dof;j++) 
		C[i][j]=0.0;	

for (i=0;i<no_elem_dof;i++)                	
	for (j=i;j<no_elem_dof;j++) 
		SMmb[j][i]=SMmt[j][i]=SMmw[j][i]=0.0;	

    C[0][0]=  1.0;
	C[2][0]= -3.0;
	C[3][0]=  2.0;

	C[1][1]=    L;
	C[2][1]= -2*L;
	C[3][1]=    L;

	C[4][2]=  1.0;
	C[8][2]= -35.0;
    C[9][2]=  84.0; 
   C[10][2]= -70.0;
   C[11][2]=  20.0;

	C[5][3]=     L;
	C[8][3]= -20*L;
	C[9][3]=  45*L;	   
   C[10][3]= -36*L;   
   C[11][3]=  10*L;

	C[6][4]=     L2/2.;
	C[8][4]=  -5*L2;
	C[9][4]=  10*L2;
   C[10][4]= -15*L2/2.;
   C[11][4]=   2*L2;
			   
 	C[7][5]=  L3/6;
	C[8][5]=  2*L2/3.0;
	C[9][5]=  L3;
   C[10][5]=  2*L3/3.0;
   C[11][5]=  L3/6;

	C[2][6]=  3.0;
	C[3][6]= -2.0;

	C[2][7]=  -L;
	C[3][7]=   L;

	C[8][8]=  35.0;
	C[9][8]= -84.0;
   C[10][8]=  70.0;
   C[11][8]= -20.0;
   	
   	C[8][9]= -15*L;
	C[9][9]=  39*L;
   C[10][9]= -34*L;
   C[11][9]=  10*L;

   	C[8][10]=  15*L2/6.0;
	C[9][10]= -7*L2;
   C[10][10]=  39.*L2/6;
   C[11][10]= -2*L2;

	C[8][11]=  L3/6;
	C[9][11]=  L3/2;
   C[10][11]= -L3/2;
   C[11][11]=  L3/6;

// Stiffness matrix Kb from paper by F.M. EL-AMIN AND M.A.KASEM Page 165

 SMmb[2][2]=   4./L2;
 SMmb[3][2]=   6./L4;
 SMmb[4][2]=  -2./rL2;
 SMmb[5][2]=  -1./rL2;
 SMmb[6][2]=  -2./(3*rL2);
 SMmb[7][2]=  -1./(2*rL2);
 SMmb[8][2]=  -2./(5*rL2);
 SMmb[9][2]=  -1./(3*rL2);
SMmb[10][2]=  -1./(7*rL2);
SMmb[11][2]=  -1./(4*rL2);

 SMmb[3][3]=   12/L4;
 SMmb[4][3]=  -3./   rL2;
 SMmb[5][3]=  -2./   rL2;
 SMmb[6][3]=  -3./(2*rL2);
 SMmb[7][3]=  -6./(5*rL2);
 SMmb[8][3]=  -1./  (rL2);
 SMmb[9][3]=  -1./(7*rL2);
SMmb[10][3]=  -3./(4*rL2);
SMmb[11][3]=  -2./(3*rL2);

 SMmb[4][4]=    1./(1*r2);
 SMmb[5][4]=    1./(2*r2);
 SMmb[6][4]=    1./(3*r2);
 SMmb[7][4]=    1./(4*r2);
 SMmb[8][4]=    1./(5*r2);
 SMmb[9][4]=    1./(6*r2);
SMmb[10][4]=    1./(7*r2);
SMmb[11][4]=    1./(8*r2);

 SMmb[5][5]=    1./(3*r2);
 SMmb[6][5]=    1./(4*r2);
 SMmb[7][5]=    1./(5*r2);
 SMmb[8][5]=    1./(6*r2);
 SMmb[9][5]=    1./(7*r2);
SMmb[10][5]=    1./(8*r2);
SMmb[11][5]=    1./(9*r2);

 SMmb[6][6]=    1./(5*r2);
 SMmb[7][6]=    1./(6*r2);
 SMmb[8][6]=    1./(7*r2);
 SMmb[9][6]=    1./(8*r2);
SMmb[10][6]=    1./(9*r2);
SMmb[11][6]=   1./(10*r2);

 SMmb[7][7]=    1./(7*r2);
 SMmb[8][7]=    1./(8*r2);
 SMmb[9][7]=    1./(9*r2);
SMmb[10][7]=   1./(10*r2);
SMmb[11][7]=   1./(11*r2);

 SMmb[8][8]=    1./(9*r2);
 SMmb[9][8]=   1./(10*r2);
SMmb[10][8]=   1./(11*r2);
SMmb[11][8]=   1./(12*r2);

 SMmb[9][9]=   1./(11*r2);
SMmb[10][9]=   1./(12*r2);
SMmb[11][9]=   1./(13*r2);

SMmb[10][10]=  1./(13*r2);
SMmb[11][10]=  1./(14*r2);

SMmb[11][11]=  1./(15*r2);

// Stiffness matrix Kt from paper by F.M. EL-AMIN AND M.A.KASEM Page 165

  SMmt[1][1]=SMmt[2][1]=SMmt[3][1]=1./r2L2;
  SMmt[5][1]=SMmt[6][1]=SMmt[7][1]=SMmt[8][1]=SMmt[9][1]=SMmt[10][1]=SMmt[11][1]=1./rL2;
 
  SMmt[2][2]=    4./(3*r2L2);
  SMmt[3][2]=    6./(4*r2L2);

  SMmt[5][2]=      1./(rL2);
  SMmt[6][2]=    4./(3*rL2);
  SMmt[7][2]=    6./(4*rL2);
  SMmt[8][2]=    8./(5*rL2);
  SMmt[9][2]=   10./(6*rL2);
 SMmt[10][2]=   12./(7*rL2);
 SMmt[11][2]=   14./(8*rL2);
 
  SMmt[3][3]=    9./(5*r2L2);

  SMmt[5][3]=      1./(rL2);
  SMmt[6][3]=    6./(4*rL2);
  SMmt[7][3]=    9./(5*rL2);
  SMmt[8][3]=   12./(6*rL2);
  SMmt[9][3]=   15./(7*rL2);
 SMmt[10][3]=   18./(8*rL2);
 SMmt[11][3]=   21./(9*rL2);

 SMmt[5][5]=SMmt[6][5]=SMmt[7][5]=SMmt[8][5]=SMmt[9][5]=SMmt[10][5]=SMmt[11][5]=1./L2;

  SMmt[6][6]=    4./(3*L2);
  SMmt[7][6]=    6./(4*L2);
  SMmt[8][6]=    8./(5*L2);
  SMmt[9][6]=   10./(6*L2);
 SMmt[10][6]=   12./(7*L2);
 SMmt[11][6]=   14./(8*L2);
  
  SMmt[7][7]=    9./(5*L2);
  SMmt[8][7]=   12./(6*L2);
  SMmt[9][7]=   15./(7*L2);
 SMmt[10][7]=   18./(8*L2);
 SMmt[11][7]=   21./(9*L2);

  SMmt[8][8]=   16./(7*L2);
  SMmt[9][8]=   20./(8*L2);
 SMmt[10][8]=   24./(9*L2);
 SMmt[11][8]=  28./(10*L2);

  SMmt[9][9]=   25./(9*L2);
 SMmt[10][9]=  30./(10*L2);
 SMmt[11][9]=  35./(11*L2);

SMmt[10][10]=  36./(11*L2);
SMmt[11][11]=  42./(12*L2);

SMmt[11][11]=  49./(13*L2);

// Stiffness matrix Kw from paper by F.M. EL-AMIN AND M.A.KASEM Page 165

  SMmw[2][2]=    4/r2;
  SMmw[3][2]=    6/r2;

  SMmw[6][2]=    4/r;
  SMmw[7][2]=    6/r;
  SMmw[8][2]=    8/r;
  SMmw[9][2]=    10/r;
 SMmw[10][2]=    12/r;
 SMmw[11][2]=    14/r;
  
  SMmw[3][3]=    12/r2;

  SMmw[6][3]=    6/r;
  SMmw[7][3]=    12/r;
  SMmw[8][3]=    13/r;
  SMmw[9][3]=    24/r;
 SMmw[10][3]=    30/r;
 SMmw[11][3]=    36/r;

  SMmw[6][6]=  4;
  SMmw[7][6]=  6;
  SMmw[8][6]=  8;
  SMmw[9][6]= 10;
 SMmw[10][6]= 12;
 SMmw[11][6]= 14;
  
  SMmw[7][7]= 12;
  SMmw[8][7]= 18;
  SMmw[9][7]= 24;
 SMmw[10][7]= 30;
 SMmw[11][7]= 36;

  SMmw[8][8]=   144./5;
  SMmw[9][8]=   240./6; 
 SMmw[10][8]=   360./7;
 SMmw[11][8]=   504./8;

  SMmw[9][9]=   400./7;
 SMmw[10][9]=   600./8;
 SMmw[11][9]=   840./9;

SMmw[10][10]=   900./9;
SMmw[11][11]=   126./10;
SMmw[11][11]=  1764./11;

// End of miserable typing! : 17/9/96

// Conventions violated ; Write a note on this in Notes.FEM
double J=Ix;
double I=Iy;
double Iw=Iz;

double E,MEU,WD,G;
material->get(E,MEU,WD,G);

double EIL=E*I*L;
double GJL=G*J*L;
double EIwbyL3=E*Iw/L3;

for (i=0;i<no_elem_dof;i++)                	
	for (j=i;j<no_elem_dof;j++) 
// Operation done on diagonal elements also.
		SMm[j][i]=(SMmb[j][i]=(EIL)*SMmb[j][i])+
				  (SMmt[j][i]=(GJL)*SMmt[j][i])+
				  (SMmw[j][i]=(EIwbyL3)*SMmw[j][i]);	

for(i=0;i<no_elem_dof;i++)
	for (j=0;j<no_elem_dof;j++)
		{
		SMs[i][j]=0;	// SMs[] is used for temporary purpose only.
		for(k=0;k<no_elem_dof;k++) 
			SMs[i][j]=SMs[i][j]+C[k][i]*SMm[k][j];
		}

for(i=0;i<12;i++)
	for (j=0;j<12;j++)
		{
		SMm[i][j]=0;
		for(k=0;k<12;k++) 
			SMm[i][j]=SMm[i][j]+SMs[i][k]*C[k][j];
		}

Fill_Upper(SMm);
//fprintf(outfile,"Element stiffness matrix in local coordinates \n");
//Print_M(SMm);
}

void Cbeam_3D_GArch::Compute_SMm(double Axial_Load)
{
// Write this by superposing HArch & VArch

double SMmb[16][16],SMmt[16][16],SMmw[16][16];
int i,j;

for (i=0;i<no_elem_dof;i++)  // Operation done on diagonal elements also.
	for (j=i;j<no_elem_dof;j++) 
		SMm[j][i]=SMmb[j][i]=SMmt[j][i]=SMmw[j][i]=0.0;	

 double L=r*IncAng/180*PI;
 double L2=L*L;
 double L4=L2*L2;
 double r2=r*r;
 double rL=r*L;
 double r2L2=r2*L2;
 double rL2=r*L2;
 double r2L=r2*L;
 double L3=L2*L;

// Stiffness matrix Kt from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402. SMmb[2][2]=  12./L4;

 SMmb[2][2]=  12./L4; 
 SMmb[5][2]=  6./L3; 
 SMmb[3][2]=  6./(5*rL2);
 SMmb[7][2]=  1./(10*rL);
 SMmb[10][2]= -12/L4;
 SMmb[13][2]=  6/L3;
 SMmb[11][2]=  -6./(5*rL2);
 SMmb[15][2]=  1./(10*rL);

 SMmb[5][5]=  4/L2;
 SMmb[3][5]=  11./(10*rL);
 SMmb[7][5]= 2./(15*r);
 SMmb[10][5]= -6/L3;
 SMmb[13][5]= 2/L2;
 SMmb[11][5]= -1./(10*rL);
 SMmb[15][5]= -1./(30*r);

 SMmb[3][3]= 13./(35*r2);
 SMmb[7][3]= 11.*L/(210*r2);
 SMmb[10][3]= -6./(5*rL2);
 SMmb[13][3]= 1./(10*rL);
 SMmb[11][3]= 9./(70*r2); 			
 SMmb[15][3]= -13.*L/(420*r2); 

 SMmb[7][7]= L2/(105*r2);
 SMmb[10][7]= -1./(10*rL);
 SMmb[13][7]= -1./(30*r);
 SMmb[11][7]= 13*L/(420*r2);
 SMmb[15][7]= -L2/(140*r2);

 SMmb[10][10]=  12/L4;   // Value given in journal is negetive 12/L4 !!!!!!!!!!!!!!
 SMmb[13][10]= -6./L3;
 SMmb[11][10]=  6./(5*rL2);		
 SMmb[15][10]=  -1./(10*rL);

 SMmb[13][13]=  4/L2;
 SMmb[11][13]= -11./(10*rL);
 SMmb[15][13]=  2./(15*r);

 SMmb[11][11]=   13./(35*r2);
 SMmb[15][11]=   -11./(210*r2);

 SMmb[15][15]=   L2/(105*r2);

// Stiffness matrix Kt from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402.

  SMmt[2][2]=   6./(5*r2L2); 
  SMmt[5][2]=   1./(10*r2L); 
  SMmt[3][2]=   6./(5*rL2);
  SMmt[7][2]=   1./(10*rL); 
  SMmt[10][2]=   6./(5*r2L2);
  SMmt[13][2]=   1./(10*r2L);
  SMmt[11][2]=  -6./(5*rL2);
  SMmt[15][2]=   1./(10*rL);

  SMmt[5][5]=    2./(15*r2);
  SMmt[3][5]=    1./(10*rL);
  SMmt[7][5]=    2./(15*r);
  SMmt[10][5]=   -1./(10*r2L);
  SMmt[13][5]=   -1./(20*r2);
  SMmt[11][5]=   -1./(10*rL);
  SMmt[15][5]=   -1./(30*r);

  SMmt[3][3]=    6./(5*L2);
  SMmt[7][3]=    1./(10*L);
  SMmt[10][3]=   -6./(5*rL2);
  SMmt[13][3]=    1./(10*rL);
  SMmt[11][3]=   -6./(5*L2);
  SMmt[15][3]=   1./(10*L);

  SMmt[7][7]=  2./15;
  SMmt[10][7]=  1./(10*rL);
  SMmt[13][7]=  -1./(10*rL);
  SMmt[11][7]=  -1./(10*L);
  SMmt[15][7]=  -1./30; 

  SMmt[10][10]=   6./(5*r2L2); 
  SMmt[13][10]=  -1./(10*r2L); 
  SMmt[11][10]=   6./(5*rL2) ;
  SMmt[15][10]=   -1./(10*rL);
 
  SMmt[13][13]=   2./(15*r2);
  SMmt[11][13]=   -1./(10*rL);
  SMmt[15][13]=   2./(15*r);

  SMmt[11][11]=   6./(5*L2);
  SMmt[15][11]=   1./(10*L);

  SMmt[15][15]=  2./15;

// Stiffness matrix Kw from paper by F.M. EL-AMIN AND D.M BROTTON Page 1402.

double r2L3=r2L2*L;
double r2L4=r2L2*L2;
double rL3=r*L3;
double rL4=r*L4;

  SMmw[2][2]=   12./(r2L4);
  SMmw[5][2]=  	6./(r2L3);
  SMmw[3][2]=  	12./(rL4);
  SMmw[7][2]=  	6./(rL3);
  SMmw[10][2]=  	-12/(r2L4);
  SMmw[13][2]=  	6/(r2L3);
  SMmw[11][2]=  	-12/(rL4);
  SMmw[15][2]=  	6/(rL3);

  SMmw[5][5]=  -4/r2L2;
  SMmw[3][5]=  6/(rL3);
  SMmw[7][5]=  4/(rL2);
  SMmw[10][5]=  -6/r2L3;  
  SMmw[13][5]=   2/r2L2; 
  SMmw[11][5]=   -6/(rL3); 
  SMmw[15][5]=   2./(rL2); 
		  
  SMmw[3][3]= 12/L4;
  SMmw[7][3]=  6/L3;
  SMmw[10][3]=  -12/rL4;
  SMmw[13][3]=  6/rL3;
  SMmw[11][3]=  -12/L4;
  SMmw[15][3]=  6/L3;
  
  SMmw[7][7]=  4/L2;
  SMmw[10][7]=  -6/rL3;
  SMmw[13][7]=  2./rL2;
  SMmw[11][7]=  -6/L3;
  SMmw[15][7]=  2/L2;

  SMmw[10][10]=   12/r2L4;
  SMmw[13][10]=   -6/r2L3;
  SMmw[11][10]=   12/rL4;
  SMmw[15][10]=  -6/rL3;
  	   
  SMmw[13][13]=  4/r2L2; 
  SMmw[11][13]=   -6/rL3;
  SMmw[15][13]=  4/rL2; 

  SMmw[11][11]=  12/L4;    
  SMmw[15][11]=  -6/L3;
 
  SMmw[15][15]=  4/L2;


/*
DOF for CBeam_2D_HArch
0 = w
1 = tw
2 = tx
3 = tx'

4 = w
5 = tw
6 = tx
7 = tx'
*/

// Ref: FE Struct. Analysis by T.Y.Yang Page 190.
// Cross checked on 17/9/96. 

 
double b=IncAng*PI/180.0;
double a=Ix/A;// ********VERIFY !!! : 22/9/1996.
L=b*r;
Init_Lower(SMm);
double E,MEU,WD,G;
material->get(E,MEU,WD,G);
//Kuu
double p=E*A/5*(1+a/(r*r));
double q=E*A/r;
double s=E*A/L;

	SMm[0][0]   = 6/L  *p;
	SMm[6][0]   = 1/2.0*p;
    SMm[1][0]   =-1/2.0*q;
    SMm[4][0]   =(L/10+a/L)*q;
	SMm[8][0]   =-6/L   *p;
	SMm[14][0]   = 1/2.0 *p;
	SMm[9][0]   = 1/2.0*q;
	SMm[12][0]   =-(L/10+a/L)*q;


	SMm[6][6]   = 2*L/3 *p;
    SMm[1][6]   =-(L/10+a/L)*q;
    SMm[4][6]   = a/2*q;
	SMm[8][6]   =-1/2.0 *p;
	SMm[14][6]   =-L/6   *p;
	SMm[9][6]   = (L/10+a/L)*q;
	SMm[12][6]   =-(L*L/60+a/2)*q;


  	SMm[1][1]   =(13*b*b/35+12*a/(L*L))*s;
	SMm[4][1]   =(11*L*b*b/210+6*a/L)*s;
	SMm[8][1]   = -1/2.0*q;
	SMm[14][1]   =(L/10+a/L)*q;
	SMm[9][1]   =(9*b*b/70-12*a/(L*L))*s;
	SMm[12][1]   =(-13*L*b*b/420+6*a/L)*s;
	

	SMm[4][4]   =(L*L*b*b/105+4*a)*s;
	SMm[8][4]   = -(L/10+a/L)*q;
	SMm[14][4]   =(L*L/60+a/2)*q;
	SMm[9][4]   =(13*L*b*b/420-6*a/L)*s;		        // corrected on 18/9/96
	SMm[12][4]   =(-L*L*b*b/140+2*a)*s;
	

	SMm[8][8]   = 6/L   *p;
	SMm[14][8]   =-1/2.0   *p;
	SMm[9][8]   = 1/2.0*q;
	SMm[12][8]   = (L/10+a/L)*q;

	SMm[14][14]   = 2*L/3.0 *p;
	SMm[9][14]   =-(L/10+a/L)*q;
	SMm[12][14]   =-a/2*q;

	SMm[9][9]   =(13*b*b/35+12*a/(L*L))*s;
	SMm[12][9]   =-(11*L*b*b/210+6*a/L)*s;

	SMm[12][12]   =(L*L*b*b/105+4*a)*s;
/*
DOF for CBeam_2D_VArch
0 = u		0 x				 
2 = v		1
h0 = w		2
h2 = tx		3
3 = ty		4
h1 = tw		5
1 = u'		6
h3 = tx'	7

4 = u		8
6 = v		9
h4 = w		10
h6 = tx		11
7 = ty		12
h5 = tw		13
5 = u'		14
h7 = tx'	15	   
*/
// Conventions violated ; Write a note on this in Notes.FEM
double J=Ix;
double I=Iy;
double Iw=Iz;

double EIL=E*I*L;
double GJL=G*J*L;
double EIwL=E*Iw*L;

for (i=0;i<no_elem_dof;i++)                	
	for (j=i;j<no_elem_dof;j++) 
// Operation done on diagonal elements also.
		SMm[j][i]+=(SMmb[j][i]=(EIL)*SMmb[j][i])+
				   (SMmt[j][i]=(GJL)*SMmt[j][i])+
				   (SMmw[j][i]=(EIwL)*SMmw[j][i]);	

Fill_Upper(SMm);

double C[16][16];

// ToDo: Insert here the Rotation Transformation matrix

for(i=0;i<no_elem_dof;i++)
	for (j=0;j<no_elem_dof;j++)
		{
		SMs[i][j]=0;	// SMs[] is used for temporary purpose only.
		for(int k=0;k<no_elem_dof;k++) 
			SMs[i][j]=SMs[i][j]+C[k][i]*SMm[k][j];
		}
		

//fprintf(outfile,"Element stiffness matrix in local coordinates \n");
//Print_M(SMm);
}




void Celem::NLA_Operations(int cnt,int count,double *disp)
{
//if (count==0)
{
	Compute_SMm(0);  
	compute_R();	
	Local_2_Global(SMm,SMs);   
}
//else  
//	for (int i=0;i<no_elem_dof;i++)
//    for (int j=0;j<no_elem_dof;j++)
//						SMs[i][j]=1e10;  // Make stiffness infinity 
	Assemble_K();  
}

void Include_NonK(Celem * elem,double E,double A,double L,double alphaz0,double alphay0,double *elem_disp,double Axial_Load)
{
double t0=alphaz0;
double t1=elem_disp[5];
double t2=elem_disp[11];
double u1=elem_disp[0];
double u2=elem_disp[6];
double n1[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],n2[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
   k0[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];

elem->Init_Lower(n1);
elem->Init_Lower(n2);
elem->Init_Lower(k0);

//for n1
double F4= t1+t2-12*t0;
double F51=4*t1-t2-3*t0;
double F52=4*t2-t1-3*t0;

// for n2
double F1=	9*t1*t1		+9*t2*t2	-2*t1*t2	-36*t1*t0	-36*t2*t0	+216*t0*t0;

double F2=	2*t1*t1		+2*t2*t2	-1*t1*t2	-3*t1*t0	-3*t2*t0	+18*t0*t0;
double F7= -2*t1*t1	    -2*t2*t2	+6*t1*t2	-2*t1*t0	-2*t2*t0	-3*t0*t0;

double F31=6*t1*t1		+1*t2*t2	+2*t1*t2	-54*t1*t0	+6*t2*t0	+54*t0*t0;
double F32=6*t2*t2		+1*t1*t1	+2*t1*t2	-54*t2*t0	+6*t1*t0	+54*t0*t0;

double F61=8*t1*t1		+3*t2*t2	-4*t1*t2	-12*t1*t0	-2*t2*t0	+27*t0*t0;
double F62=8*t2*t2		+3*t1*t1	-4*t1*t2	-12*t2*t0	-2*t1*t0	+27*t0*t0;

t0=alphay0;
t1=elem_disp[4];
t2=elem_disp[10];

//for n1
double G4= t1+t2-12*t0;
double G51=4*t1-t2-3*t0;
double G52=4*t2-t1-3*t0;

// for n2
double G1=	9*t1*t1		+9*t2*t2	-2*t1*t2	-36*t1*t0	-36*t2*t0	+216*t0*t0;

double G2=	2*t1*t1		+2*t2*t2	-1*t1*t2	-3*t1*t0	-3*t2*t0	+18*t0*t0;
double G7= -2*t1*t1	    -2*t2*t2	+6*t1*t2	-2*t1*t0	-2*t2*t0	-3*t0*t0;

double G31=6*t1*t1		+1*t2*t2	+2*t1*t2	-54*t1*t0	+6*t2*t0	+54*t0*t0;
double G32=6*t2*t2		+1*t1*t1	+2*t1*t2	-54*t2*t0	+6*t1*t0	+54*t0*t0;

double G61=8*t1*t1		+3*t2*t2	-4*t1*t2	-12*t1*t0	-2*t2*t0	+27*t0*t0;
double G62=8*t2*t2		+3*t1*t1	-4*t1*t2	-12*t2*t0	-2*t1*t0	+27*t0*t0;


n1[1][0]=n2[7][6]=-(n1[7][0]=n1[6][1]=F4/(10*L));
n1[1][1]=n1[2][2]=n1[7][7]=n1[8][8]=-(n1[7][1]=n1[8][2]=-6.0*(u2-u1)/(5*L*L));
n1[4][4]=n1[5][5]=n1[10][10]=n1[11][11]=2*(u2-u1)/15.0;
n1[2][0]=n1[8][6]=-(n1[8][0]=n1[6][2]=-G4/(10*L));
n1[5][0]=-(n1[6][5]=F51/30);
n1[4][0]=-(n1[6][4]=G51/30);
n1[11][0]=-(n1[11][6]=F52/30);
n1[10][0]=-(n1[10][6]=G52/30);
n1[5][1]=n1[11][1]=n1[8][4]=n1[10][8]=(u2-u1)/(10*L);
n1[4][2]=n1[10][2]=n1[11][7]=n1[7][5]=-n1[5][1];
n1[10][4]=n1[11][5]=-(u2-u1)/30;

n2[1][1]=n2[7][7]=-(n2[7][1]=-(F1/100+G2/25)/L);
n2[2][1]=n2[8][7]=-(n2[8][1]=n2[7][2]=F4*G4/(100*L));
n2[2][2]=n2[8][8]=-(n2[8][2]=-(G1/100+F2/25)/L);
n2[5][1]=-(n2[7][5]=-(F31+G2)/300);
n2[4][1]=F4*G51/300;
n2[5][2]=-(n2[8][5]=G4*F51/300);
n2[4][2]=-(n2[8][4]=(G31+F2)/300);
n2[5][5]=(F61/300+G2/225)*L;
n2[5][4]=F51*G51*L/900;

n2[4][4]=(G61/300+F2/225)*L;
n2[7][4]=-F4*G51/300;
n2[11][1]=-(n2[11][7]=-(F32+G2)/300);
n2[10][1]=-(n2[10][7]=-F4*G52/300);
n2[11][2]=-(n2[11][8]=G4*F52/300);
n2[10][2]=-(n2[10][8]=(G32+F2)/300);
n2[11][5]=(F7-G2/3)*L/300;
n2[10][5]=F51*G52*L/900;
n2[11][4]=F52*G51*L/900;
n2[10][4]=(G7-F2/3)*L/900;
n2[11][11]=(F62/300+G2/225)*L;
n2[11][10]=F52*G52*L/900;
n2[10][10]=(G62/300+F2/225)*L;

// ************** K0
t0=alphaz0;
t1=elem_disp[5];
t2=elem_disp[11];

F4= t1+t2-12*t0;
F51=4*t1-t2-3*t0;
F52=4*t2-t1-3*t0;

t0=alphay0;
t1=elem_disp[4];
t2=elem_disp[10];

G4= t1+t2-12*t0;
G51=4*t1-t2-3*t0;
G52=4*t2-t1-3*t0;

k0[1][0]=n2[7][6]=-(k0[7][0]=k0[6][1]=F4/(10*L));
k0[1][1]=k0[2][2]=k0[7][7]=k0[8][8]=-(k0[7][1]=k0[8][2]=-6*(u2-u1)/(5*L*L));
k0[4][4]=k0[5][5]=k0[10][10]=k0[11][11]=2*(u2-u1)/15;
k0[2][0]=k0[8][6]=-(k0[8][0]=k0[6][2]=-G4/(10*L));
k0[5][0]=-(k0[6][5]=F51/30);
k0[4][0]=-(k0[6][4]=G51/30);
k0[11][0]=-(k0[11][6]=F52/30);
k0[10][0]=-(k0[10][6]=G52/30);
k0[5][1]=k0[11][1]=k0[8][4]=k0[10][8]=(u2-u1)/(10*L);
k0[4][2]=k0[10][2]=k0[11][7]=k0[7][5]=-k0[5][1];
k0[10][4]=k0[11][5]=-(u2-u1)/30;

elem->Fill_Upper(n1);
elem->Fill_Upper(n2);
elem->Fill_Upper(k0);
for (int i=0;i<12;i++)
for (int j=0;j<12;j++)
	{
	if (NLM==Direct) 
			SMm[i][j]+=(n1[i][j]+n2[i][j])*E*A;
		else
			{
			SMm[i][j]+=(0.5*n1[i][j]+1.0/3.0*n2[i][j])*E*A;
//			if (!eq(u2-u1,ZERO) && !eq(Axial_Load,ZERO))
//			SMm[i][j]+=Axial_Load*k0[i][j]*L/((u2-u1)*E*A);
			}
	}
}

void Cbeam_3D::NLA_Operations(int cnt,int count,double *disp)
{
int gdir_for[MAX_NO_ELEM_DOF];
double elem_disp[12],Aelem_disp[12],elem_F[12];
double Axial_Load;
double  x1,y1,z1,x2,y2,z2;
//#define prn 1
      //fprintf(outfile," ELEMENT NO = %d \n",elem->no);
      compute_gdir_vector(gdir_for);
     // find initial geometry at start of each step
if ((count==0) && (cnt==0))
	{
	  nodeman->getnode(enode[0]->no,x1,y1,z1);
	  nodeman->getnode(enode[1]->no,x2,y2,z2);			        
      dx0=x2-x1;  							//    horizontal component of length
      dy0=y2-y1;	        				    //    vertical  component of length
	  dz0=z2-z1;
  	}

      nodeman->getnode(enode[0]->no,x1,y1,z1);
	  nodeman->getnode(enode[1]->no,x2,y2,z2);			        
      L=sqrt(sqr(x2-x1)+sqr(y2-y1)+sqr(z2-z1));		//    length						

//#ifdef prn 
fprintf(outfile," Length = %lf coords = %lf,%lf,%lf  %lf,%lf,%lf\n",L,x1,y1,z1,x2,y2,z2); 
//#endif

      double L0=sqrt(sqr(dx0)+sqr(dy0)+sqr(dz0));		//    length						

      double alphaz0= atan2(dy0,dx0);			// 	  angle
      double alphay0= atan2(dx0,dz0);			// 	  angle
      double alphax0= atan2(dz0,dy0);			// 	  angle

      double alphaz=  atan2(y2-y1,x2-x1);			// 	  angle
      double alphay=  atan2(x2-x1,z2-z1);			// 	  angle
      double alphax=  atan2(z2-z1,y2-y1);			// 	  angle

      if ((eq(y2-y1,0) && eq(x2-x1,0)) || (eq(dy0,0) && eq(dx0,0))) { alphaz0=0; alphaz=0; }
	  if ((eq(x2-x1,0) && eq(z2-z1,0)) || (eq(dx0,0) && eq(dz0,0))) { alphay0=0; alphay=0; }
	  if ((eq(z2-z1,0) && eq(y2-y1,0)) || (eq(dz0,0) && eq(dy0,0))) { alphax0=0; alphax=0; }


	  // store global displacement (which is initially [0])  to elem displacement
	  for(int i=0;i<no_elem_dof;i++)
			{
			int dir = gdir_for[i];
			if (dir>-1) 
					{
						elem_disp[i]=disp[dir]; 
						Aelem_disp[i]=disp[dir];
					}
				else 
					{
						elem_disp[i]=0; 
						Aelem_disp[i]=0;
					}
		  	}


  	  double dispe60=Aelem_disp[6]-Aelem_disp[0]; // find x deformation
      double dispe71=Aelem_disp[7]-Aelem_disp[1]; // find y deformation
      double dispe82=Aelem_disp[8]-Aelem_disp[2]; // find z deformation

fprintf(outfile," ELEM %d \n",no);
fprintf(outfile," ax = %lf %lf ,ay = %lf %lf ,az = %lf %lf \n",alphax0,alphax,alphay0,alphay,alphaz0,alphaz);

      // find displacements corresponding to member straining modes

      elem_disp[0]=0.0;			// axial 1
      elem_disp[1]=0.0;			// lateral y 1
	  elem_disp[2]=0.0;			// lateral z 1
	  elem_disp[3]=elem_disp[3]-(alphax-alphax0);			// torsion 1
	  elem_disp[4]=elem_disp[4]-(alphay-alphay0);	// theta y 1
      elem_disp[5]=elem_disp[5]-(alphaz-alphaz0);	// theta z 1

      elem_disp[6]=0;
      //
      //1/(L0+L)*((2*dx0+dispe60)*dispe60+(2*dy0+dispe71)*dispe71
      //				+(2*dz0+dispe82)*dispe82);	// axial  2

      elem_disp[7]=0.0;		// lateral y 2
      elem_disp[8]=0.0;		// lateral z 2

	  elem_disp[9] =elem_disp[9]-(alphax-alphax0);	// torsion 2  ??
      elem_disp[10]=elem_disp[10]-(alphay-alphay0);	// theta y 2			
      elem_disp[11]=elem_disp[11]-(alphaz-alphaz0);  // theta z 2

for (i=0;i<12;i++)  if (eq(elem_disp[i],0)) elem_disp[i]=0;  

#ifdef prn 
fprintf(outfile," Element Disp in global dirns \n");
for (i=0;i<6;i++)  fprintf(outfile," %30.10lf ",elem_disp[i]);
fprintf(outfile,"\n");
for (i=6;i<12;i++)  fprintf(outfile," %30.10lf ",elem_disp[i]);
fprintf(outfile,"\n");
#endif




  //elem_disp[9]-=elem_disp[3];
  //elem_disp[3]=0.0;

  compute_R();
  Change_2_Local(elem_disp);
  elem_disp[6]=1/(L0+L)*((2*dx0+dispe60)*dispe60+(2*dy0+dispe71)*dispe71
      				+(2*dz0+dispe82)*dispe82);	
  
	double E,MEU,WD,G;
  	material->get(E,MEU,WD,G);

	double yy=(alphay-alphay0);
	double zz=(alphaz-alphaz0);

	Axial_Load=(E*A/L)/(L0+L)*((2*dx0+dispe60)*dispe60+(2*dy0+dispe71)*dispe71
    				+(2*dz0+dispe82)*dispe82); // NLIncr or Lincr


/*
	if (((cnt==0) && (NLM==NLincr)) || (NLM==Direct)) Axial_Load=-AML[0];
	 	else //if (count==0) 
		{
	Axial_Load=A*E*(dispe60/L+0.5*sqr(dispe71/L)+0.5*sqr(dispe82/L))
	+1.0/30*(2*sqr(Aelem_disp[4]-yy)-(Aelem_disp[4]-yy)*(Aelem_disp[10]-yy)
	+2*sqr(Aelem_disp[10]-yy)
	+2*sqr(Aelem_disp[5]-zz)-(Aelem_disp[5]-zz)*(Aelem_disp[11]-zz)
	+2*sqr(Aelem_disp[11]-zz));

//	fprintf(outfile,"\nAxial Load for elem %d = %lf \n",no,Axial_Load);
	     }
*/

	 // find axial force in element


	double temp=L;
	double temp11,temp12,temp13;
    double temp21,temp22,temp23;
	
	temp11=enode[0]->x[0];
    temp12=enode[0]->x[1];
    temp13=enode[0]->x[2];
    
    temp21=enode[1]->x[0];
    temp22=enode[1]->x[1];
    temp23=enode[1]->x[2];

    enode[0]->x[0]=0;
    enode[0]->x[1]=0;
    enode[0]->x[2]=0;
    
    enode[1]->x[0]=dx0;
    enode[1]->x[1]=dy0;
    enode[1]->x[2]=dz0;
    L=L0;

	Compute_SMm(Axial_Load);  // Geometric stiffness is included in SMm


  	L=temp;
    enode[0]->x[0]=temp11;
    enode[0]->x[1]=temp12;					 
    enode[0]->x[2]=temp13;
    
    enode[1]->x[0]=temp21;
    enode[1]->x[1]=temp22;
    enode[1]->x[2]=temp23;

//    compute_R();	// find R for new geometry

#ifdef prn 
	fprintf(outfile,"\nTransformation Matrix\n");
	for (i=0;i<3;i++)
	{
	for (int j=0;j<3;j++)
		fprintf(outfile," %10.4lf ",R[i][j]);
	fprintf(outfile,"\n");
	}
#endif

// No more > Include_NonK(this,E,A,L0/*L*/,zz,yy,Aelem_disp,Axial_Load);

#ifdef prn 
fprintf(outfile," Element Stiffness Matrix in local \n");
Print_M(SMm);
#endif

    Local_2_Global(SMm,SMs);   // Transform SMm(local) to SMs(global) 

#ifdef prn 
fprintf(outfile," Element Stiffness Matrix in Global \n");
Print_M(SMs);
#endif


#ifdef prn 
fprintf(outfile," Element Disp in local dirns \n");
for (i=0;i<no_elem_dof;i++)  fprintf(outfile," %30.25lf ",elem_disp[i]);
fprintf(outfile,"\n");
#endif

/*
if (((no==1) || (no==20)) && (cnt==1) && (count==0)) 
	{
	fprintf(outfile," SMm %d\n ",no);
	Print_M(SMm);	
	}
*/

for (i=0;i<no_elem_dof;i++)  elem_F[i]=0.0;

if (NLM==NLincr)
	{
	// find internal member forces in global directions
    for (i=0;i<no_elem_dof;i++)
	      for (int j=0;j<no_elem_dof;j++)
		 	  	 {
		 	  	 elem_F[i]=elem_F[i]+SMm[i][j]*elem_disp[j]; 
//   fprintf(outfile," elem_F[%d][%d] %15.3lf =  %15.3lf %15.3lf, \n",i,j,elem_F[i],SMm[i][j],elem_disp[j]);
				 }

#ifdef prn 
fprintf(outfile," Element Forces in local dirns \n");
for (i=0;i<no_elem_dof;i++)  fprintf(outfile," %lf ",elem_F[i]);
fprintf(outfile,"\n");
#endif

	Change_2_Global(elem_F);

//#ifdef prn 
fprintf(outfile," Element Forces in global dirns \n");
for (i=0;i<6;i++)  fprintf(outfile," %30.10lf ",elem_F[i]);
fprintf(outfile,"\n");
for (i=6;i<12;i++)  fprintf(outfile," %30.10lf ",elem_F[i]);
fprintf(outfile,"\n");
//#endif

    for (i=0;i<no_elem_dof;i++)
 	  		if (gdir_for[i]>-1) AE[gdir_for[i]]+=elem_F[i];

	}

    Assemble_K();  
}


