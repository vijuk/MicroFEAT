/*
Last modified at 11.15 pm on 29/6/96

*******************************************************************************
		Description  :          Space Frame Program
								for Static & Dynamic Analysis
								(Uses Straight and Curved beam elements).

		Coding       :          Viju K S
		Start date   :          04-06-1996
		End date     :          ??-??-1996
*******************************************************************************/

#include "stdafx.h"
//#include "MicroFEAT.h"
//#include "mainfrm.h"
//#include "Microdoc.h"
//#include "Microvw.h"
 
#ifdef _DEBUG
#undef THIS_FILE
static char BASED_CODE THIS_FILE[] = __FILE__;
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include "space.h"

// Global variables 
 
 int Sol=-1;
 BOOL dyna=0;
 BOOL NLA=0;
 int NLM = Direct;
 int Nsteps = 10;
 
CString InpFileName,OutFileName;

 FILE *inpfile,*outfile;
 double   *AE=NULL;
 double   *d=NULL;
 double   *eig_val=NULL; // double[no_dof]
 double   **K=NULL;	     	 // Stiffness matrix : double[no_dof][sbw]
 double   **M=NULL;		     // Mass matrix      : double[no_dof][sbw]
 double   **eig_vect=NULL;      // Eigen Vectors    : double[no_dof][no_dof]

 double *SK=NULL;

 double   SMm[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
 double   SMs[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];                
 double   MMm[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
 double   MMs[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];	  
 double   SMmP[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
 double   SMsP[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];  
 
 double   R[3][3];    				  

 int no_nodes=0;
 int no_rest_nodes=0;
 int no_elem_groups=0;
 int no_elems[MAX_NO_ELEM_GROUPS];
 int no_loaded_nodes[10]={0,0,0,0,0  ,0,0,0,0,0};
 int no_loaded_elems=0;
 int sbw=0;
 int no_dof=0;

 
 Cnodeman *nodeman=NULL;			// node manager
 Celemman *elemman=NULL;			// element manager
 Cmaterialman *materialman=NULL; 	// material manager
 Cloadman *loadman=NULL;			 	// load manager

// *********START THREE DIMENSIONAL ANALYSIS GLOBAL VARIABLES HERE

// IMP : B matrix is not used to find stiffness matrix 
//		 It is used only to find out stresses only.


double c1,c2,c3;
int no_materials,no_dist_loads;

int  *ColHeight,*DiagDOF;  
double *stress;		
double *Q;		// Nodal Load vector = Qface+Qbody (constant for all load cases)
double *P;		// Nodal Load vector = Q*ELM[loadcase_no]+Qnodal[loadcase_no]

double face_centre_rst[7][3]={	0,	0,	0,  		// centre
							  	1,	0,	0,			// +r
							   -1,	0,	0,			// -r
							  	0,	1,	0,			// +s
							  	0, -1,	0,			// -s
							  	0,	0,	1,			// +t
							  	0,	0, -1,			// -t
							 };

double Gauss_Weight[5][5],Gauss_Point[5][5];  //change to 4,4
double rst[3];


double JM[3][3]; 				// Jacobian Matrix
double Inv_JM[3][3];			// Inverse Jacobian Matrix
double N[MAX_NO_ELEM_DOF];		// Shape functions
double B[6][MAX_NO_ELEM_DOF];	// Strain displacement matrix


double C[6][6];					// Constitutive matrix
double CB[6][MAX_NO_ELEM_DOF];	// Stress displacement matrix

double DN_rst[3][MAX_NO_ELEM_NODES];			// Derivative of shape function wrt. r,s,t
double DN_xyz[MAX_NO_ELEM_NODES][3];    		// Derivative of shape function wrt. x,y,z

double DET_J;	   				// Determinant Jacobian

//double DET_J_C;					// Determinant Jacobian at Centre (incompatiable)
//double JM_C[3][3];			// Jacobian at centre for incompatiable modes
//double Inv_JM_C[999][999];		// Inverse of JM_C

int per[]={1,2,0};				// To perform permutation of indexes in loops
int stress_Oface[7];			// Face numbers for stress output for each elem
int no_sky_dof;
double Qbody[MAX_NO_ELEM_DOF];
double Qface[MAX_NO_ELEM_DOF];

int no_load_cases;// no of load cases to be analysed

int elem_LM[10];

// *********END THREE DIMENSIONAL ANALYSIS GLOBAL VARIABLES HERE

 void disp_error(int err_no);

// ******* CStructure definitions

void Cstructure::Read_Control()
{

/*
CFileDialog dial(1);
dial.DoModal();
InpFileName=dial.GetPathName();
dial.DoModal();
OutFileName=dial.GetPathName();

CString InpFileName="solid.inp";
CString OutFileName="solid.OUT";
*/

if ((inpfile=fopen(InpFileName,"r"))==NULL)
	{
	printf("Input file not found\n");
	exit(1);
	}

	outfile=fopen(OutFileName,"w");

    // Read Modulas of elasticity & poisson's ratio for structure 
/*
	if (!D3)
	{
	fscanf(inpfile,"%lf %lf\n",&E_global,&G_global);
	fprintf(outfile,"E= %12.3lf   G = %12.3lf\n",E_global,G_global);
	}
*/
	fscanf(inpfile,"%d %d %d %d %d",&Sol,&dyna,&NLA,&NLM,&Nsteps);
	fprintf(outfile," Sol = %d , Dyna = %d , NLA = %d Nsteps = %d\n",Sol,dyna,NLA,Nsteps);
    // Read no. of (nodes,element groups & restrained nodes )
	fscanf(inpfile,"%d %d %d",&no_nodes,&no_elem_groups,&no_rest_nodes);
	
	fprintf(outfile,
		"No. of nodes            = %d\n"
		"No. of element groups   = %d\n"
		"No. of restrained nodes = %d\n",
			 no_nodes,no_elem_groups,no_rest_nodes);

// Read no. of material groups and no, of distributed loads,no. of load cases
//D3
		fscanf(inpfile,"%d %d %d",&no_materials,&no_dist_loads,&no_load_cases);
		fprintf(outfile,
		"No. of materials            = %d\n"
		"No. of distributed loads    = %d\n"
		"No. of load cases			 = %d\n",
			 no_materials,no_dist_loads,no_load_cases);		

	/* Read no. of (loaded nodes & loaded elements) */

	fscanf(inpfile,"%d %d",&no_loaded_nodes[0],&no_loaded_elems);
	fprintf(outfile,"No. of loaded nodes    = %d\n"
		"No. of loaded elements = %d\n",no_loaded_nodes[0],no_loaded_elems);

    // Read no. of elements for each element group

	int elem_group_no[MAX_NO_ELEM_GROUPS];
	for (int i=0;i<no_elem_groups;i++)
		{
		fscanf(inpfile,"%d",&elem_group_no[i]);
		fscanf(inpfile,"%d",&no_elems[i]);
		fprintf(outfile,"Group %d :No. of elements = %d \n"
				,elem_group_no[i],no_elems[i]);
		}

	fprintf(outfile,"Note: Group 1 = beam_3D , Group 2 = beam_3D_curved \n");

	/* Removed on 25/9/96.
	switch (elem_group_no[0]){
	ASSIGN_NO_NODE_DOF_MACRO();
	default:{
		fprintf(outfile,"Element group not available in library.\n");
		exit(1);
		}
		}
	*/

//D3
for(i=1;i<=no_materials;i++)  materialman->add(i);

	// Create linked list for nodes
	for(i=1;i<=no_nodes;i++) nodeman->add(i);

//(D3)
for(i=1;i<=no_dist_loads;i++) loadman->add(i);

	// Create linked list for elements
	if (no_elem_groups<MAX_NO_ELEM_GROUPS)
	{
	for(i=0;i<no_elem_groups;i++) 
				elemman->add_elem_group(elem_group_no[i],no_elems[i]);
	}
	else
	{
	fprintf(outfile,"Error : Maximum number of element groups exceeded.\n");
	fprintf(outfile,"Maximum number allowed = %d\n",MAX_NO_ELEM_GROUPS);
	}
}


// Allocate and Initialize global matrices
BOOL Cstructure::Alloc_Init_Memory(int no_dof,int sbw)
{
int i;


if (!(AE      =new double[no_dof])) return 0;
if (!(d       =new double[no_dof])) return 0;

/*
delete AE;
delete d;
AE=0;
d=0;

return 1;
*/

if ((Sol<=BANDED) || (Sol==ALL)) 
	{
	if (!(K       =new double* [no_dof])) return 0;
	for (i=0;i<no_dof;i++) if (!(K[i]=new double[sbw])) return 0;

	if (dyna)
		{
			if (!(M       =new double* [no_dof])) return 0;
			for (i=0;i<no_dof;i++) 
				if (!(M[i]=new double[sbw])) return 0;
		}
	}


//DYNAMIC MATRICES
if (dyna) {
	if (!(eig_val =new double[no_dof])) return 0;
    if (!(eig_vect=new double* [no_dof])) return 0;
				for (i=0;i<no_dof;i++) 
   					if (!(eig_vect[i]=new double[no_dof])) return 0;
	}

// Initialize global arrays to zero

for (i=0;i<no_dof;i++) {
						AE[i]=0.0;
						d[i]=0.0;
						int j;
						if ((Sol<=BANDED) || (Sol==ALL)) 
							for (j=0;j<sbw;j++) K[i][j]=0.0;
					
						//DYNAMIC MATRICES
						if (dyna) {
							eig_val[i]=0.0;
							for (j=0;j<no_dof;j++) eig_vect[i][j]=0.0;
						    for (j=0;j<sbw;j++) M[i][j]=0.0;
							}     
					
 						}

if ((Sol==SKYLINE) || (Sol==ALL)) 
	{
	if (!(SK       =new double[no_sky_dof])) return 0;
	for (i=0;i<no_sky_dof;i++) SK[i]=0.0;
	}


if (!(Q       =new double[no_dof])) return 0;
if (!(P       =new double[no_dof])) return 0;

for (i=0;i<no_dof;i++) P[i]=Q[i]=0.0;

return 1; // memory allocation success
}


Cstructure::~Cstructure()
{

fclose(inpfile);
fclose(outfile);

if (AE) delete AE;
if (d) delete d;

if ((Sol<=BANDED) || (Sol==ALL)) if (K) 
	{
	for (int i=0;i<no_dof;i++) 	
						delete K[i];
	delete K;
	}

// FREE ALL DYNAMICS MATRICES.
if (dyna) {
			if (M) 
				{
				for (int i=0;i<no_dof;i++) delete M[i];
				delete M;
				}

			if (eig_vect)
				{
				for (int i=0;i<no_dof;i++) delete eig_vect[i];
				delete eig_vect;
				}
			
			if (eig_val) delete eig_val;
    	  }

if ((Sol==SKYLINE) || (Sol==ALL)) 
{
if (SK) delete SK;
if (DiagDOF) delete DiagDOF;
if (ColHeight) delete ColHeight;
}


if (Q) delete Q;
if (P) delete P;

if (nodeman) delete nodeman;
if (elemman) delete elemman;
if (loadman) delete loadman;
if (materialman) delete materialman;
}


void disp_error(int err_no)
{
fprintf(outfile,"LOGIC/INPUT ERROR(%2d):*******%%%%%%%%+++++++++@@@@@@@@\n",err_no);
}


void Cstructure::Display_Displacements()
{ 
/*
fprintf(outfile,"Deflection at each dof.\n");
for(int i=0;i<no_dof;i++)
	fprintf(outfile,"Defl. @ %d =  %16.8lf\n",i,d[i]);
*/

/*
int gdir_for[MAX_NO_ELEM_DOF];
Celem * elem=elemman->first;
while (elem)
	{
	elem->compute_gdir_vector(gdir_for);
	for (i=0;i<elem->no_elem_dof;i++)
		{
		int k=gdir_for[i];
		if (k>-1)
		 fprintf(outfile," ELEM %d : DOF %d = %12.6lf \n",elem->no,i,d[k]);
		}
	elem=elem->next;
	}
*/

int onode1,onode2;fscanf(inpfile,"%d %d",&onode1,&onode2);
Cnode * node=nodeman->first;
while (node)
	{
	for (int i=0;i<node->no_node_dof;i++)
	   {
		int k=node->ID[i];
		if (k>-1)
		if ((onode1==node->no)||(onode2==node->no)) 
		fprintf(outfile," NODE %d : DOF %d = %12.6lf \n",node->no,i,d[k]);
	   }
	node=node->next;
	}
}


/* SOLUTION OF EQUATIONS using Guass elemination */
BOOL Cstructure::Band_Solve()
{
int 	n,p,i,j,k;
double 	temp;

for (i=0;i<no_dof;i++) d[i]=AE[i];

for (n=0;n<no_dof;n++)
{
if (eq(K[n][0],0)) 
		{
			fprintf(outfile,"Stiffness matrix is not positive definite.\n");
			return 0; // Failed 
		}
i=n+1;
for (p=1;p<sbw;p++,i++)
		if (i<no_dof)
		{
		temp=K[n][p]/K[n][0]; /* divide 2..dof by diaganal */
		j=0;
		for (k=p;k<sbw;k++,j++)
			 K[i][j] -= temp*K[n][k];
		K[n][p]=temp;
		d[i]-=temp*d[n];
		}
d[n]/=K[n][0];
}

//Back substitution
for (n=(p=no_dof-1)-1;n>=0;p=n--)      
		for (k=1;k<sbw;k++,p++)		d[n]-=K[n][k]*d[p];

//for (i=0;i<no_dof;i++) fprintf(outfile," %5d = %10.8lf \n",d[i]);
return 1;  // Success 
}

/* SOLUTION OF EQUATIONS using Guass elemination */
BOOL Cstructure::Simple_Solve()
{
for (int i=0;i<no_dof;i++) d[i]=AE[i];

for (i=0;i<no_dof;i++)
{
if (eq(K[i][i],0.0)) 
		{
			fprintf(outfile,"Stiffness matrix is not positive definite.\n");
			return 0; // Failed 
		}

for (int k=i+1;k<no_dof;k++)
		{
		double t=K[k][i]/K[i][i];
		for (int j=0;j<no_dof;j++) K[k][j]-=K[i][j]*t;
		d[k]-=d[i]*t;
		}
}


// Back substitution phase 
d[no_dof-1]/=K[no_dof-1][no_dof-1];
for (int k=no_dof-2;k>=0;k--)
		{
		double t=0;
		for (int j=no_dof-1;j>=k+1;j--) 
			t=t+K[k][j]*d[j];
		d[k]=(d[k]-t)/K[k][k];
		}

//for ( i=0;i<no_dof;i++) fprintf(outfile," %5d = %10.8lf \n",d[i]);
return 1;
}


int eq(double x,double y)
{
if (fabs(x-y)<2*ZERO) return 1; else return 0;
}


void Cstructure::Initialize_Gauss()	 // thrd20
{

Gauss_Weight[0][0]=2;
Gauss_Weight[1][0]=0;
Gauss_Weight[2][0]=0;
Gauss_Weight[3][0]=0;

Gauss_Weight[0][1]=.3351800554016;
Gauss_Weight[1][1]=Gauss_Weight[0][1];
Gauss_Weight[2][1]=0;
Gauss_Weight[3][1]=0;

Gauss_Weight[0][2]=0.5555555555556;
Gauss_Weight[1][2]=Gauss_Weight[0][2];
Gauss_Weight[2][2]=0.8888888888889;
Gauss_Weight[3][2]=0;

Gauss_Weight[0][3]=0.3478546451375;
Gauss_Weight[1][3]=0.6521451548625;
Gauss_Weight[2][3]=Gauss_Weight[1][3];
Gauss_Weight[3][3]=Gauss_Weight[0][3];

Gauss_Point[0][0]=0;
Gauss_Point[1][0]=0;
Gauss_Point[2][0]=0;
Gauss_Point[3][0]=0;

Gauss_Point[0][1]=-0.7587869106393;
Gauss_Point[1][1]=-Gauss_Point[0][1];
Gauss_Point[2][1]=0;
Gauss_Point[3][1]=0;

Gauss_Point[0][2]=-0.7745966692415;
Gauss_Point[1][2]=-Gauss_Point[0][2];
Gauss_Point[2][2]=0;
Gauss_Point[3][2]=0;

Gauss_Point[0][3]=-0.86113631115941;
Gauss_Point[1][3]=-0.3399810435849;
Gauss_Point[2][3]=-Gauss_Point[1][3];
Gauss_Point[3][3]=-Gauss_Point[0][3];
//

Gauss_Point[0][0]=  0.000000000000000  ;Gauss_Weight[0][0]= 2.000000000000000  ;

Gauss_Point[0][1]= -0.577350269189626  ;Gauss_Weight[0][1]= 1.000000000000000  ;
Gauss_Point[1][1]= +0.577350269189626  ;Gauss_Weight[1][1]= 1.000000000000000  ;

Gauss_Point[0][2]= -0.774596669241483  ;Gauss_Weight[0][2]= 0.555555555555556  ;
Gauss_Point[1][2]=  0.000000000000000  ;Gauss_Weight[1][2]= 0.888888888888889  ; 					
Gauss_Point[2][2]= +0.774596669241483  ;Gauss_Weight[2][2]= 0.555555555555556  ;

Gauss_Point[0][3]= -0.861136311594053  ;Gauss_Weight[0][3]= 0.347854845137454  ;
Gauss_Point[1][3]= -0.339981043584856  ;Gauss_Weight[1][3]= 0.652145154862546  ;
Gauss_Point[2][3]= +0.339981043584856  ;Gauss_Weight[2][3]= 0.652145154862546  ;
Gauss_Point[3][3]= +0.861136311594053  ;Gauss_Weight[3][3]= 0.347854845137454  ;

}

void Cstructure::Triangularize_SK()
{
int k,ic,klt;
for(int n=0;n<no_dof;n++)
		{
		int kn=DiagDOF[n];	// addr of n th diag
		int kl=kn+1;		
		int ku=DiagDOF[n+1]-1;  // addr of n+1 th diag
		int kh=ku-kl;			// colheight for n row - 1 

		if (kh>0)			// if colheight > 1
			{
			 k=n-kh-1;			 // row no.of lowest non-0 elem in col
			 ic=0;				 // 0 > -1
			 klt=ku;			 // addr of lowest " non-0 elem
			for (int j=0;j<kh;j++)
				{
				ic++;		// ic increasing 0 to kh
				klt--;		// klt decreasing ku to kl
			    k++;
				int ki=DiagDOF[k];		   // addr of k th diag
				int ND=DiagDOF[k+1]-ki-1;  // col height of k th col 
				
				if (ND>0)  // 0 > -1 
					{
				    int kk=min(ic,ND);
					double c=0;
					for (int s=1;s<=kk;s++)
						c+=SK[ki+s]*SK[klt+s];
					SK[klt]-=c;
					}
				}
			}

		if (kh>=0)
				{
				k=n;		// k  > k+1
				double b=0;
				if (k>0)
					for (int kk=kl;kk<=ku;kk++)
						{
						k--;
						int ki=DiagDOF[k];
						double c=SK[kk]/SK[ki];
						b+=c*SK[kk];
						SK[kk]=c;	  // L
						}
				SK[kn]-=b;			  // D
				}
		if (SK[kn]<0) ; //RE_MSG(4);
		}

/*
fprintf(outfile,"Printing SK[] /n");
for (int i=0;i<no_sky_dof;i++) 
		fprintf(outfile,"%d = %12.3lf \n",i,SK[i]);
*/
return;
}



BOOL Cstructure::Sky_Solve()
{
// Reduce Right Hand Side Load Vector

for(int n=0;n<no_dof;n++)
	{
	int kl=DiagDOF[n]+1;
	int ku=DiagDOF[n+1]-1;
	if (ku-kl>=0)
		{
		int k=n;
		double c=0;
		for (int kk=kl;kk<=ku;kk++)
			{
			k--;
			c+=SK[kk]*P[k];
			}
			P[n]-=c;
		}
	}

// Back Substitution
for (n=0;n<no_dof;n++)
	{
	int k=DiagDOF[n];
	P[n]/=SK[k];
	}

if (no_dof==1) return TRUE;  // only one dof

	n=no_dof;
	for (int l=1;l<no_dof;l++)
		{
		int kl=DiagDOF[n]+1;
		int ku=DiagDOF[n+1]-1;
		if (ku-kl>=0)
			{
			int k=n;
			for (int kk=kl;kk<=ku;kk++)
				{
				k--;
				P[k]-=SK[kk]*P[n];
				}
			}
		n--;
		}

/*
fprintf(outfile,"Deflections \n");
for (int i=0;i<no_dof;i++) 
		fprintf(outfile," DOF %d = %12.6lf \n",i,P[i]);
*/

/*int gdir_for[MAX_NO_ELEM_DOF];
Celem * elem=elemman->first;
while (elem)
	{
	elem->compute_gdir_vector(gdir_for);
	for (i=0;i<elem->no_elem_dof;i++)
		{
		int k=gdir_for[i];
		if (k>-1)
		fprintf(outfile," ELEM %d : DOF %d = %12.6lf \n",elem->no,i,P[k]);
		}
	elem=elem->next;
	}
*/

double TIP[13][3];  // remove 
Cnode * node=nodeman->first;
while (node)
	{
	for (int i=0;i<node->no_node_dof;i++)
		{
		int k=node->ID[i];
		if (k>-1)
			{

			if (node->no==1) TIP[0][i]=P[k];
			if (node->no==9) TIP[1][i]=P[k];
			if (node->no==13) TIP[2][i]=P[k];
			if (node->no==21) TIP[3][i]=P[k];
			if (node->no==22) TIP[4][i]=P[k];
			if (node->no==29) TIP[5][i]=P[k];
			if (node->no==33) TIP[6][i]=P[k];
			if (node->no==34) TIP[7][i]=P[k];
			if (node->no==41) TIP[8][i]=P[k];
			if (node->no==42) TIP[9][i]=P[k];
			if (node->no==49) TIP[10][i]=P[k];
			if (node->no==53) TIP[11][i]=P[k];
			if (node->no==54) TIP[12][i]=P[k];

			fprintf(outfile," NODE NO %3d : DOF %3d = %12.6lf \n",node->no,i,P[k]);
			}
		}
	node=node->next;
	}

// This block is for a special case : remove after use.
for (int j=0;j<13;j++)
for (int i=0;i<3;i++)
    fprintf(outfile," ***** TIPNODE %d(%1d) = %12.8lf \n",j,i,TIP[j][i]);
fprintf(outfile,"\n");
for (int i=0;i<3;i++)
	fprintf(outfile," ###### CENTRAL TIP NODE(%1d)  = %12.8lf \n",i,TIP[5][i]);

return TRUE;
}

int Cstructure::Compute_DiagDOF_sky()
{
for (int i=0;i<no_dof+1;i++) DiagDOF[i]=0;

DiagDOF[0]=0;
DiagDOF[1]=1;
int nsbw=0;                        
//fprintf(outfile,"DiagDOF[]\n");

for (i=1;i<no_dof;i++)
		{
		if (ColHeight[i]>sbw) sbw=ColHeight[i];
		DiagDOF[i+1]=DiagDOF[i]+ColHeight[i]+1;
//		fprintf(outfile,"%3d ",DiagDOF[i]);
		}
nsbw++;
int no_sky_dof=DiagDOF[no_dof];
fprintf(outfile,"\n No.sky.dof =  %d  sbw = %d \n",no_sky_dof,sbw);
return no_sky_dof;
}


Cstructure::Cstructure()
{

nodeman=new Cnodeman;
elemman=new Celemman;

// 3D added on 21/11/96 for D3
loadman=new Cloadman;
materialman=new Cmaterialman;

this->Initialize_Gauss();
this->Read_Control();               // Read control data  
materialman->Read_Properties();     // Read properties for all materials 
nodeman->Read_Coords();             // Read coordinates for all nodes
elemman->Read_Properties();         // Read properties for all elements
no_dof = nodeman->Create_IDarray(); // Create ID array and returns no_dof

if ((Sol==SKYLINE) || (Sol==ALL)) 
	{
		DiagDOF       =new int[no_dof+1];
		ColHeight     =new int[no_dof];
	}

if ((Sol==SKYLINE)||(Sol==ALL)) 
		{
			elemman->Compute_ColHeight();		 // Find height of column above diagonal
			no_sky_dof=Compute_DiagDOF_sky();
		}

if ((Sol<=BANDED) ||(Sol==ALL)) // IMP : this also computes Length of elem. Change!
		sbw  =elemman->compute_SBW();       // Compute semi-bandwidth

if (Sol==SIMPLE) sbw=no_dof;  // temporary to test dynamic part

if (!(this->Alloc_Init_Memory(no_dof,sbw)))// Allocate memory to global matrices
	{
		fprintf(outfile,"Failure in allocating memory!\n");
		fprintf(outfile,"No. of degrees of freedom = %d )!\n",no_dof);
		fprintf(outfile,"Shutdown all other programs and try again or\n");
		exit(1);
	}


//!D3) 
// removed nodeman->Read_Loads();              // Read forces acting on node in S-axis
elemman->Read_Loads();			    // Read forces acting on elem in S-axis

//D3
loadman->Read_Loads();				 // Read values for all loads 

// Assembles Global Stiffness matrix

/*
if (!NLA) 
	{
	if ((Sol<=BANDED) || (Sol==ALL))  elemman->SetUp_BandedSystem();			  
	if ((Sol==SKYLINE) || (Sol==ALL)) 
			{
			elemman->SetUp_SkySystem();
			this->Triangularize_SK();
			}
	}
Enable this block & remove from loop below : Qbody assembled here not used.
*/


loadman->Read_Multipliers(); // Read element load multipliers

for (int loadcase_no=0;loadcase_no<no_load_cases;loadcase_no++)
	{
	loadman->Multiply_Q_toget_P(); 	//  P{} = Q{} * elem_LM[load_case_no]
	loadman->Read_Qnodal_2_P(loadcase_no);  	//  P{} = P{} + Qnodal{}

	// Assembles Global Stiffness matrix which is a function of load
if (!NLA)  
		{
		if ((Sol<=BANDED) || (Sol==ALL))  elemman->SetUp_BandedSystem();			  
		if ((Sol==SKYLINE) || (Sol==ALL)) 
				{
				elemman->SetUp_SkySystem();
				this->Triangularize_SK();
				}

		if ((Sol==SKYLINE) || (Sol==ALL))  if (!this->Sky_Solve()) RE_MSG(4);
		if ((Sol==BANDED)  || (Sol==ALL))  if (!this->Band_Solve()) RE_MSG(4);	 // Solve equations in banded form
		if ((Sol==SIMPLE)  || (Sol==ALL))  if (!this->Simple_Solve()) RE_MSG(4); // Solve equations in banded form
		}

	if (NLA) 	NLA_Iteration(Nsteps);     // Iterate until convergence

if (Sol!=SKYLINE) this->Display_Displacements();      // Display results 
	}

if (dyna) this->compute_EigenVectors();
}

// Caution: memeory not allocated for variable stress[]



