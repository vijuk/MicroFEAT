#ifndef _INC_SPACE_H
#define _INC_SPACE_H


#include "constant.h"
#include "node.h"
#include "elem.h"
#include "nodeman.h"
#include "elemman.h"
#include <math.h>

#define  sqr(x)     ((x)*(x))
int eq(double x,double y);
void disp_error(int i);

extern int Sol;  // this variable decides the solution method of equations
extern BOOL dyna;
extern BOOL NLA;
extern int NLM; 
extern int Nsteps;

extern CString InpFileName,OutFileName;

extern FILE *inpfile,*outfile;
extern double   *AE,*d,*eig_val; // double[no_dof]
extern double   **K;	     	 // Stiffness matrix : double[no_dof][sbw]
extern double   **M;		     // Mass matrix      : double[no_dof][sbw]
extern double   **eig_vect;      // Eigen Vectors    : double[no_dof][no_dof]

extern double   SMm[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
extern double   SMs[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];                
extern double   MMm[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
extern double   MMs[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];	  
extern double   SMmP[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];
extern double   SMsP[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF];  
extern double   R[3][3];    				


extern int no_nodes,no_elems[MAX_NO_ELEM_GROUPS],no_rest_nodes;
extern int no_loaded_nodes[10],no_loaded_elems;
extern int sbw,no_dof;


extern Cnodeman *nodeman;
extern Celemman *elemman;
extern Cmaterialman *materialman;
extern Cloadman *loadman;


class Cstructure 		   // Class containing general structure fns.
	{
	public:		
	Cstructure();    
    ~Cstructure();
	void Read_Control();            // Read control data  
	BOOL Alloc_Init_Memory(int no_dof,int sbw); // Alloc memory & initialize
//	void SetUp_BandedSystem();				// Assembles banded stiffness matrix
	
	// added for 3D analysis
	int Compute_DiagDOF_sky();  //Fill DiagDOF array and return no_sky_dof;

	BOOL Simple_Solve();			// Solve equations in banded form
	BOOL Band_Solve();			// Solve equations in banded form
	BOOL Sky_Solve();           // Solve Equations in sky line form

	void Triangularize_SK();
	void Initialize_Gauss();	 // thrd20
	
	
	void Display_Displacements();   // Display results 
	BOOL compute_EigenVectors();
    void NLA_Iteration(int Nsteps);
	};






// *********START THREE DIMENSIONAL ANALYSIS GLOBAL VARIABLES HERE

// IMP : B matrix is not used to find stiffness matrix 
//		 Itis used only to find out stresses only.
extern double face_centre_rst[7][3];
/*={	0,	0,	0,  		// centre
							  	1,	0,	0,			// +r
							   -1,	0,	0,			// -r
							  	0,	1,	0,			// +s
							  	0, -1,	0,			// -s
							  	0,	0,	1,			// +t
							  	0,	0, -1,			// -t
							 }
*/

extern double c1,c2,c3;
extern double rst[3];

extern int no_materials,no_dist_loads;
extern double Gauss_Weight[5][5],Gauss_Point[5][5];

extern double JM[3][3]; 				// Jacobian Matrix
extern double Inv_JM[3][3];			    // Inverse Jacobian Matrix
extern double N[MAX_NO_ELEM_DOF];		// Shape functions
extern double B[6][MAX_NO_ELEM_DOF];	// Strain displacement matrix

extern double *Q;		// Nodal Load vector = Qface+Qbody (constant for all load cases)
extern double *P;		// Nodal Load vector = Q*ELM[loadcase_no]+Qnodal[loadcase_no]

extern double *SK,*stress;
extern int *ColHeight,*DiagDOF;

extern double C[6][6];					// Constitutive matrix
extern double CB[6][MAX_NO_ELEM_DOF];	// Stress displacement matrix

extern double DN_rst[3][MAX_NO_ELEM_NODES];			// Derivative of shape function wrt. r,s,t
extern double DN_xyz[MAX_NO_ELEM_NODES][3];    		// Derivative of shape function wrt. x,y,z

extern double DET_J;	   				// Determinant Jacobian

//double DET_J_C;					// Determinant Jacobian at Centre (incompatiable)
//double JM_C[3][3];			// Jacobian at centre for incompatiable modes
//double Inv_JM_C[999][999];		// Inverse of JM_C

extern int per[];				// To perform permutation of indexes in loops
extern int stress_Oface[7];			// Face numbers for stress output for each elem
extern int no_sky_dof;
extern double Qbody[MAX_NO_ELEM_DOF];
extern double Qface[MAX_NO_ELEM_DOF];
extern int no_load_cases,elem_LM[10];
// *********END THREE DIMENSIONAL ANALYSIS GLOBAL VARIABLES HERE


#endif








