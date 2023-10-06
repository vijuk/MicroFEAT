// Modify this file when adding new elements. 
// Follow the 7 NE_Step for adding a new element to LIBRARY

#ifndef _INC_CONSTANT_H
#define _INC_CONSTANT_H

// CONSTANTS

// Solution methods
#define  SIMPLE 	0
#define  BANDED     1
#define  SKYLINE    2
#define  FRONT      3
#define  ALL	  	4

// Nonlinearity methods
#define Direct      0
#define Lincr       1
#define NLincr      2

//Time domain analysis methods
#define NEWMARK 1
#define WILSONT 2
#define CENTRAL 3
#define HOBOLT  4

#define  MAX_NO_ELEM_GROUPS   50
#define  MAX_NO_ELEM_NODES    20           
#define  MAX_NO_ELEM_DOF      60
#define  MAX_NO_NODE_DOF	   6 // this affect input of ID array and nodal loads
								 // change also reqd in preprocessor's MAX_NO_NODE_DOF 

// ******** NE_Step(1): Declare element group numbers here

#define  ID_Cbeam_3D           1
#define  ID_Cbeam_3D_curved	   2
#define  ID_Cbeam_2D_VArch	   3
#define  ID_Cbeam_2D_HArch	   4
#define  ID_Cbeam_2D_HArch_12  5
#define  ID_Cbeam_3D_GArch 	   6
#define  ID_Ctriangle_Constant 	   20  // added on 08/10/96
#define  ID_Cquad_iso8             21  // added on 08/11/96

#define  ID_Csolid20               31  // added on 08/11/96
#define  ID_Ctruss_3D              41  
#define  ID_Cgap_3D                42  

// END NE_Step(1)

// ******** NE_Step(2): Delcare No. of element dof for each element class here
#define  NEDOF_Cbeam_3D               12
#define  NEDOF_Cbeam_3D_curved        12
#define  NEDOF_Cbeam_2D_VArch          8
#define  NEDOF_Cbeam_2D_HArch          8
#define  NEDOF_Cbeam_2D_HArch_12      12
#define  NEDOF_Cbeam_3D_GArch         16
#define  NEDOF_Ctriangle_Constant      6
#define  NEDOF_Cquad_iso8         	   8

#define  NEDOF_Csolid20         	  60

#define  NEDOF_Ctruss_3D               6  
#define  NEDOF_Cgap_3D                 6    

// END NE_Step(2)

// ******** NE_Step(3): Delcare No. of node dof for each element class here
#define  NNDOF_Cbeam_3D               6
#define  NNDOF_Cbeam_3D_curved        6
#define  NNDOF_Cbeam_2D_VArch         4
#define  NNDOF_Cbeam_2D_HArch         4
#define  NNDOF_Cbeam_2D_HArch_12      6
#define  NNDOF_Cbeam_3D_GArch         8
#define  NNDOF_Ctriangle_Constant     2
#define  NNDOF_Cquad_iso8	          2

#define  NNDOF_Csolid20 	          3

#define  NNDOF_Ctruss_3D              3  
#define  NNDOF_Cgap_3D                3  
// END NE_Step(3)

// ********* NE_Step(4): Declare No. of nodes for each element class here
#define  NENODES_Cbeam_3D           2
#define  NENODES_Cbeam_3D_curved    2
#define  NENODES_Cbeam_2D_VArch     2
#define  NENODES_Cbeam_2D_HArch     2
#define  NENODES_Cbeam_2D_HArch_12  2
#define  NENODES_Cbeam_3D_GArch     2
#define  NENODES_Ctriangle_Constant      3
#define  NENODES_Cquad_iso8		         8

#define  NENODES_Csolid20		         20
#define  NENODES_Ctruss_3D               2  
#define  NENODES_Cgap_3D               2  
// END NE_Step(4)


/*
0 = u					displacement in X-dir
1 = v                   displacement in Y-dir
2 = w					displacement in Z-dir
3 = tx					Rotation in X-dir
4 = ty					Rotation in Y-dir
5 = tw					Rotation in Z-dir
6 = d(u)/dx				1st derivative of u
7 = d(tx)/dx			1st derivative of tx
8 = d(ty)/dx			1st derivative of ty
9 = d(tw)/dx            1st derivative of tw

10 = d2(u)/dx2          2nd derivative of u
11 = d2(tx)/dx2			2nd derivative of tx
12 = d2(ty)/dx2			2nd derivative of ty
13 = d2(tw)/dx2         2nd derivative of tw
*/


// ********* Declare dof ID's here ( add only if extra dof is present )
#define  T_X            0
#define  T_Y            1
#define  T_Z            2

#define  R_X            3
#define  R_Y            4
#define  R_Z            5

#define  T_X_1X			6
#define  R_X_1X   		7
#define  R_Y_1X   		8
#define  R_W_1X			9

#define  T_X_2X			10
#define  R_X_2X   		11
#define  R_Y_2X   		12
#define  R_W_2X			13
//


// ********* NE_Step(5): Declare Nodal DOFMAP for element here 
// If all nodes have similar dof a only one declaration is needed

#define NDOFMAP_Cbeam_3D                 1<<T_X | \
										 1<<T_Y | \
										 1<<T_Z | \
										 1<<R_X | \
										 1<<R_Y | \
										 1<<R_Z


#define NDOFMAP_Cbeam_3D_curved		     1<<T_X \
										|1<<T_Y \
										|1<<T_Z \
										|1<<R_X \
										|1<<R_Y \
										|1<<R_Z 
									
#define NDOFMAP_Cbeam_2D_VArch 		     1<<T_X \
										|1<<T_Y \
										|1<<R_Z \
										|1<<T_X_1X 

#define NDOFMAP_Cbeam_2D_HArch			 1<<T_X \
										|1<<T_Y \
										|1<<T_Z \
										|1<<R_X \
										|1<<R_Y \
										|1<<R_Z  

#define NDOFMAP_Cbeam_2D_HArch_12		 1<<T_X \
										|1<<T_Y \
										|1<<T_Z \
										|1<<R_X \
										|1<<R_Y \
										|1<<R_Z 

#define NDOFMAP_Cbeam_3D_GArch			 1<<T_X \
										|1<<T_Y \
										|1<<T_Z \
										|1<<R_X \
										|1<<R_Y \
							    		|1<<R_Z 

#define NDOFMAP_Ctriangle_Constant		 1<<T_X \
										|1<<T_Y 

#define NDOFMAP_Cquad_iso8				 1<<T_X \
										|1<<T_Y 



#define NDOFMAP_Csolid20			 	 1<<T_X \
										|1<<T_Y \
										|1<<T_Z 

#define NDOFMAP_Ctruss_3D                1<<T_X | \
										 1<<T_Y | \
										 1<<T_Z  

#define NDOFMAP_Cgap_3D                  1<<T_X | \
										 1<<T_Y | \
										 1<<T_Z  
// END NE_Step(5)


//unsigned  int NDOFA_Cbeam_3D_curved[]   ={NDOFMAP_Cbeam_3D ,NDOFMAP_Cbeam_3D }; 
//#define X_NDOFMAP_MACRO(CELEM) extern unsigned  int NDOFA_##CELEM ={NDOFMAP_##CELEM  ,NDOFMAP_##CELEM  } 

#define _NDOFMAP_MACRO(CELEM)\
	unsigned  int EDMA_##CELEM[]={NDOFMAP_##CELEM  ,NDOFMAP_##CELEM  };\

// ******* NE_Step(6): Declare No of nodal DOF Arrays for all nodes of element type.

#define DECLARE_NDOFMAP_MACRO() \
 _NDOFMAP_MACRO(Cbeam_3D) \
 _NDOFMAP_MACRO(Cbeam_3D_curved) \
 _NDOFMAP_MACRO(Cbeam_3D_GArch) \
 _NDOFMAP_MACRO(Cbeam_2D_VArch) \
 _NDOFMAP_MACRO(Cbeam_2D_HArch) \
 _NDOFMAP_MACRO(Cbeam_2D_HArch_12)\
 _NDOFMAP_MACRO(Ctruss_3D)\
 _NDOFMAP_MACRO(Cgap_3D)\
 unsigned  int EDMA_Ctriangle_Constant[]={NDOFMAP_Ctriangle_Constant,\
 										  NDOFMAP_Ctriangle_Constant, \
 										  NDOFMAP_Ctriangle_Constant}; \
 unsigned  int EDMA_Cquad_iso8[]        ={NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8,\
										  NDOFMAP_Cquad_iso8};\
 unsigned  int EDMA_Csolid20[]          ={NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  NDOFMAP_Csolid20,\
										  };


//elems with variable dof can  be defined below
//unsigned  int NDOFA_Cbeam_3D[] ={NDOFMAP_Cnew_elem_1  ,NDOFMAP_Cnew_elem_2, ... }

//CELEM::enode_dofmap=NDOFA_##CELEM;
// END NE_Step(6)
 

#define _CREATE_NEW_ELEM_MACRO(CELEM) \
  case ID_##CELEM: \
  add(i,new CELEM(ID_##CELEM,NEDOF_##CELEM,NENODES_##CELEM,EDMA_##CELEM)); \
  break;

// Used in element manger.

// ******* NE_Step(7): Register element to the element manager
#define NEW_FROM_ELEM_LIB_MACRO() \
switch (elem_group_no) \
	{ \
_CREATE_NEW_ELEM_MACRO(Cbeam_3D) \
_CREATE_NEW_ELEM_MACRO(Cbeam_3D_curved) \
_CREATE_NEW_ELEM_MACRO(Cbeam_2D_VArch) \
_CREATE_NEW_ELEM_MACRO(Cbeam_2D_HArch_12) \
_CREATE_NEW_ELEM_MACRO(Cbeam_2D_HArch) \
_CREATE_NEW_ELEM_MACRO(Cbeam_3D_GArch) \
_CREATE_NEW_ELEM_MACRO(Ctriangle_Constant) \
_CREATE_NEW_ELEM_MACRO(Csolid20) \
_CREATE_NEW_ELEM_MACRO(Ctruss_3D) \
_CREATE_NEW_ELEM_MACRO(Cgap_3D) \
default:fprintf(outfile,"Error: Incorrect element group no.(%d)\n",elem_group_no); \
	} 

// _CREATE_NEW_ELEM_MACRO(Cquad_iso8) \  Add this when Cquad_iso_8 is used. on 22/11/96

#define  PI         3.141592653589793240
#define  ZERO       1e-8
extern   FILE * outfile;
#define fpout(s) fprintf(outfile,s)

// Define DEBUG VERSION error messages.
#define DE_MSG(_Message_No){\
switch (_Message_No)\
{\
case 1  : fpout("Debug: Element node has zero degrees of freedom.\n");\
case 3  : fpout("Debug: Element dofmap already allocated.\n");\
default : fpout("Debug: Unknown Error occured!.\n");\
}\
fclose(outfile);\
/*exit(_Message_No);*/}

// Define RELEASE VERSION error messages.			   
#define RE_MSG(_Message_No)\
{switch (_Message_No){\
case   1: fpout("Run Time: Error 1:reading ID matrix.\n");\
case   2: fpout("Run Time: Error 2:searching nodes.\n");\
default : fpout("Run Time: Unknown Error occured!\n");\
}\
fclose(outfile);\
/*exit(_Message_No);*/}

#endif


