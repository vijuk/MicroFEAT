// NON LINEAR ANALYSIS USING DIRECT METHOD
#include <stdafx.h>
#include "space.h"
     
void Cstructure::NLA_Iteration(int Nsteps)
{
int gdir_for[MAX_NO_ELEM_DOF];
double elem_disp[12],disp[100],store_AE[100];

int count =0, stop =0; 
for (int i=0;i<no_dof;i++) 
			{ disp[i]=0; d[i]=0; store_AE[i]=0;}

double  x1,y1,z1,x2,y2,z2;
do {
fprintf(outfile,"************* COUNT = %d \n",count);

	  Celem *elem=elemman->first;
	  while (elem)
	  {

if (count==0)
	{	  
	  // find initial geometry	
	  nodeman->getnode(elem->enode[0]->no,x1,y1,z1);
	  nodeman->getnode(elem->enode[1]->no,x2,y2,z2);			        
      
      elem->dx0=x2-x1;  							//    horizontal component of length
      elem->dy0=y2-y1;	        				    //    vertical  component of length
	  elem->dz0=z2-z1;
   
	}

    double L0=sqrt(sqr(elem->dx0)+sqr(elem->dy0)+sqr(elem->dz0));		//    length						

      double alphaz0= atan2(elem->dy0,elem->dx0);			// 	  angle
      double alphay0= atan2(elem->dx0,elem->dz0);			// 	  angle
      double alphax0= atan2(elem->dz0,elem->dy0);			// 	  angle
	  
	  // store global displacement (which is initially [0])  to elem displacement
	  elem->compute_gdir_vector(gdir_for);
	  for(int i=0;i<elem->no_elem_dof;i++)
			{
			int dir = gdir_for[i];
			if (dir>-1) 
					elem_disp[i]=disp[dir]; 
				else 
					elem_disp[i]=0; 
		  	}

//
	 // Update Coordinates of element

/*
	  for (i=0;i<3;i++)
		{
	  	elem->enode[0]->x[i]+=elem_disp[i];
    	elem->enode[1]->x[i]+=elem_disp[6+i];
		}

      nodeman->getnode(elem->enode[0]->no,x1,y1,z1);
	  nodeman->getnode(elem->enode[1]->no,x2,y2,z2);			        
*/
      
      elem->L=sqrt(sqr(x2-x1)+sqr(y2-y1)+sqr(z2-z1));		//    length						

      double alphaz= atan2(y2-y1,x2-x1);			// 	  angle
      double alphay= atan2(x2-x1,z2-z1);			// 	  angle
      double alphax= atan2(z2-z1,y2-y1);			// 	  angle

//

  	  double dispe60=elem_disp[6]-elem_disp[0]; // find axial deformation
      double dispe71=elem_disp[7]-elem_disp[1]; // find lateral deformation
      double dispe82=elem_disp[8]-elem_disp[2]; // find lateral deformation

      // find displacements corresponding to member straining modes
      elem_disp[0]=0.0;			// axial 1
      elem_disp[1]=0.0;			// lateral y 1
	  elem_disp[2]=0.0;			// lateral z 1
	  elem_disp[3]=0.0;			// torsion 1
	  elem_disp[4]=0.0; //elem_disp[4]-(alphay-alphay0);	// theta y 1
      elem_disp[5]=elem_disp[5]-(alphaz-alphaz0);	// theta z 1

      elem_disp[6]=1/(L0+elem->L)*((2*elem->dx0+dispe60)*dispe60+(2*elem->dy0+dispe71)*dispe71);
      			//	+(2*elem->dz0+dispe82)*dispe82);	// axial  2

      elem_disp[7]=0.0;		//  lateral y 2
      elem_disp[8]=0.0;		// lateral z 2

	  elem_disp[9] =0; //elem_disp[9]-(alphay-alphay0);	// torsion 2  ??
      elem_disp[10]=0; //elem_disp[10]-(alphay-alphay0);	// theta y 2			
      elem_disp[11]=elem_disp[11]-(alphaz-alphaz0);  // theta z 2

//for (i=0;i<elem->no_elem_dof;i++) fprintf(outfile,"elem_disp[%d] =  %20.8lf \n",i,elem_disp[i]);

	  // remove
	  double E,MEU,WD,G;

	  elem->material->get(E,MEU,WD,G);
	  double A=elem->A;
	  double Ix=elem->Ix;
	  double Iy=elem->Iy;
	  double Iz=elem->Iz;

	if (count!=0) elem->AML[0]=((E*A)/elem->L)*elem_disp[6];  // find axial force in element
//     elem->AML[0]=0;

	 elem->Compute_SMm(0);  // Geometric stiffness is included in SMm
   
      // find internal member forces in local axis
//if (count!=0)
	 for (i=0;i<elem->no_elem_dof;i++)
	  {
	  elem->AML[i]=0.0;
 	  //if (gdir_for[i]>-1) 
      for (int j=0;j<elem->no_elem_dof;j++)
 	  if (gdir_for[j]>-1) 
     		elem->AML[i]=elem->AML[i]+SMm[i][j]*d[gdir_for[j]]; // elem_disp[j];
	  }

	for (i=0;i<elem->no_elem_dof;i++) 
		if (gdir_for[i]>-1) fprintf(outfile,"d[%d] =  %20.8lf \n",i,d[gdir_for[i]]);

	for (i=0;i<elem->no_elem_dof;i++) 
		fprintf(outfile,"AML[%d] =  %20.8lf \n",i,elem->AML[i]);

    elem->compute_R();	// find R for new geometry
    elem->Local_2_Global(SMm,SMs);   // Transform SMm(local) to SMs(global) 
	elem->Change_2_Global(elem->AML);
	elem->compute_P();
    elem->Assemble_K();  
	elem=elem->next;
} 

	// The incremental displacement d[] due to load imbalance is computed      
	// using stiffness matrix in the current displacement pattern
   
	Band_Solve();
     
    // disp[] is updated 
	for (int i=0;i<no_dof;i++) 
		{
			disp[i]=d[i];    // D = D + delta(D)
			store_AE[i]+=d[i];
			if (fabs(d[i]) > 1e-10) stop=0; else stop=1;
			for (int j=0;j<no_dof;j++) 
			K[i][j]=0.0;
		}

count++;
} while (!stop && count<100);  
fprintf(outfile,"Final deflections \n");
for (i=0;i<no_dof;i++) fprintf(outfile," %d = %12.8lf \n",i,store_AE[i]);
}

// matout(disp,no_dof,1);
// end main()

/*void matout(double y[100][100],int nrow,int ncolumn)
      {
      for (int i=0;i<nr;i++)
	  	{
      	for (int j=0;j<nc;j++)
      		fprintf(out," %lf ",y[i][j]);
		fprintf(out,"\n");
		}
      }
*/

       



/*
       void beamstifforce(force,Iy,A,L,elemo)
       implicit double precision *8(a-dx,o-z)
       dimension elemo[6][6]
       do 5 i=1,6
       do 5 j=1,6
5      elemo[i][j]=0.0
       elemo(2,2)=36
       elemo(2,3)=3*L
       elemo(2,5)=-36
       elemo(2,6)=3*L
       elemo(3,3)=4*(L**2)
       elemo(3,5)=-3*L
       elemo(3,6)=-1*(L**2)
       elemo(5,5)=36
       elemo(5,6)=-3*L
       elemo[6][6]=4*(L**2) 
       do 6 i=1,6
       do 6 j=1,6
6      elemo[j][i]=elemo[i][j]
       do 7 i=1,6
       do 7 j=1,6
7      elemo[i][j]=(force/(30*L))*elemo[i][j]
       return
       end
*/
