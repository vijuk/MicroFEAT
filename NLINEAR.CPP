// NON LINEAR ANALYSIS INCREMENTAL AND DIRECT METHODS
#include <stdafx.h>
#include "space.h"

void Cstructure::NLA_Iteration(int Nsteps)
{

double *disp     = new double[no_dof];
double *pd    = new double[no_dof];
double *store_AE = new double[no_dof];
double *stepd    = new double[no_dof];
double *tempd    = new double[no_dof];

//#define prn 
// NLM = 0 means Direct Incremental method
// NLM = 1 means Linear Incremental method
// NLM = 2 means Non-Linear Incremental method

int steps=Nsteps;
if (NLM==Direct) steps=1;

for (int i=0;i<no_dof;i++) 	
	{ 
	disp[i]=0; 
	d[i]=0; 
	pd[i]=0;
	store_AE[i]=AE[i];
	stepd[i]=0;
	}


int cnt=0;
ITER: 
int count =0, stop =0; 
if (NLM==NLincr)
for (i=0;i<no_dof;i++) 	{
						d[i]=0;
						pd[i]=0;
						}

do {
//fprintf(outfile,"*CNT = %d  COUNT = %d \n",cnt,count);
for (i=0;i<no_dof;i++) 	AE[i]=0;
Celem *elem=elemman->first;
while (elem)
	{
	elem->NLA_Operations(cnt,count,disp);
	elem=elem->next;
	}
	 
for (i=0;i<no_dof;i++)
	AE[i]=store_AE[i]/Nsteps*(cnt+1)-AE[i];


/*
if (count!=0) // Zeroing all forces on non Cbeam_3D elems for fixator
{
int gdir_for[MAX_NO_ELEM_DOF];
Celem *elem=elemman->first;
while (elem)
	{
	if (elem->elem_ID!=ID_Cbeam_3D) 
     	{
     	elem->compute_gdir_vector(gdir_for);
		for (int i=0;i<elem->no_elem_dof;i++)
			if (gdir_for[i]>-1) AE[gdir_for[i]]=0;
		}
	elem=elem->next;
	}
}
*/

Band_Solve();

// Note: This statement is not legal. 
//for (i=0;i<no_dof;i++) if (fabs(ceil(d[i]*1e10))<1e-7) d[i]=0;

//for (i=0;i<no_dof;i++) d[i]=int(d[i]*1e8)/1e8;

//#endif

//if (count!=0) // Zeroing all displacements on non Cbeam_3D elems for fixator
/*{
int gdir_for[MAX_NO_ELEM_DOF];      
Celem *elem=elemman->first;
while (elem)
	{
	if (elem->elem_ID!=ID_Cbeam_3D) 
     	{
     	elem->compute_gdir_vector(gdir_for);
		for (int i=0;i<elem->no_elem_dof;i++)
			if (gdir_for[i]>-1) 
					{
						tempd[gdir_for[i]]=d[gdir_for[i]];
						d[gdir_for[i]]=0;
					}

		}
	elem=elem->next;
	}
}
*/

if (NLM==NLincr)	// Update Coordinates of element 
	{
	Cnode *node=nodeman->first;
	while (node) {
	    for (i=0;i<3;i++)
			if (node->ID) if (node->ID[i]>-1) node->x[i]+=d[node->ID[i]];
		node=node->next;
		}
	}

// for (i=0;i<6;i++) fprintf(outfile,"*** %lf %lf\n",d[i],d[i+no_dof-6]);

/*
{
int gdir_for[MAX_NO_ELEM_DOF];
Celem *elem=elemman->first;
while (elem)
	{
	if (elem->elem_ID!=ID_Cbeam_3D) 
     	{
     	elem->compute_gdir_vector(gdir_for);
		for (int i=0;i<elem->no_elem_dof;i++)
			if (gdir_for[i]>-1) 
					{
						d[gdir_for[i]]=tempd[gdir_for[i]];
					}

		}
	elem=elem->next;
	}
}
*/

// disp[] is updated 
stop=TRUE;
for (i=0;i<no_dof;i++) 
		{
			if (fabs(pd[i]-d[i])>1e-3) stop=FALSE;
			pd[i]=d[i];    // D = D + delta(D)
			disp[i]+=d[i]; 
			if (count==0) stepd[i]=0; 
			stepd[i]+=d[i]; 
			for (int j=0;j<sbw;j++) K[i][j]=0.0;
		}					   




/*
if (count!=0) // Zeroing all displacements on non Cbeam_3D elems for fixator
{
int gdir_for[MAX_NO_ELEM_DOF];
Celem *elem=elemman->first;
while (elem)

	{
	if (elem->elem_ID!=ID_Cbeam_3D) 
     	{
     	elem->compute_gdir_vector(gdir_for);
		for (int i=0;i<elem->no_elem_dof;i++)
			if (gdir_for[i]>-1) disp[gdir_for[i]]=0;
		}
	elem=elem->next;
	}
}
*/

count++;
//if (count>2) exit(0);
 } while (!stop && count<25);  

fprintf(outfile,"AE D step disp\n");
    for (i=0;i<no_dof;i++)
		fprintf(outfile,"%d  %30.20lf %30.20lf  %30.20lf %30.20lf \n",i,AE[i],d[i],stepd[i],disp[i]);

//fprintf(outfile,"cnt =  %d count %d \n",cnt,count);

if (NLM==NLincr) 
  {
	//for (i=0;i<no_dof;i++) disp[i]+=d[i]; // Updated displacements
	cnt++;
	if (steps>cnt) goto ITER;  // $$
  }

//fprintf(outfile,"Final deflections ***** \n");
for (i=0;i<no_dof;i++) 
{
d[i]=disp[i];
//fprintf(outfile," %d = %12.8lf \n",i,disp[i]);
}



//Export output to GRAPHER

/*
FILE *grffile=fopen("GRAF.DAT","w");
i=0;
Cnode *node=nodeman->first;
while (node) {
	if (!node->ID) continue; // Unused nodes have ID=NULL.
	double u=0;
	if (node->ID[0]>-1) u=disp[node->ID[0]];
	double v=0;
	if (node->ID[1]>-1) v=disp[node->ID[1]];

	fprintf(grffile,"%5d %lf %lf \n",i,u,v);
	node=node->next;
	i++;
	}
fclose(grffile);
*/

FILE *deffile=fopen("d:\\fem\\feat\\ifm.out","a");
Celem *elem=elemman->first;
while (elem)
	{
	fprintf(deffile,"%3d %3d %10.3lf %10.3lf %10.3lf  %3d %10.3lf %10.3lf %10.3lf 2 \n",
	elem->no,elem->enode[0]->no,elem->enode[0]->x[0],elem->enode[0]->x[1],elem->enode[0]->x[2],
	         elem->enode[1]->no,elem->enode[1]->x[0],elem->enode[1]->x[1],elem->enode[1]->x[2]);
	elem=elem->next;
	}
fclose(deffile);

delete(disp);
delete(pd);
delete(store_AE);
delete(stepd);
delete(tempd);
}
