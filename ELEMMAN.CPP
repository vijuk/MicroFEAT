#include "stdafx.h"
#include "space.h"	 
#include "constant.h"

DECLARE_NDOFMAP_MACRO();
 

// ******** Celemman definition
Celemman::Celemman()
{
first=0;
last=0;
}

Celemman::~Celemman()
{
if (first) delete first;
}

Celem* Celemman::getelem(int n)
{
Celem* T;
T=first;
while (T)
		{
			if (T->no==n) break;
			T=T->next;
		}
if (T==0) {
			fprintf(outfile,"Input Error. Element manager failed.");
			fclose(outfile);
			exit(1);
		  }
return T;
}



void Celemman::Read_Properties()
{
Celem *elem;
elem=first;
while (elem) 
		{
		elem->Read_Properties();
		elem=elem->next;
		}						   	
}

void Celemman::Compute_ColHeight()
{
Celem *elem;
elem=first;
while (elem) 
		{
		elem->Update_Colum_Heights();
		elem=elem->next;
		}						   	
}



int Celemman::compute_SBW()
{
double x1,y1,z1,x2,y2,z2;
int sbw=0,elem_sbw;
Celem *Telem;

Telem=first;
while (Telem)   {

/* Compute length of each member */
nodeman->getnode(Telem->enode[0]->no,x1,y1,z1);
nodeman->getnode(Telem->enode[1]->no,x2,y2,z2);
Telem->L = sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2));

// Compute semi band width = max(max_dof-min_dof+1 for all elements)
elem_sbw=Telem->get_elem_sbw();
if (elem_sbw>sbw) sbw= elem_sbw;
Telem=Telem->next;
}

fprintf(outfile,"Semibandwidth = %d\n",sbw);
return sbw;
}



void Celemman::add(int n,Celem* new_elem)
{
if (new_elem==0) 
		fprintf(outfile,"Not enough memory to add element!\n");
Celem *T;
if (!first) {
		last=first=T=new_elem;
		T->no=n;
		}
else
		{
		T=first;			        // Added on 24/7/96
		while (T->next) T=T->next;
		last=T->next=new_elem;
		T->next->no=n;
		}
}

void Celemman::add_elem_group(int elem_group_no,int no_elems_for_group)
{
for (int i=1;i<=no_elems_for_group;i++) 
					NEW_FROM_ELEM_LIB_MACRO();
}

void Celemman::Read_Loads()          // Read loads on element wrt. L-axis
{
Celem* Telem;
int node_no;
fprintf(outfile,"Member forces\n");
for (int j=1;j<=no_loaded_elems;j++)
		{
		fscanf(inpfile,"%d ",&node_no);
		fprintf(outfile,"Member : %d \n",node_no);
		Telem=getelem(node_no);
		Telem->Read_Loads();
		}
}

void Celemman::Post_Analysis()
{
Celem *elem=first;
while (!elem) 
	{
		elem->Compute_Stresses();
		elem->Compute_Principal_Stresses();	 
		elem=elem->next;
	}
}

void Celemman::SetUp_SkySystem()
{
Celem *elem=first;
while (elem) {
elem->Initialize_C_B_SMs_Qbody(); 				// INITIALIZE C,SMm,Qbody,B
elem->Compute_C(c1,c2,c3);
elem->Compute_B_Matrix();
elem->Compute_SMs_Qbody(c1,c2,c3); // Csolid08 = 2 ; Csolid20 = ?
elem->Compute_CB();
elem->Assemble_SK();
elem->Compute_Qface();
elem->Assemble_Q();
if ((Sol<=BANDED) || (Sol==ALL)) elem->Assemble_K();
elem=elem->next;
}
}


void Celemman::SetUp_BandedSystem() 
{  
Celem* elem;                              
elem=first;
while (elem)  
// Assembling of global stiffness and mass matrices 
	{
	elem->Compute_SMm(0);	     // Compute element stiffness matrix
	if (dyna) elem->compute_MMm(); // Compute element mass matrix
	elem->compute_R();	// Compute Local->Global transformation matrix
	elem->Local_2_Global(SMm,SMs);//Compute global element stiffness matrix 
	if (dyna) elem->Local_2_Global(MMm,MMs);// global element mass matrix
	elem->compute_P(); // Compute nodal load vector
	elem->Assemble_K();

	/*
	fprintf(outfile,"Element %d: > K \n",elem->no);
	for(int i=0;i<no_dof;i++)
		{
		for(int j=0;j<sbw;j++)
			fprintf(outfile," %8.4lf ",K[i][j]);
		fprintf(outfile,"\n");
		}
	if (dyna)
	{
	fprintf(outfile,"Element %d: > M \n",elem->no);
	for(int i=0;i<no_dof;i++)
		{
		for(int j=0;j<sbw;j++)
			fprintf(outfile," %8.4lf ",M[i][j]);
		fprintf(outfile,"\n");
		}
	}
	*/
	elem=elem->next;
	}
}
