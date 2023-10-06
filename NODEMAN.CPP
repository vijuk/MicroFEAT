#include "stdafx.h"
#include "space.h"

//************** Cnodeman definitions 

Cnodeman::Cnodeman()
{
first=0;
last=0;
}

Cnodeman::~Cnodeman()
{
if (first) delete first;
}

void Cnodeman::add(int n)
{
Cnode *T;
if (!first) {			  // if first 
			last=first=T=new Cnode;
			T->no=n;
			}
else					  // if not first 
        {
		T=first;
		while (T->next) T=T->next;
		last=T->next=new Cnode;
		T->next->no=n;
		}
}


Cnode* Cnodeman::getnode(int n)
{
Cnode* node;

node=first;
while (node) 
	{
	if (node->no==n) break;
	node=node->next;
	}

if (node==0) RE_MSG(2);
return node; // if requested node is not found 0 is returned
}

Cnode* Cnodeman::getnode(int n,double &rx,double &ry,double &rz)
{
Cnode* node;
node=first;
while (node)
{
if (node->no==n)
	{
	rx=node->x[0];
	ry=node->x[1];
	rz=node->x[2];
	break;
	}
node=node->next;
}
if (node==0) RE_MSG(2);
return node; // if requested node is not found 0 is returned
}

// Added on 25/09/96.
Cnode *Cnodeman::assign_elem_node(Celem *elem,Cnode* &node,int enode_no,int node_no)
{
// Search for structure node "node_no".
node=getnode(node_no);
// vks: if (!node) print error message
// "enode_no" is the node no. in single element node ( not structure node ).

//Compute all possible dof for common node between elems.
node->no_node_dof=max(node->no_node_dof,elem->no_enode_dof[enode_no]);
return node;
}


void Cnodeman::Read_Coords()
{
Cnode *node;
node=first;
while (node)
		{
			node->Read_Coords(inpfile);
			Generate_intermediate_nodes();
			node=node->next;
		}
}

int Cnodeman::Create_IDarray()
{
int node_no,i,j;
// Read ID array
Cnode* node=first;

// Allocate Memory and Initialize ID[no_node_dof]=0.
while (node)
//for(node_no=1;node_no<=no_nodes;++node_no)
		{
		//node=getnode(node->no);
		node->AllocInit_ID_memory();
		node=node->next;
		}

for(i=1;i<=no_rest_nodes;++i)
	{
	fscanf(inpfile,"%d",&node_no);
	//fprintf(outfile,"Node No. %d  ",node_no);
	node=getnode(node_no);
	
	for(j=0;j<node->no_node_dof;++j)
		{
		fscanf(inpfile,"%d",&node->ID[j]);
	//	fprintf(outfile,"%d ",node->ID[j]);

		// if -1 is input, arrest all subsequent nodal dof.
		if (node->ID[j]==-1) {
							 for(int k=j;k<node->no_node_dof;++k)
							 			node->ID[k]=1;
							 break;
							 }
		if (node->ID[j]>1) RE_MSG(1);
		}

	for(j=node->no_node_dof;j<MAX_NO_NODE_DOF;++j)      // added on 5/12/96
		{
		int pseudoID;
		fscanf(inpfile,"%d",&pseudoID);  // not used.
		}

	//fprintf(outfile,"\n");
	}

/* Reconstruction of ID array */
int	no_dof=0;

node=first;
while (node)
//	for(node_no=1;node_no<=no_nodes;++node_no)
	{
	//node=getnode(node_no);
	for(j=0;j<node->no_node_dof;++j)
		{
		if (node->ID[j]==0)
				node->ID[j]=++no_dof-1;
			else
				node->ID[j]=-1; // restraint direction
		}
	node=node->next;
	}

//fprintf(outfile,"ID matrix \n");
node=first;
while (node)
//for(node_no=1;node_no<=no_nodes;++node_no)
	{
////fprintf(outfile,"NODE %3d : ",node->no);
	if (node->no_node_dof==0) fprintf(outfile,"\n********Warning from Node Manager: Node %3d is unused!\n",node->no);
//	for(j=0;j<node->no_node_dof;++j)
//			fprintf(outfile,"%3d ",node->ID[j]);
//	fprintf(outfile,"\n");
	node=node->next;
	}

fprintf(outfile,"\nTotal no. of global dof. = %d \n",no_dof);
return (no_dof) ;
}

/*
void Cnodeman::Read_Loads()
{
Cnode *node;
int node_no;
for (int j=1;j<=no_loaded_nodes[0];j++)
		{
		fscanf(inpfile,"%d         ",&node_no);
		fprintf(outfile,"Node No. = %d  ",node_no);
		node=getnode(node_no);
	
	for (int i=0;i<node->no_node_dof;i++)       
			{
			fscanf(inpfile,"%lf",&AE[node->ID[i]]);
			fprintf(outfile,"%6.3lf ",AE[node->ID[i]]);
			}

	for(i=node->no_node_dof;i<MAX_NO_NODE_DOF;++j)      // added on 5/12/96
		{
		int pseudoLOAD;
		fscanf(inpfile,"%d",&pseudoLOAD);  // not used.
		}
	
		fprintf(outfile,"\n");
		}
}
*/

void Cnodeman::Generate_intermediate_nodes()
{
/*
	k=0;
	for (i=1;i<=no_elems_this;i++)
	{
		// Automatic generation of intermediate node numbers
		k++;
	 	while (this->no>k)
			{
			Celem* Telem=elemman->getelem(k);
			
			// Increase all elem node numbers by KI
			
			for(j=0;j<no_elem_nodes;j++)
			 	Telem->enode[j]->no=enode[j]->no+KI; 
			
			k++;
			if (elem_no>no_elems_this) break;
			}
	}
*/
}
//************** Cloadman definitions 

Cloadman::Cloadman()
{
first=0;
last=0;
}

Cloadman::~Cloadman()
{
if (first) delete first;
}


void Cloadman::add(int n)
{
Cload *T;
if (!first) {			  // if first 
			last=first=T=new Cload;
			T->no=n;
			}
else					  // if not first 
        {
		T=first;
		while (T->next) T=T->next;
		last=T->next=new Cload;
		T->next->no=n;
		}
}


Cload* Cloadman::getload(int n)
{
Cload* load;

load=first;
while (load) 
	{
	if (load->no==n) break;
	load=load->next;
	}

if (load==0) RE_MSG(2);
return load; // if requested load is not found 0 is returned
}

/*
Cload* Cloadman::getload(int n)
{
Cload* load;
load=first;
while (load)
{
if (load->no==n)	break;
load=load->next;
}
return load; // if requested load is not found 0 is returned
}
*/

// Added on 25/09/96.
Cload *Cloadman::assign_elem_load(Celem *elem,Cload* &load,int eload_no,int load_no)
{
// Search for structure load "load_no".
load=getload(load_no);
// "eload_no" is the load no. in single element load ( not structure load ).

return load;
}

void Cloadman::Read_Loads()
{
Cload *load;
load=first;
while (load)
		{
			load->Read_Loads(inpfile);
			load=load->next;
		}
}

// Added on 22/11/96
void Cloadman::Read_Multipliers() // Read element load multipliers
{

for (int i=0;i<10;i++) elem_LM[i]=1;

/*
for (i=0;i<no_load_cases;i++)
	{
fscanf(inpfile,"%lf",&elem_LM[i]);
fprintf(outfile,"Element Load Multipliers for load case %d = %lf ",i,elem_LM[i]);
	}
fprintf(outfile,"\n");

for (i=0;i<no_load_cases;i++)
	{
	fscanf(inpfile,"%d",&no_loaded_nodes[i]);
	fprintf(outfile,"No. of loaded nodes for load case %d = %d ",i,no_loaded_nodes[i]);
	}
fprintf(outfile,"\n");
*/
}


void Cloadman::Read_Qnodal_2_P(int load_case_no) // Read element load multipliers
{
int node_no;

for (int j=1;j<=no_loaded_nodes[load_case_no];j++)
		{
		fscanf(inpfile,"%d         ",&node_no);
		fprintf(outfile,"Node No. = %d  ",node_no);
		Cnode *node=nodeman->getnode(node_no);

		for (int i=0;i<node->no_node_dof;i++)       
			{
			fscanf(inpfile,"%lf",&AE[node->ID[i]]);
			fprintf(outfile,"%6.3lf ",AE[node->ID[i]]);
			
			// Modify AE is redundant . P and AE : resolve ambiguity.
			P[node->ID[i]]+=AE[node->ID[i]];
			}

		for (i=node->no_node_dof;i<MAX_NO_NODE_DOF;i++)       
			{
			double pseudoLOAD;
			fscanf(inpfile,"%lf",&pseudoLOAD);
			}

		fprintf(outfile,"\n");
		}
}



void Cloadman::Multiply_Q_toget_P() 	//  P{} = Q{} * elem_LM[load_case_no]
{
for (int i=0;i<no_dof;i++)
		P[i]=/*P[i]+*/Q[i]*elem_LM[i];
}


//************** Cnodeman definitions 

Cmaterialman::Cmaterialman()
{
first=0;
last=0;
}

Cmaterialman::~Cmaterialman()
{
if (first) delete first;
}

void Cmaterialman::add(int n)
{
Cmaterial *T;
if (!first) {			  // if first 
			last=first=T=new Cmaterial;
			T->no=n;
			}
else					  // if not first 
        {
		T=first;
		while (T->next) T=T->next;
		last=T->next=new Cmaterial;
		T->next->no=n;
		}
}


Cmaterial* Cmaterialman::getmaterial(int n)
{
Cmaterial* material;

material=first;
while (material) 
	{
	if (material->no==n) break;
	material=material->next;
	}

if (material==0) RE_MSG(2);
return material; // if requested material is not found 0 is returned
}

Cmaterial* Cmaterialman::getmaterial(int n,double &rE,double &rMEU,double &rWD)
{
Cmaterial* material;
material=first;
while (material)
{
if (material->no==n)
	{
	rE=material->E;
	rMEU=material->MEU;
	rWD=material->WD;
	break;
	}
material=material->next;
}
return material; // if requested material is not found 0 is returned
}

// Added on 25/09/96.
Cmaterial *Cmaterialman::assign_elem_material(Celem *elem,Cmaterial* &material,int ematerial_no,int material_no)
{
// Search for structure material "material_no".
material=getmaterial(material_no);
// "ematerial_no" is the material no. in single element material ( not structure material ).

return material;
}


void Cmaterialman::Read_Properties()
{
Cmaterial *material;
material=first;
while (material)
		{
			material->Read_Properties(inpfile);
			material=material->next;
		}
}

