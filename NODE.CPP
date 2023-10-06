#include "stdafx.h"
#include "space.h"


// ******** CNode definition
Cnode::Cnode()
{
no=0;
next=0;
no_node_dof=0;    // Value decided while connecting to elements.
}

Cnode::~Cnode()
{
if (ID) delete(ID);
if (next) delete(next);
}

void Cnode::AllocInit_ID_memory()
{
if (no_node_dof)
	{
	ID=new int[no_node_dof];
	for (int i=0;i<no_node_dof;i++) ID[i]=0;
	}
	else ID=0;
}


void Cnode::Read_Coords(FILE *inpfile)
{
	fscanf(inpfile,"%d",&no);
	fscanf(inpfile,"%lf %lf %lf\n",&x[0],&x[1],&x[2]); 
//	fprintf(outfile,"Node No. %3d  x = %6.3lf y = %6.3lf z = %6.3lf\n",no,x[0],x[1],x[2]);
}


BOOL Cnode::is_restrained()
{
for (int i=0;i<no_node_dof;i++) 
			if (ID[i]==-1) return 1;
return 0;
}	


// ****************************** Cload definition********************************

Cload::Cload()
{
no=0;
next=0;
}

Cload::~Cload()
{
if (next) delete(next);
}

void Cload::Read_Loads(FILE *inpfile)
{
	fscanf(inpfile,"%d %d %d %lf %lf \n",&no,&loadface,
							&loadtype,&magnitude,&zhead); 
	fprintf(outfile,"load No. %d  Load face = %d Load Type = %6.3lf "
	,"Magnitude = %6.3lf Head = %6.3lf\n",no,loadface,loadtype,magnitude,zhead);
}


 // no_dloads = no of distributed load sets

// Load_Type = 1  UDL 
//				  z_head = 0
// Load_Type = 2  Hydrostatically varying loads 
//				  Magnitude = weight density of fluid

// z_head = Global z cordinate of the surface of the fluid causing pressure loading
// loadface 	= face on which load acts


// ************************** Cmaterial definition********************************

Cmaterial::Cmaterial()
{
no=0;
next=0;
}

Cmaterial::~Cmaterial()
{
if (next) delete(next);
}

void Cmaterial::Read_Properties(FILE *inpfile)
{
// no. = Material no.,  E =Young's Modulas,  MEU = Poisson's ratio,  WD = Weight Density
fscanf(inpfile,"%d %lf %lf %lf\n",&no,&E,&MEU,&WD); 
fprintf(outfile,"Material No. %d  E = %6.3lf MEU = %6.3lf WD = %6.3lf\n",no,E,MEU,WD);
}

// Return material constants
void Cmaterial::get(double &pE,double &pMEU,double &pWD,double &pG)
{
pE=E;
pMEU=MEU;
pWD=WD;
pG=E/(2.0*(1.0+MEU));
}


