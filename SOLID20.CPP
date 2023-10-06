#include "stdafx.h"
#include "space.h"
 // Common functions

void Csolid20::Compute_C(double &c1,double &c2, double &c3)
	{

	double E=material->E;
	double MEU=material->MEU;

	double  CF=E/(1.0+MEU)/(1.0-2.0*MEU);
			C[0][0]=CF*(1.0-MEU);
			C[0][1]=CF*MEU;
			C[0][2]=C[0][1];
			C[1][0]=C[0][1];
			C[1][1]=C[0][0];
			C[1][2]=C[0][1];
			C[2][0]=C[0][1];
			C[2][1]=C[0][1];
			C[2][2]=C[0][0];
			C[3][3]=CF*(1.0-2.0*MEU)/2.0;
			C[4][4]=C[3][3];
			C[5][5]=C[3][3];

			c1=C[0][0];
			c2=C[0][1];
			c3=C[3][3];

/*
fprintf(outfile,"C Matrix \n");

	for (int i=0;i<6;i++)
		{
		for (int j=0;j<6;j++)
				fprintf(outfile,"%6.3lf ",C[i][j]);
		fprintf(outfile,"\n");
		}

*/
	}

// routine calculates the stress-displacement matrix CB at the centres
// of specified faces . called by THRD08 & THRD20
void Csolid20::Compute_CB_at_face_centre(/*int stress_Oface[7],double JM[3][3],double Inv_JM[3][3],
		double CB[7*6][60],int no_elem_dof*/) // [60] = no_elem_dof
{
int face_no;   // ??
double TC[6][60];
double e[3][3];
double tr[6][6];

int indexmap[6][2]={{0,1},{2,1},{2,0},{0,1},{2,1},{2,0}}; // why not [7][2]

int output_face=stress_Oface[face_no];

// Direction cosines l1,m1,n1 etc.  for x',y',z' axes 
// by normalising the vectors [ J11' J12' J13' ] etc.


// find "normalized" vector e[1] for face i  ref. Page 369
int index=indexmap[output_face][0];
double temp=sqrt(sqr(Inv_JM[0][index])+sqr(Inv_JM[1][index])+sqr(Inv_JM[2][index]));
// l(1)= e[1][1], m(1) = e[1][2] etc.

e[2][0]=Inv_JM[0][index]/temp;
e[2][1]=Inv_JM[1][index]/temp;
e[2][2]=Inv_JM[2][index]/temp;

// find "normalized" vector e[2] for face i  
index=indexmap[output_face][1];
temp=sqrt(sqr(JM[index][0])+sqr(JM[index][1])+sqr(JM[index][2]));

e[0][0]=JM[index][0]/temp;
e[0][1]=JM[index][1]/temp;
e[0][2]=JM[index][2]/temp;

// find "normalized" vector e[3] for face i  ; e[3] = e[1] x e[2]

e[1][0]=e[2][1]*e[0][2]-e[2][2]*e[0][1];
e[1][1]=e[2][2]*e[0][0]-e[2][0]*e[0][2];
e[1][2]=e[2][0]*e[0][1]-e[2][1]*e[0][0];

for(int i=0;i<3;i++)
for(int j=0;j<3;j++)
		tr[i][j]=sqr(e[i][i]);

tr[0][3]=e[0][0]*e[0][1]*2;
tr[0][4]=e[0][1]*e[0][2]*2;
tr[0][5]=e[0][0]*e[0][2]*2;

tr[1][3]=e[1][0]*e[1][1]*2;
tr[1][4]=e[1][1]*e[1][2]*2;
tr[1][5]=e[1][0]*e[1][2]*2;

tr[2][3]=e[2][0]*e[2][1]*2;
tr[2][4]=e[2][1]*e[2][2]*2;
tr[2][5]=e[2][0]*e[2][2]*2;

tr[3][0]=e[0][0]*e[1][0];
tr[3][1]=e[0][1]*e[1][1];
tr[3][2]=e[0][2]*e[1][2];

tr[3][3]=e[0][0]*e[1][1]+e[0][1]*e[1][0];
tr[3][4]=e[0][1]*e[1][2]+e[0][2]*e[1][1];
tr[3][5]=e[0][0]*e[1][2]+e[0][2]*e[1][0];

tr[4][0]=e[1][0]*e[2][0];
tr[4][1]=e[1][1]*e[2][1];
tr[4][2]=e[1][2]*e[2][2];

tr[4][3]=e[1][0]*e[2][1]+e[1][1]*e[1][0];
tr[4][4]=e[1][1]*e[2][2]+e[1][2]*e[2][1];
tr[4][5]=e[1][0]*e[2][2]+e[1][2]*e[2][0];

tr[5][0]=e[2][0]*e[0][0];
tr[5][1]=e[2][1]*e[0][1];
tr[5][2]=e[2][2]*e[0][2];

tr[5][3]=e[2][0]*e[0][1]+e[2][1]*e[0][0];
tr[5][4]=e[2][1]*e[0][2]+e[2][2]*e[0][1];
tr[5][5]=e[2][0]*e[0][2]+e[2][2]*e[0][0];

for(i=0;i<6;i++)
for(int j=0;j<no_elem_dof;j++)
	{
	TC[i][j]=0;
	for(int k=0;k<6;k++)
		TC[i][j]+=tr[i][k]*CB[6*(face_no-1)+k][j];
	}

for(i=0;i<6;i++)
	for(int j=0;j<no_elem_dof;j++)
		CB[6*(face_no-1)+i][j]=TC[i][j]; // CB for face centre 
}			

// routine computes the principal stresses at specified points
// called by thrd08 & thrd20

void Csolid20::Compute_Principal_Stresses(/*int no_stress_faces,int ns,int stress_Oface[7],
			double stress[42],double sp[21]*/)
{
double sg[6]; 
return;

for(int n=1;n<=no_stress_faces;n++)
{
int k=3*n-3;
int ii=k*2;
double aa,bb,cc;
double sp[999];
if (stress_Oface[n]!=0) 
	{
	aa=(stress[ii+1]+stress[ii+2])/2.0;
	bb=(stress[ii-1]+stress[ii-2])/2.0;
	cc= sqrt(sqr(bb)+sqr(stress[ii+4]));
	sp[k+1]=aa+cc;
	sp[k+2]=aa-cc;
	sp[k+3]=0;
	if (bb!=0) 
		sp[k+3]=28.648*atan2(stress[ii+4],bb);
	}

aa=(stress[ii+1]+stress[ii+2]+stress[ii+3])/3;

for(int i=1;i<=3;i++) 
	{
	sg[i]=stress[ii+1]-aa;
	sg[i+3]=stress[ii+i+3];
	}

double c1,c2,t,a;

c1=sg[1]*sg[1]+sg[2]*sg[2]+0.5*sg[3]*sg[3]+sg[4]*sg[4]+sg[5]*sg[5]+sg[6]*sg[6];
c2=sg[1]*(sg[2]*sg[3]-sg[5]*sg[5]+sg[4]*(sg[5]*sg[6]-sg[4]*sg[3])+sg[6]*sg[5]-sg[2]*sg[6]);
t=sqrt(c1/1.5);
a=c2*1.414214/(t*t*t);
if (a<-1) a=1;
if (a>1) a=1;
a=acos(a)/3;
t=t*1414214;

sp[k+1]=t*cos(a);
sp[k+2]=t*cos(a+2.0944);
sp[k+3]=t*cos(a-2.0944);

for(i=2;i<=3;i++) 
	if (sp[k+1]<=sp[k+i])
		{
		c2=sp[k+1];
		sp[k+1]=sp[k+i];
		sp[k+i]=c2;
		}

if (sp[k+2]>=sp[k+3])
	{
	c2=sp[k+2];
	sp[k+2]=sp[k+3];
	sp[k+3]=c2;
	for(i=1;i<=3;i++)
		sp[k+i]=sp[k+i]+aa;
	}
}
}



// computes load vector due to element applied loads 
void Csolid20::Compute_Qface(/*double force,double Gauss_Weight[],double Area,double lmn[3],double lx,double ly*/)
{
// Corner node numbers for 6 element faces
int corner_node[6][8]={	1,2,6,5,9,18,13,17,    // for Csolid20
						4,3,7,8,11,19,15,20,
						2,3,7,6,10,19,14,18,
						1,4,8,5,12,20,16,17,
						6,7,8,5,14,15,16,13,
						2,3,4,1,10,11,12,9
					   };

/*
int corner_node[6][4]=	{					   // for Csolid08
							{1,4,2,1},       
							{6,2,2,3},
							{3,4,7,3},
							{6,7,7,8},
							{8,4,5,8},
							{6,5,5,1}
						};

*/
int normal_axis[6]={1,1,2,2,3,3};			  


int no_face_nodes=8;  // for Csolid20
int Gauss_Order=3; // ??


for(int i=0;i<no_face_nodes*3;i++)	Qface[i]=0.0;

int nn;
for(int kk=0;kk<2;kk++)  // 2 = no_loads
if (eload[kk]!=NULL)
		{
		nn=eload[kk]->no;

		// Take 'nn'th load structure

		// integrate over the surface
		
		int	tz=normal_axis[eload[kk]->loadface];
		int tx=per[tz];
		int ty=per[tx];

		// Value of natural cordinate normal to the element face
		rst[tz]=((eload[kk]->loadface % 2) == 1) ? 1.0 : -1.0;


		for(int lx=0;lx<Gauss_Order;lx++)
			{
			rst[tx]=Gauss_Point[lx][Gauss_Order-1];
			for(int ly=0;ly<Gauss_Order;ly++)
				{
				rst[ty]=Gauss_Point[ly][Gauss_Order-1];
 				Compute_N();		// Shape functions at r,s,t
				Compute_DN_rst(); // Derivative of shape functions at r,s,t wrt r,s,t
 				Compute_JM(/*DN_rst*/);	// Jacobian at r,s,t
//				double e[3],Area;  // ??
				Compute_lmn_and_area(/*JM,e,Area,lx,ly*/);	// compute direction cosines normal to surface and area
				Compute_forces_at_GP(/*N*/); // compute pressure at GP
//				double force,lmn[3];

				Compute_Qface(/*force,Gauss_Weight,Area,lmn*/); // Find load components and store in Q_face
				}
		}

}
}





void Csolid20::Contribute_SMs(double WrstJ,double WrstJ_C,double cc1,double cc2,double cc3)  // bases 1>0
{
// no_elem_nodes= 8 or 11 for  Csolid08
// no_elem_nodes= 20  Csolid20

double uj,vj,wj,uu,vv,ww,uv,vu,uw,wu,vw,wv;


		// ADD CONTRIBUTION TO STIFFNESS MATRIX
				for(int i=0;i<no_elem_nodes;i++)  
					{
					int k1=3*i;
					int k2=k1+1;
					int k3=k2+1;

					double ui=DN_xyz[i][0];
					double vi=DN_xyz[i][1];
					double wi=DN_xyz[i][2];

					for(int j=0;j<no_elem_nodes;j++)
						{
						
						int m1=3*j;
						int m2=m1+1;
						int m3=m2+1;
						
						uj=DN_xyz[j][0];
						vj=DN_xyz[j][1];
						wj=DN_xyz[j][2];
						
						uu=ui*uj;
						vv=vi*vj;
						ww=wi*wj;
						uv=ui*vj;
						vu=vi*uj;
						uw=ui*wj;
						wu=wi*uj;
						vw=vi*wj;
						wv=wi*vj;
						
						double  c1,c2,c3;

/*						if (i>7 || j>7) // if incompatiable nodes included
							{
							c1=WrstJ_C*cc1;
							c2=WrstJ_C*cc2;
							c3=WrstJ_C*cc3;
							}
							else
*/							{
							c1=WrstJ*cc1;
							c2=WrstJ*cc2;
							c3=WrstJ*cc3;
							}
						
						SMs[k1][m1]+=c1*uu+c3*(vv+ww);
						SMs[k2][m2]+=c1*vv+c3*(ww+uu);
						SMs[k3][m3]+=c1*ww+c3*(uu+vv); 
						
						SMs[k1][m2]+=c2*uv+c3*vu;
						SMs[k1][m3]+=c2*uw+c3*wu; 

						SMs[k2][m3]+=c2*vw+c3*wv;

						if(i!=j)
							{
							SMs[k2][m1]+=c2*vu+c3*uv;
							SMs[k3][m1]+=c2*wu+c3*uw;
							SMs[k3][m2]+=c2*wv+c3*vw;
							}
						}
					}
}


void Csolid20::Assemble_Q() 
{
int gdir_for[MAX_NO_ELEM_DOF];
compute_gdir_vector(gdir_for);
for(int j=0;j<no_elem_dof;j++)
	{
		int kk=gdir_for[j];
		if (kk!=0) Q[kk]+=Qface[j]-Qbody[j];
	}
}

				void Csolid20::Compute_lmn_and_area(/*int tx,int ty,int tz,double JM[3][3],double e[3],double &Area*/) // 1 > 0
				{
					// compute direction cosines normal to surface and area
				int tx,ty;
				double e[3],Area;

					e[0]=(JM[tx][1]*JM[ty][2]-JM[tx][2]*JM[ty][1]);
					e[1]=(JM[tx][2]*JM[ty][0]-JM[tx][0]*JM[ty][2]);
					e[2]=(JM[tx][0]*JM[ty][1]-JM[tx][1]*JM[ty][0]);

					Area=sqrt(sqr(e[0])+sqr(e[1])+sqr(e[2]));
					if (Area<=0) fprintf(outfile,"Area non positive !\n");
					for (int i=0;i<3;i++) e[i]/=Area;
				}


				double Csolid20::Compute_forces_at_GP(/*double N[],Cload *load*/)
				{
					double force=0;
					Cload * load;
					// compute pressure at guass point
					if (load->loadtype==1/*UNIFORM_PRESSURE*/) force= load->magnitude;				// UDL
					if (load->loadtype==2/*HYDROSTATIC_PRESSURE*/)						
						{
						double zz=0;
						int no_face_nodes=8;				 
						for(int i=0;i<no_face_nodes;i++)
							zz+=N[i]*enode[i]->x[2];   // z cordinate
						zz=load->zhead-zz; // Find out head causing pressure at gauss point

						if (zz<0) 
								force=0;
							else
								force = load->magnitude*zz;
						}
				return (force);
				}

// functions for 20 noded elems
// computes the load vector due to element applied loads
// called by thrd20
/*void Csolid20::Compute_Qface(double force,double Gauss_Weight[],double Area,double lmn[3],int lx,int ly)
{
	double total=force*Gauss_Weight[lx]*Gauss_Weight[ly]*Area;
	// Load vectors for nodes on the face
	for(int i=0;i<no_face_nodes;i++)
		{
		int node_index=corner_node[loadface][i];
		double q=total*N[node_index];
		for (k=0;k<3;k++)
			Qface[3*ni+k]+=q*lmn[k];  // components in x,y,z directions
		}
}
*/
void Csolid20::Compute_N()
{

double RP=1.0+rst[0];
double RM=1.0-rst[0];
double SP=1.0+rst[1];
double SM=1.0-rst[1];
double TP=1.0+rst[2];
double TM=1.0-rst[2];

// SHAPE FUNCTIONS

 N[0]=RP*SM*TM*(RP+SM+TM-5.0)*0.125;
 N[1]=RP*SP*TM*(RP+SP+TM-5.0)*0.125;
 N[2]=RM*SP*TM*(RM+SP+TM-5.0)*0.125;
 N[3]=RM*SM*TM*(RM+SM+TM-5.0)*0.125;

 N[4]=RP*SM*TP*(RP+SM+TP-5.0)*0.125;
 N[5]=RP*SP*TP*(RP+SP+TP-5.0)*0.125;
 N[6]=RM*SP*TP*(RM+SP+TP-5.0)*0.125;
 N[7]=RM*SM*TP*(RM+SM+TP-5.0)*0.125;

 N[8]=RP*SP*SM*TM*0.25;
 N[9]=RP*RM*SP*TM*0.25;
N[10]=RM*SP*SM*TM*0.25;
N[11]=RP*RM*SM*TM*0.25;

N[12]=RP*SP*SM*TP*0.25;
N[13]=RP*RM*SP*TP*0.25;
N[14]=RM*SP*SM*TP*0.25;
N[15]=RP*RM*SM*TP*0.25;

N[16]=RP*SM*TP*TM*0.25;
N[17]=RP*SP*TP*TM*0.25;
N[18]=RM*SP*TP*TM*0.25;
N[19]=RM*SM*TP*TM*0.25;

//fprintf(outfile," %d :Shape function at rst = %10.3lf %10.3lf %10.3lf \n",no,rst[0],rst[1],rst[2]);
//for (int i=0; i<no_elem_nodes; i++) fprintf(outfile," %10.3lf ",N[i]);
}



void Csolid20::Compute_DN_rst()
{

double RP=1.0+rst[0];
double RM=1.0-rst[0];
double SP=1.0+rst[1];
double SM=1.0-rst[1];
double TP=1.0+rst[2];
double TM=1.0-rst[2];

 DN_rst[0][0]= SM*TM*(2.0*RP+SM+TM-5.0)*0.125;
 DN_rst[0][1]= SP*TM*(2.0*RP+SP+TM-5.0)*0.125;
 DN_rst[0][2]=-SP*TM*(2.0*RM+SP+TM-5.0)*0.125;
 DN_rst[0][3]=-SM*TM*(2.0*RM+SM+TM-5.0)*0.125;
 DN_rst[0][4]= SM*TP*(2.0*RP+SM+TP-5.0)*0.125;
 DN_rst[0][5]= SP*TP*(2.0*RP+SP+TP-5.0)*0.125;
 DN_rst[0][6]=-SP*TP*(2.0*RM+SP+TP-5.0)*0.125;
 DN_rst[0][7]=-SM*TP*(2.0*RM+SM+TP-5.0)*0.125;

 DN_rst[0][8]= SP*SM*TM*0.25;
 DN_rst[0][9]=-rst[0]*SP*TM*0.5;

DN_rst[0][10]=-DN_rst[0][8];
DN_rst[0][11]=-rst[0]*SM*TM*0.5;
DN_rst[0][12]= SP*SM*TP*0.25;
DN_rst[0][13]=-rst[0]*SP*TP*0.5;
DN_rst[0][14]=-DN_rst[0][12];
DN_rst[0][15]=-rst[0]*SM*TP*0.5;
DN_rst[0][16]=SM*TP*TM*0.25;
DN_rst[0][17]=SP*TP*TM*0.25;
DN_rst[0][18]=-DN_rst[0][17];
DN_rst[0][19]=-DN_rst[0][16];

DN_rst[1][0]=-RP*TM*(RP+2.0*SM+TM-5.0)*0.125;
DN_rst[1][1]= RP*TM*(RP+2.0*SP+TM-5.0)*0.125;
DN_rst[1][2]= RM*TM*(RM+2.0*SP+TM-5.0)*0.125;
DN_rst[1][3]=-RM*TM*(RM+2.0*SM+TM-5.0)*0.125;
DN_rst[1][4]=-RP*TP*(RP+2.0*SM+TP-5.0)*0.125;
DN_rst[1][5]= RP*TP*(RP+2.0*SP+TP-5.0)*0.125;
DN_rst[1][6]= RM*TP*(RM+2.0*SP+TP-5.0)*0.125;
DN_rst[1][7]=-RM*TP*(RM+2.0*SM+TP-5.0)*0.125;

DN_rst[1][8]=-RP*rst[1]*TM   *0.5;
DN_rst[1][9]= RP*RM*TM       *0.25;

DN_rst[1][10]=-RM*rst[1]*TM  *0.5;
DN_rst[1][11]=-DN_rst[1][9];
DN_rst[1][12]=-RP*rst[1]*TP  *0.5;
DN_rst[1][13]= RP*RM*TP      *0.25;
DN_rst[1][14]=-RM*rst[1]*TP  *0.5;
DN_rst[1][15]=-DN_rst[1][13];

DN_rst[1][16]=-RP*TP*TM      *0.25;
DN_rst[1][17]=-DN_rst[1][16];
DN_rst[1][18]= RM*TP*TM      *0.25;
DN_rst[1][19]=-DN_rst[1][18];

DN_rst[2][0]=-RP*SM*(RP+SM+2.0*TM-5.0)*0.125;
DN_rst[2][1]=-RP*SP*(RP+SP+2.0*TM-5.0)*0.125;
DN_rst[2][2]=-RM*SP*(RM+SP+2.0*TM-5.0)*0.125;
DN_rst[2][3]=-RM*SM*(RM+SM+2.0*TM-5.0)*0.125;

DN_rst[2][4]= RP*SM*(RP+SM+2.0*TP-5.0)*0.125;
DN_rst[2][5]= RP*SP*(RP+SP+2.0*TP-5.0)*0.125;
DN_rst[2][6]= RM*SP*(RM+SP+2.0*TP-5.0)*0.125;
DN_rst[2][7]= RM*SM*(RM+SM+2.0*TP-5.0)*0.125;

 DN_rst[2][8]=-RP*SP*SM*0.25;
 DN_rst[2][9]=-RP*RM*SP*0.25;
DN_rst[2][10]=-RM*SP*SM*0.25;
DN_rst[2][11]=-RP*RM*SM*0.25;

DN_rst[2][12]=-DN_rst[2][8];
DN_rst[2][13]=-DN_rst[2][9];
DN_rst[2][14]=-DN_rst[2][10];
DN_rst[2][15]=-DN_rst[2][11];

DN_rst[2][16]=-RP*SM*rst[2]*0.5;
DN_rst[2][17]=-RP*SP*rst[2]*0.5;
DN_rst[2][18]=-RM*SP*rst[2]*0.5;
DN_rst[2][19]=-RM*SM*rst[2]*0.5;

//fprintf(outfile,"\n %d :Shape DERIVATIVES at rst = %10.3lf %10.3lf %10.3lf \n",no,rst[0],rst[1],rst[2]);
//for (int i=0; i<no_elem_nodes; i++) 
//	fprintf(outfile," %10.3lf %10.3lf %10.3lf \n",DN_rst[0][i],DN_rst[1][i],DN_rst[2][i]);

}

void Csolid20::Compute_JM()
{
// JACOBIAN MATRIX JM
for (int i=0;i<3;i++)
for (int j=0;j<3;j++)
	{
	JM[i][j]=0.0;
	for (int k=0;k<no_elem_nodes;k++)
		JM[i][j]+=DN_rst[i][k]*enode[k]->x[j];
	}
}

void Csolid20::Compute_Inv_JM()		// verified/2
{
// Invert Jacobian Matrix
for (int i=0;i<3;i++)
{
int j=per[i];
int k=per[j];
Inv_JM[i][i]= JM[j][j]*JM[k][k]-JM[k][j]*JM[j][k];
Inv_JM[i][j]= JM[k][j]*JM[i][k]-JM[i][j]*JM[k][k];
Inv_JM[j][i]= JM[j][k]*JM[k][i]-JM[j][i]*JM[k][k];
}
}

void Csolid20::Compute_DET_J(/*JM,Inv_JM*/) // verified/1
{
DET_J=JM[0][0]*Inv_JM[0][0]+JM[0][1]*Inv_JM[1][0]+JM[0][2]*Inv_JM[2][0];
//fprintf(outfile,"ELEM %d |JM| at rst (%6.3lf %6.3lf %6.3lf) = %10.4lf\n",no,rst[0],rst[1],rst[2],DET_J);
if (DET_J<=0) fprintf(outfile,"Elem No %d :Determinant Jacobian is non positive!\n",no);
}

void Csolid20::Compute_DN_xyz(/*Inv_JM,DN_srt,DET_J*/)
{
// Matrix of xyz derivatives of shape functions
//fprintf(outfile,"DN_xyz ******************* %d\n",no);
for (int i=0;i<3;i++)
for (int j=0;j<no_elem_nodes;j++)
	{
	DN_xyz[j][i]=0;
	for (int k=0;k<3;k++)
		DN_xyz[j][i]+=Inv_JM[i][k]*DN_rst[k][j];
	DN_xyz[j][i]/=DET_J;
	//fprintf(outfile,"DN_xyz [%d][%d] = %lf \n",j,i,DN_xyz[j][i]);
	}
}

void Csolid20::Compute_Stresses()
{
/*
for(int elem_no=1;elem_no<=no_elems_this;elem_no++)
		scanf("",ID[l],l=1,60,Gauss_Order,MatGroup,(xx(n,1),xx(n,2),xx(n,3),elem_node_numbers[n],n=1,20)
		,Load_No[1],Load_No[2],no_stress_faces,(stress_Oface[j],j=1,no_stress_faces);
	scanf("",CB[i][j],i=1,ns),j=1,ns)
	/
*/
//	int no_stress_faces;
	for (int i=0;i<no_stress_faces*6;i++) stress[i]=0;
	int ii=(no_stress_faces-1)*6;
	int gdir_for[MAX_NO_ELEM_DOF];
	compute_gdir_vector(gdir_for);

	for(i=0;i<no_elem_dof;i++) 
		{
		int k=gdir_for[i];
		if (k!=0)
			{
			for(int l=0;l<6*no_stress_faces;l++) 
				stress[l]+=CB[l][i]*P[k];
			}
//		double sps;
		Compute_Principal_Stresses(/*no_stress_faces,ns,stress_Oface,stress,sps*/);
		int	m=7;
		int ki=1;
		do{
		int	m5=m+5;
		int	m1=(m+1)/2;
		int	m2=(m+5)/2;
		i++;
		int 	num=stress_Oface[ki];
		m=m+6;
		}
		while (m!=6*no_stress_faces);
	}
}

void Csolid20::Compute_B_Matrix()
{
double e[3]; // ?
    	// compute B matrix
		for(int l=0;l<no_stress_faces;l++)
			{
			int ll=stress_Oface[l];
			for  (int i=0;i<3;i++)	e[i]=face_centre_rst[ll][i];

Compute_DN_rst(/*rst*/);  // Derivative of shape functions at r,s,t wrt r,s,t
Compute_JM(/*DN_rst*/);	   // Jacobian at r,s,t
Compute_Inv_JM(/*JM*/);		  // Inverse Jacobain at r,s,t
Compute_DET_J(/*JM,Inv_JM*/);  // Determinant Jacobian at r,s,t
Compute_DN_xyz(/*Inv_JM,DN_srt,DET_J*/); // Derivative of shape function wrt. xyz.

			int 	l3=6*l-1;
			for(int k=0;k<no_elem_nodes;k++)
				{
				int k1=3*k;
				int k2=k1+1;
				int k3=k2+1;

				B[l3+1][k1]=DN_xyz[k][0];
				B[l3+2][k2]=DN_xyz[k][1];
				B[l3+3][k3]=DN_xyz[k][2];
				B[l3+4][k1]=DN_xyz[k][1];
				B[l3+4][k2]=DN_xyz[k][0];

				B[l3+5][k2]=DN_xyz[k][2];
				B[l3+5][k3]=DN_xyz[k][1];
				B[l3+6][k1]=DN_xyz[k][2];
				B[l3+6][k3]=DN_xyz[k][0];
				}
			}
}


    void Csolid20::Initialize_C_B_SMs_Qbody()
	{
	for(int i=0;i<6;i++)
	for(int j=0;j<6;j++)
		C[i][j]=0;
	for(i=0;i<no_elem_dof;i++)
			{
			Q[i]=0;
			for(int j=0;j<no_elem_dof;j++)
				SMs[i][j]=0;
			}
	for(i=0;i<no_stress_faces;i++)
	for(int j=0;j<no_elem_dof;j++)
			B[i][j]=0;

	}


void Csolid20::Compute_CB()
{
			for (int i=0;i<no_stress_faces;i++)
				{
				int ii=6*i;
				for (int j=0;j<6;j++)
				for (int k=0;k<no_elem_dof;k++)
					{
					CB[ii+j][k]=0.0;
					for (int l=0;l<6;l++)
						CB[ii+j][k]+=(C[j][l]*B[ii+l][k]);
					}

				int k;
				if ((k=stress_Oface[i])>0)
					{
					for (int j=0;j<3;j++) rst[j]=face_centre_rst[k][j];
Compute_DN_rst(/*rst*/);  // Derivative of shape functions at r,s,t wrt r,s,t
Compute_JM(/*DN_rst*/);	   // Jacobian at r,s,t
Compute_Inv_JM(/*JM*/);		  // Inverse Jacobain at r,s,t
Compute_CB_at_face_centre(/*stress_Oface,JM,Inv_JM,CB,i,60*/);
					}
				}
}

void Csolid20::Fill_Upper_SMs() // Fill Lower SMs
{
for (int i=0;i<no_elem_dof-1;i++)
	for (int j=i+1;j<no_elem_dof;j++)
		SMs[j][i]=SMs[i][j];   
}



void Csolid20::Contribute_Qbody(double WrstJ)
{
// Gravity Load Vector
double temp=WrstJ*material->WD;
for(int i=1;i<=no_elem_nodes;i++)	Qbody[3*i-1]+=temp*N[i-1];   // 8 or 20
}



// General functions for Solids

void Csolid20::Update_Colum_Heights()
{
// ColHeight = no. of elements in each column below the skyline 
// and above the diagonal element

int gdir_for[MAX_NO_ELEM_DOF];

int LS=100000;  // a LARGE number

compute_gdir_vector(gdir_for); 

// Find the least active global dof no. for the element
for (int i=0;i<no_elem_dof;i++)
		{
		int j=gdir_for[i];
		if ((j>-1) && (j<LS)) LS=j;
		}

//fprintf(outfile,"Updated ColHeight[] for elem: %d \n",no);
for (i=0;i<no_elem_dof;i++)
	{
	int j=gdir_for[i];
	if ((j>-1) && (j-LS>ColHeight[j])) ColHeight[j]=j-LS;
	//fprintf(outfile,"%6d ",ColHeight[j]);
	}

//fprintf(outfile,"\n");
}

void Csolid20::Assemble_SK()
{
int gdir_for[MAX_NO_ELEM_DOF];
compute_gdir_vector(gdir_for);
for(int i=0;i<no_elem_dof;i++)
	{
	int ii=gdir_for[i];
	if (ii>=0) 
		{
		for(int j=0;j<no_elem_dof;j++)
			{
			int jj=gdir_for[j];
			if (jj>=0) 
				{
				int mi=DiagDOF[jj];
				int ij = jj- ii ;
				if (ij>=0)
					SK[mi+ij]+=SMs[i][j];
				}
			}
		}
	}
}

double Csolid20::Compute_WrstJ(int lx,int ly,int lz,double DET_J,int Gauss_Order,int icycle)  // icycle=0 default
{
double Wrst;
if (Gauss_Order==2) 
		{
		if (icycle>8)	
				Wrst=0.8864265927977;
			else 
				Wrst=Gauss_Weight[0][1];
		}
		else
	Wrst=Gauss_Weight[lx][Gauss_Order-1]*Gauss_Weight[ly][Gauss_Order-1]*Gauss_Weight[lz][Gauss_Order-1];
return (Wrst*DET_J);
}


void Csolid20::Compute_SMs_Qbody(double cc1,double cc2,double cc3)
{

double xw=-0.7958224257542;
int icycle=0;
int mx=0;
int Gauss_Order=3 ; //??
for(int lx=0;lx<Gauss_Order;lx++)
for(int ly=0;ly<Gauss_Order;ly++)
for(int lz=0;lz<Gauss_Order;lz++)
	{
		rst[0]=Gauss_Point[lx][Gauss_Order-1];
		rst[1]=Gauss_Point[ly][Gauss_Order-1];
		rst[2]=Gauss_Point[lz][Gauss_Order-1];

		icycle++;
		BOOL flag=0;

		do {
			if (flag)
				{
				mx++; icycle++;
				for (int i=0;i<3;i++) rst[i]=0;
				switch (mx)
					{
					case 1: rst[0]= xw;
					case 2: rst[0]=-xw;
					case 3: rst[1]= xw;
					case 4: rst[1]=-xw;
					case 5: rst[2]=-xw;
					case 6: rst[2]= xw;
					}
				}


		Compute_N();/*rst*/		  // Shape functions at r,s,t
		Compute_DN_rst(/*rst*/);  // Derivative of shape functions at r,s,t wrt r,s,t
		Compute_JM(/*DN_rst*/);	   // Jacobian at r,s,t
		Compute_Inv_JM(/*JM*/);		  // Inverse Jacobain at r,s,t
		Compute_DET_J(/*JM,Inv_JM*/);  // Determinant Jacobian at r,s,t
		Compute_DN_xyz(/*Inv_JM,DN_srt,DET_J*/); // Derivative of shape function wrt. xyz.
		double WrstJ=Compute_WrstJ(lx,ly,lz,DET_J,Gauss_Order,icycle);
		double WrstJ_C;
		Contribute_Qbody(WrstJ);

		Contribute_SMs(WrstJ,WrstJ_C,cc1,cc2,cc3);
		if ((Gauss_Order==2) && (icycle>=8) && (icycle<14)) flag=1; else break;
		} while (1);
	}	

Fill_Upper_SMs();  // Actually fills Lower SMs
/*fprintf(outfile," STIFFNESS FOR ELEM (%d) \n",no);
for (int i=0;i<no_elem_dof;i++)
	{
	for (int j=0;j<no_elem_dof/4;j++)
			fprintf(outfile," %12.2lf  ",SMs[i][j]);
	fprintf(outfile," Q = %6.3lf\n",Qbody[i]);
	}
*/
}

/*
void Csolid20::Compute_SMs_Qbody()
{

double rst[3]; // ??
for(int lx=1;lx<=2;lx++)
		{
		rst[1]=Gauss_Point[lx];
		for(int ly=1;ly<=2;ly++)
			{
			rst[2]=Gauss_Point[ly];
			for(int lz=1;lz<=2;lz++)
				{
				rst[3]=Gauss_Point[lz];

//Compute_N(/*rst*);		// Shape functions at r,s,t
//Compute_DN_rst(/*rst*); // Derivative of shape functions at r,s,t wrt r,s,t
//Compute_JM(/*DN_rst*);	// Jacobian at r,s,t
//Compute_Inv_JM(/*JM*);	// Inverse Jacobain at r,s,t
//Compute_DET_J(/*JM,Inv_JM*);	// Determinant Jacobian at r,s,t

//Compute_JM_C();		// Jacobian for incompatiable modes at centre
//Compute_Inv_JM_C(JM_C);	// Inverse Jacobian for incompatiable modes at centre
//Compute_DET_JM_C(JM_C,Inv_JM_C); // Determinat JAcobain for incompatiable modes at centre

//Compute_DN_xyz(/*Inv_JM,DN_srt,DET_J*); // Derivative of shape function wrt. xyz.
//double WrstJ   =Gauss_Weight[lx]*Gauss_Weight[ly]*Gauss_Weight[lz]*DET_J;

//double WrstJ_C =Gauss_Weight[lx]*Gauss_Weight[ly]*Gauss_Weight[lz]*DET_J_C;
	
//Compute_Q_body();
//Contribute_SMs();
//}
//Condense_SMs_B();  // Condense incompatiable modes
//}
*/

void Csolid20::Read_Properties()	// Modified on 25/09/96.
{
Celem::Read_Properties();

// Gauss Order = Numerical Integration order
int Gauss_Order=3;  // ToDo: Read from user for generality

	// KI = element data generator
	// Load_No(1),, Load_No(2) = load set id. no.'s for udl acting on this element
	// no_stress_faces =  no. of faces on which stresses are to be calculated.
	// incomp_mode =1 inclusion of incompatible modes incomp_mode=0 for suppression

int Load_No[2];// ??
int KI;
fscanf(inpfile,"%d %d %d %d",&Load_No[0],&Load_No[1],&KI,&no_stress_faces);
/*fprintf(outfile,"Load 1 = %d Load 2 = %d KI = %d no_stress_faces = %d \n",
					Load_No[0],Load_No[1],KI,no_stress_faces);
*/
for(int j=0;j<no_stress_faces;j++) scanf("%d",&stress_Oface[j]);	

//incomp_mode = INCOMP_MODE; // ToDo: Read from user for generality
}

