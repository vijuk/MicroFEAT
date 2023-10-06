#ifndef _INC_ELEM_H
#define _INC_ELEM_H	

#define 	ELEM3D_FUNCTIONS\
	public:\
        virtual void Compute_C(double &c1,double &c2, double &c3){};\
		virtual void Contribute_SMs(double WrstJ,double WrstJ_C,double cc1,double cc2,double cc3){};\
		virtual void Assemble_Q(){};\
		virtual void Compute_B_Matrix(){};\
		virtual void Compute_CB_at_face_centre(){};\
		virtual void Compute_Principal_Stresses(){};\
		virtual void Initialize_C_B_SMs_Qbody(){};\
		virtual void Contribute_Qbody(double WrstJ){};\
		virtual void Update_Colum_Heights(){};\
		virtual void Assemble_SK(){};\
		virtual void Compute_lmn_and_area(/*double JM[3][3],double e[3],double &Area*/){};\
		virtual void Compute_Qface(/*double force,double Gauss_Weight[],double Area,double lmn[3]*/){};\
		virtual void Compute_N(){};\
		virtual void Compute_DN_rst(){};\
		virtual void Compute_JM(/*double DN_rst[MAX_NO_ELEM_NODES][3]*/){};\
		virtual void Compute_Inv_JM(/*double JM[3][3]*/){};\
		virtual void Compute_DET_J(/*double JM[3][3],double Inv_JM[3][3]*/){};\
		virtual void Compute_DN_xyz(/*double Inv_JM[3][3],double DN_rst[MAX_NO_ELEM_NODES][3],double DET_J*/){};\
		virtual void Compute_Stresses(){};\
		virtual void Fill_Upper_SMs(){};\
		virtual void Compute_CB(){};\
		virtual void Compute_SMs_Qbody(double cc1,double cc2,double cc3){};\
		virtual double Compute_forces_at_GP(){return 0;};\
		virtual void Compute_lmn_and_area(int tx,int ty,int tz,double JM[3][3],double e[3],double &Area){};\
		virtual double Compute_WrstJ(int lx,int ly,int lz,double DET_J,int Gauss_Order,int icycle=0){return 0;}



//Elem 0
class Celem   				        // base class for element classes ( no = 0)
		{
		public:
		 Celem();
		 Celem(int eID,int ned,int nen,unsigned EDMA[]);
		~Celem();
		virtual void Assemble_K();
		virtual void Compute_SMm(double Axial_Load){};
		virtual void compute_MMm(){};
		virtual void compute_P(){};
		virtual void Read_Properties();
		virtual void Read_Loads();

		virtual void compute_R(){};
		virtual void Local_2_Global(double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
				double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]){};
		virtual void Change_2_Global(double *X){};
	    virtual void Change_2_Local(double *X){};
	
		virtual BOOL is_loaded(); 
		virtual void compute_gdir_vector(int gdir_for[]); /* on 9/9/96 */
		virtual int get_elem_sbw();
        virtual void NLA_Operations(int cnt,int count,double *disp);
		int elem_ID;
	 	int no;
		double L,A,Ix,Iy,Iz;	     // !move to child class

		double *AML,PreTension;	// Applied Member Loads and Pretension

    	Cnode **enode;		         // on 27/09/96
    	int  no_elem_dof,no_elem_nodes; /* on 9/9/96 */ 
		int * no_enode_dof;          // on 27/09/96 Store dof for member nodes.
		unsigned int *enode_dofmap;  // dof for element nodes in order.
	    Celem *next;

		Cmaterial *material; 	// on 21/11/96 Material Group No.
		double dx0,dy0,dz0;  // Used in non linear analysis

		ELEM3D_FUNCTIONS;

//		protected:
		
		void Init_Lower(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		void Fill_Upper(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		void Print_M(double A[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);

		};		

class Ctruss_3D :public Celem
		{
		public:
        Ctruss_3D(int eID,int ned,int nen,unsigned EDMA[]):Celem(eID,ned,nen,EDMA){};
		virtual void Read_Properties();
		virtual void Compute_SMm(double Axial_Load); // SMs formed directly here
		virtual void compute_MMm();
		virtual void compute_P();
		virtual void compute_R(){};  // not used
		virtual void Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]){}; 
		virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
		ELEM3D_FUNCTIONS;
		};


class Cgap_3D :public Ctruss_3D
		{
		public:
        Cgap_3D(int eID,int ned,int nen,unsigned EDMA[]):Ctruss_3D(eID,ned,nen,EDMA){};
		virtual void Read_Properties();
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm(){};
		virtual void compute_P(){};
		virtual void compute_R(){};
		virtual void Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]){};
		virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
		ELEM3D_FUNCTIONS;
		};


class Cbeam_3D : public Celem       // 3D beam element ( ID = 1)
		{
		public:
		Cbeam_3D();
		Cbeam_3D(int eID,int ned,int nen,unsigned EDMA[]);
		virtual void Read_Properties();
	
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm();
		virtual void compute_P();
		virtual void compute_R();
		virtual void Local_2_Global(
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
	    Cnode *knode;  	
		virtual void Change_2_Global(double *X);
		virtual void Change_2_Local(double *X);
        virtual void NLA_Operations(int cnt,int count,double *disp);

		ELEM3D_FUNCTIONS;
		
		private:

		void Include_SMmP(double Axial_Load);
//		static unsigned  int dofdec;	 
		};

class Cbeam_3D_curved: public Cbeam_3D  // 3D curved beam element ( ID = 2 )
		{
		public:
	      Cbeam_3D_curved(int eID,int ned,int nen,unsigned EDMA[]):Cbeam_3D(eID,ned,nen,EDMA){};
		  double r,IncAng;  // Radius & Included Angle of an arc element
		//virtual void Assemble_K();
		  virtual void Read_Properties();
		  virtual void Compute_SMm(double Axial_Load);
		  virtual void compute_R();
		  virtual void Local_2_Global(
		  	double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],
		  	double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
		  ELEM3D_FUNCTIONS;
		};

 class Cbeam_2D_VArch: public Cbeam_3D_curved  // 2D Vertical Arch beam element ( ID = 3)
		{
		public:
	      Cbeam_2D_VArch(int eID,int ned,int nen,unsigned EDMA[]):Cbeam_3D_curved(eID,ned,nen,EDMA){};
	//   ~beam_3D_curved();
//		  double       r,IncAng;  // Radius & Included Angle of an arc element
//		  virtual void Read_Properties();
// dof no. & dirn's are different >> seperate functions needed ( avoid !)
		 	//virtual void Assemble_K();
		 virtual void Compute_SMm(double Axial_Load);
 		  virtual void compute_R();
		  virtual void compute_P();
		  virtual void Local_2_Global(double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
		  ELEM3D_FUNCTIONS;
		};

 class Cbeam_2D_HArch: public Cbeam_2D_VArch  // 2D Horizontal Arch element ( ID =4 )
		{
		public:
		friend class  Ctriangle_Constant;

	      Cbeam_2D_HArch(int eID,int ned,int nen,unsigned EDMA[]):Cbeam_2D_VArch(eID,ned,nen,EDMA){};
	//	  beam_3D_curved();
	//   ~beam_3D_curved();
//		  double r,IncAng;  // Radius & Included Angle of an arc element
//		  virtual void Read_Properties();
		  //virtual void Assemble_K();
		  virtual void Compute_SMm(double Axial_Load);
//		  virtual void compute_R();
//    	  virtual void Local_2_Global(double (&const SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
		ELEM3D_FUNCTIONS;
		};
 class Cbeam_2D_HArch_12: public Cbeam_2D_VArch  // 2D Higher order Horizontal Arch element ( ID =5 )
		{
		public:
	      Cbeam_2D_HArch_12(int eID,int ned,int nen,unsigned EDMA[]):Cbeam_2D_VArch(eID,ned,nen,EDMA){};
	//	  beam_3D_curved();
	//   ~beam_3D_curved();
//		  double r,IncAng;  // Radius & Included Angle of an arc element
//		  virtual void Read_Properties();
		//virtual void Assemble_K();
		  virtual void Compute_SMm(double Axial_Load);
//		  virtual void compute_R();
 //   	  virtual void Local_2_Global(double (&const SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
 virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
			ELEM3D_FUNCTIONS;
		};

class Cbeam_3D_GArch: public Cbeam_3D_curved  // 3D General Arch element ( ID =6 )
		{
		public:
		Cbeam_3D_GArch(int eID,int ned,int nen,unsigned EDMA[]):Cbeam_3D_curved(eID,ned,nen,EDMA){};
	//	  //	  beam_3D_curved();
	//   ~beam_3D_curved();
	//	  double r,IncAng;  // Radius & Included Angle of an arc element
	//	  virtual void Read_Properties();
		//virtual void Assemble_K();
		  virtual void Compute_SMm(double Axial_Load);
	//	  virtual void compute_R();
	//	  virtual void Local_2_Global(double (&const SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
				  ELEM3D_FUNCTIONS;
	};

class Ctriangle_Constant : public Celem      // 3D beam element ( ID = 1)
		{
		public:
		Ctriangle_Constant();
		Ctriangle_Constant(int eID,int ned,int nen,unsigned EDMA[]):Celem(eID,ned,nen,EDMA){};
	
		//virtual void Assemble_K();
		virtual void Read_Properties();
		virtual void Read_Loads(){};		  
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm(){};
		virtual void compute_P();
		virtual void compute_R();		  // R= unity
		
		virtual void Local_2_Global(	  // remove
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);

virtual void NLA_Operations(int cnt,int count,double *disp){Celem::NLA_Operations(cnt,count,disp);};
				ELEM3D_FUNCTIONS;
		double Thickness;
		};

/*
class Cquad_iso4 : public Ctriangle_Constant
	{
	public :
		Cquad_iso4();
		Cquad_iso4(int eID,int ned,int nen,unsigned EDMA[]):Ctriangle_Constant(eID,ned,nen,EDMA){};
			virtual void Read_Properties();
			virtual void Assemble_K();
		virtual void Read_Loads(){};		  
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm(){};
		virtual void compute_P();
		virtual void compute_R();		  // R= unity
		virtual void Local_2_Global(	  // remove
		double (&const SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
};

*/

class Cquad8 : public Ctriangle_Constant
	{
	public :
		Cquad8();
		Cquad8(int eID,int ned,int nen,unsigned EDMA[]):Ctriangle_Constant(eID,ned,nen,EDMA){};
		virtual void Read_Properties();
		virtual void Assemble_K();
		virtual void Read_Loads(){};		  
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm(){};
		virtual void compute_P();
		virtual void compute_R();		  // R= unity
		virtual void Local_2_Global(	  // remove
		double (&SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		double S[100];

		virtual void Compute_C(double &c1,double &c2, double &c3);
		virtual void Contribute_SMs(double WrstJ,double WrstJ_C,double cc1,double cc2,double cc3);
		virtual void Assemble_Q();
		virtual void Compute_B_Matrix();
		virtual void Compute_CB_at_face_centre();
		virtual void Compute_Principal_Stresses();
		virtual void Initialize_C_B_SMs_Qbody(); 
		virtual void Contribute_Qbody(double WrstJ);
		virtual void Update_Colum_Heights();
		virtual void Assemble_SK();
		virtual void Compute_lmn_and_area(/*int tx,int ty,int tz,double JM[3][3],double e[3],double &Area*/);
		virtual void Compute_Qface(/*double force,double Gauss_Weight[],double Area,double lmn[3]*/);
		virtual void Compute_N();
		virtual void Compute_DN_rst();
		virtual void Compute_JM(/*double DN_rst[MAX_NO_ELEM_NODES][3]*/);
		virtual void Compute_Inv_JM(/*double JM[3][3]*/);				   
		virtual void Compute_DET_J(/*double JM[3][3],double Inv_JM[3][3]*/);
		virtual void Compute_DN_xyz(/*double Inv_JM[3][3],double DN_rst[MAX_NO_ELEM_NODES][3],double DET_J*/);
		virtual void Compute_Stresses();
		virtual void Fill_Upper_SMs();
		virtual void Compute_CB();
		virtual void Compute_SMs_Qbody(double cc1,double cc2,double cc3);
		virtual double Compute_forces_at_GP(/*double N[]*/);
		virtual double Compute_WrstJ(int lx,int ly,int lz,double DET_J,int Gauss_Order,int icycle=0);
private :
		void Compute_Stiffness_PartI(); // products of Ni'Nj' and W is assembled ino SMm[]
		void Compute_Stiffness_PartII();
};

/*
class CCsolid08 : public Celem
	{
	public :
		CCsolid08();
		CCsolid08(int eID,int ned,int nen,unsigned EDMA[]):Chexa_iso4(eID,ned,nen,EDMA){};
		virtual void Assemble_K();
		virtual void Read_Properties();
		virtual void Read_Loads(){};		  
		virtual void Compute_SMm(double Axial_Load);
		virtual void compute_MMm(){};

		virtual void compute_P();
		virtual void compute_R();		  
		virtual void Local_2_Global(	  
		double (&const SMm)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF],double (&SMs)[MAX_NO_ELEM_DOF][MAX_NO_ELEM_DOF]);
		void mult1();
		void mult2();
		void DMAT();
		void jacob();
		double S[100];
		int no_stress_faces;  	// no. of elem faces on which 
								// stress output is required.
		Cload *loads;
		int MatGroupID;
};
*/


class Csolid20 : public Celem
	{
	public :
		Csolid20();
		Csolid20(int eID,int ned,int nen,unsigned EDMA[]):Celem(eID,ned,nen,EDMA)
					{ eload[0]=0; eload[1]=0 ;};
		virtual void Read_Properties();
// 		virtual void compute_P();
		virtual void Compute_C(double &c1,double &c2, double &c3);
		virtual void Contribute_SMs(double WrstJ,double WrstJ_C,double cc1,double cc2,double cc3);
		virtual void Assemble_Q();
		virtual void Compute_B_Matrix();
		virtual void Compute_CB_at_face_centre();
		virtual void Compute_Principal_Stresses();
		virtual void Initialize_C_B_SMs_Qbody(); 
		virtual void Contribute_Qbody(double WrstJ);
		virtual void Update_Colum_Heights();
		virtual void Assemble_SK();
		virtual void Compute_lmn_and_area(/*int tx,int ty,int tz,double JM[3][3],double e[3],double &Area*/);
		virtual void Compute_Qface(/*double force,double Gauss_Weight[],double Area,double lmn[3]*/);
		virtual void Compute_N();
		virtual void Compute_DN_rst();
		virtual void Compute_JM(/*double DN_rst[MAX_NO_ELEM_NODES][3]*/);
		virtual void Compute_Inv_JM(/*double JM[3][3]*/);				   
		virtual void Compute_DET_J(/*double JM[3][3],double Inv_JM[3][3]*/);
		virtual void Compute_DN_xyz(/*double Inv_JM[3][3],double DN_rst[MAX_NO_ELEM_NODES][3],double DET_J*/);
		virtual void Compute_Stresses();
		virtual void Fill_Upper_SMs();
		virtual void Compute_CB();
		virtual void Compute_SMs_Qbody(double cc1,double cc2,double cc3);
		virtual double Compute_forces_at_GP(/*double N[]*/);
		virtual double Compute_WrstJ(int lx,int ly,int lz,double DET_J,int Gauss_Order,int icycle=0);

		int no_stress_faces;  	// no. of elem faces on which 
		Cload *eload[2];		// Two element loads only allowed now ?
 	};	

#endif 

// END 3D analysis
