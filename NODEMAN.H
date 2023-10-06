#ifndef _INC_NODEMAN_H
#define _INC_NODEMAN_H
#include "space.h"
class Cnodeman     // node manager 
		{
		public:
		Cnodeman();
		~Cnodeman();

		void add(int n);

		Cnode *first,*last;

        Cnode* getnode(int n,double &rx,double &ry,double &rz);
		Cnode* getnode(int n);
        Cnode* assign_elem_node(Celem *elem,Cnode * &node,int enode_no,int node_no);
		void Read_Coords();
		int  Create_IDarray();
	//	void Read_Loads();		
		void Generate_intermediate_nodes(); // added on 22/11/96		 	
	};

class Cmaterialman    // node manager 
	{
	public:
	Cmaterialman();
	~Cmaterialman();
	Cmaterial *first,*last;
	void Read_Properties();
	void add(int n);					// ??
	Cmaterial *getmaterial(int n);
Cmaterial* Cmaterialman::getmaterial(int n,double &rE,double &rMEU,double &rWD);
	Cmaterial *assign_elem_material(Celem *elem,Cmaterial* &material,int ematerial_no,int material_no);
	};

class Cloadman      // node manager 
	{
	public:
	Cloadman();
	~Cloadman();

	Cload *first,*last;
	void add(int n);               // ?? 
	Cload* getload(int n);
	Cload * assign_elem_load(Celem *elem,Cload* &load,int eload_no,int load_no);
	void Read_Loads();
	void Read_Multipliers();
	void Read_Qnodal_2_P(int load_case_no);
	void Multiply_Q_toget_P();
	};
#endif

