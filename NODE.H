#ifndef _INC_NODE_H
#define _INC_NODE_H

class Cnode {              // a single node
			public:
			Cnode();
			~Cnode();
			void AllocInit_ID_memory();
     		BOOL is_restrained();
			void Read_Coords(FILE *inpfile);

			int no;
		//	double x,y,z;   // remove;
			double x[3];
			int *ID;
			Cnode* next;
			int no_node_dof;
			};

class Cmaterial{              // a single material
			public:
			Cmaterial();
			~Cmaterial();
			virtual void Read_Properties(FILE *inpfile);
			void get(double &pE,double &pMEU,double &pWD,double &pG);
			int no;
			double E,MEU,WD;
			Cmaterial* next;
			};

class Cload  {              // a single load
			public:
			Cload();
			~Cload();
			void Read_Loads(FILE *inpfile);

			int no;
		   	
		   	int loadtype;
			double magnitude;
			double zhead;
			int loadface;

			Cload* next;
			};

#endif
