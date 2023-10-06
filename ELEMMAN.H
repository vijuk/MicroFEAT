#ifndef _INC_ELEMMAN_H
#define _INC_ELEMMAN_H
class Celemman		   // element manager 
		{
		public:
		Celemman();
		~Celemman();
		void add(int n,Celem* new_elem);
		void add_elem_group(int elem_group_no,int no_elems);
		void Read_Loads();
		void Read_Properties();
		int  compute_SBW();
		void Compute_ColHeight();
		void Post_Analysis();
		void SetUp_SkySystem();
		void SetUp_BandedSystem();
		Celem* getelem(int n);
		Celem *first,*last;
		};
#endif

