#ifndef _INC_LINE_H
#define _INC_LINE_H
#include <afxtempl.h>
#include "stdafx.h"
#include "vm.h"

class CLine
{   
	// data variables
	public:
	CArray <VmPoint3D, VmPoint3D&>	m_pt_array;

	// a single node
	public:
		CLine(){};
		CLine(CArray <VmPoint3D, VmPoint3D &> pt_array)
		{
			m_pt_array.Copy(pt_array);
		}

		CLine & operator = (CLine & line)
		{
			this->m_pt_array.Copy(line.m_pt_array);
			return *this;
		}

		~CLine(){};
		void OnDraw(CDC* pDC);
};

#endif
