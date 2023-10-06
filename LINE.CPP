#include <stdafx.h>
#include "line.h"

void CLine::OnDraw(CDC* pDC)
{
	CPoint	point;

    if (m_pt_array.GetSize() > 1)
	{
		point.x = (int)m_pt_array[0].x;
		point.y = (int)m_pt_array[0].y;

		pDC->MoveTo(point);
		for (int i=1; i< m_pt_array.GetSize(); i++)
		{
			point.x = (int)m_pt_array[i].x;
			point.y = (int)m_pt_array[i].y;
			pDC->LineTo(point);
		}
	}
}
