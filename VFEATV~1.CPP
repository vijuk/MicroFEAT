// vfeatView.cpp : implementation of the CVfeatView class
//

#include "stdafx.h"
#include "vfeat.h"

#include "vfeatDoc.h"
#include "vfeatView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CVfeatView

IMPLEMENT_DYNCREATE(CVfeatView, CView)

BEGIN_MESSAGE_MAP(CVfeatView, CView)
	//{{AFX_MSG_MAP(CVfeatView)
	ON_COMMAND(VA_ID_ANALYSE, OnAnalyse)
 	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_COMMAND(VM_ID_LINE, OnLine)
	ON_COMMAND(VM_ID_POLYLINE, OnPolyline)
	ON_COMMAND(VM_ID_RECTANGLE, OnRectangle)
	ON_COMMAND(VM_ID_TEXT, OnText)
	ON_COMMAND(VM_ID_ELLIPSE, OnEllipse)
	ON_WM_LBUTTONDBLCLK()
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CVfeatView construction/destruction

CVfeatView::CVfeatView()
{
	// TODO: add construction code here
	// m_pt_array.SetSize(100);
}

CVfeatView::~CVfeatView()
{
}

BOOL CVfeatView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CVfeatView drawing

void CVfeatView::OnDraw(CDC* pDC)
{
	CVfeatDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

	// TODO: add draw code for native data here
	
	pDoc->OnDraw(pDC);
	
}

/////////////////////////////////////////////////////////////////////////////
// CVfeatView printing

BOOL CVfeatView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CVfeatView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CVfeatView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

/////////////////////////////////////////////////////////////////////////////
// CVfeatView diagnostics

#ifdef _DEBUG
void CVfeatView::AssertValid() const
{
	CView::AssertValid();
}

void CVfeatView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CVfeatDoc* CVfeatView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CVfeatDoc)));
	return (CVfeatDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CVfeatView message handlers


void CVfeatView::OnAnalyse() 
{
	// TODO: Add your command handler code here
	
}



static BOOL	_mousePressed;
static BOOL	_mouseMoved;

CPoint GetScreenPoint(
		VmPoint3D	pt
		)
{
	CPoint	spt;
	spt.x = int(pt.x);
	spt.y = int(pt.y);
	return spt;
}

void CVfeatView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	VmPoint3D	pt3;
	// TODO: Add your message handler code here and/or call default
	_mousePressed = TRUE;	
	m_press_point = point;

	pt3.x = point.x;
	pt3.y = point.y;
	pt3.z = 0;
	_mouseMoved = FALSE;


	m_pt_array.Add(pt3);
	CClientDC	dc(this);

	int count = m_pt_array.GetSize();
	for (int i=0;i<count-1;i++)
	{
		dc.MoveTo(GetScreenPoint(m_pt_array[i]));
	    dc.LineTo(GetScreenPoint(m_pt_array[i+1]));
	}

	CView::OnLButtonDown(nFlags, point);
}

void CVfeatView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	VmPoint3D	pt3;

	// TODO: Add your message handler code here and/or call default
	_mousePressed = FALSE;	
	CView::OnLButtonUp(nFlags, point);

	CVfeatDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);

 	{
		pt3.x = point.x;
		pt3.y = point.y;
		pt3.z = 0.0;

		m_pt_array.Add(pt3);
	}
	
    if (m_context == VmContextLINE)
	{
		pDoc->AddLine(m_pt_array);
		m_pt_array.RemoveAll();
	}
}

void CVfeatView::OnMouseMove(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	VmPoint3D			pt3;
	static	CPoint		move_pt;

    if (_mousePressed)
	{
		CClientDC	dc(this);

		dc.SetROP2(R2_NOTMERGEPEN);
		if (!_mouseMoved)
		{
			_mouseMoved = TRUE;
		}
		else
		{
			
			dc.MoveTo(m_press_point);
			dc.LineTo(move_pt);
		}
		
		dc.MoveTo(m_press_point);
		dc.LineTo(point);

		move_pt = point;
		pt3.x = point.x;
		pt3.y = point.y;
		pt3.z = 0;

	    if (m_context == VmContextFREEHAND)
		{
			m_pt_array.Add(pt3);
			m_press_point = point;
		}
	}
	CView::OnMouseMove(nFlags, point);
}

void CVfeatView::OnLine() 
{
	// TODO: Add your command handler code here
	StartContext(VmContextLINE);
}

void CVfeatView::OnPolyline() 
{
	// TODO: Add your command handler code here
	StartContext(VmContextPOLYLINE);
}

void CVfeatView::OnRectangle() 
{
	// TODO: Add your command handler code here
	
}

void CVfeatView::OnText() 
{
	// TODO: Add your command handler code here
	
}

void CVfeatView::OnEllipse() 
{
	// TODO: Add your command handler code here
	
}

void CVfeatView::StartContext(int	context)
{
	m_context = context;
}

void CVfeatView::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	CVfeatDoc* pDoc = GetDocument();
	pDoc->AddLine(m_pt_array);
	m_pt_array.RemoveAll();

	CView::OnLButtonDblClk(nFlags, point);
}
