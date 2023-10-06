// vfeatView.h : interface of the CVfeatView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_VFEATVIEW_H__A427EE30_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
#define AFX_VFEATVIEW_H__A427EE30_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
#include <afxtempl.h>

typedef enum 
{
	VmContextDEFAULT,
	VmContextLINE,
	VmContextPOLYLINE,
	VmContextRECTANGLE,
	VmContextELLIPSE,
	VmContextTEXT,
	VmContextFREEHAND,
	VmContextSELECT,
	VmContextABORT
} VmContext;


class CVfeatView : public CView
{
protected: // create from serialization only
	CVfeatView();
	
	CPoint							m_press_point;
	CPoint							m_release_point;

	int								m_context;
	CArray <VmPoint3D, VmPoint3D &>	m_pt_array;


	DECLARE_DYNCREATE(CVfeatView)

// Attributes
public:
	CVfeatDoc* GetDocument();

// Operations
public:
	void StartContext(int	context);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CVfeatView)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CVfeatView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CVfeatView)
	afx_msg void OnAnalyse();
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLine();
	afx_msg void OnPolyline();
	afx_msg void OnRectangle();
	afx_msg void OnText();
	afx_msg void OnEllipse();
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in vfeatView.cpp
inline CVfeatDoc* CVfeatView::GetDocument()
   { return (CVfeatDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_VFEATVIEW_H__A427EE30_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
