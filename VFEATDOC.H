// vfeatDoc.h : interface of the CVfeatDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_VFEATDOC_H__A427EE2E_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
#define AFX_VFEATDOC_H__A427EE2E_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "stdafx.h"
#include "line.h"

class CVfeatDoc : public CDocument
{
protected: // create from serialization only
	CVfeatDoc();
	DECLARE_DYNCREATE(CVfeatDoc)

// Attributes
public:
	CArray <CLine, CLine &>	m_lines;
// Operations
public:
	void	AddLine(CArray <VmPoint3D, VmPoint3D &> & line);
	void OnDraw(CDC* pDC);

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CVfeatDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CVfeatDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CVfeatDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_VFEATDOC_H__A427EE2E_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
