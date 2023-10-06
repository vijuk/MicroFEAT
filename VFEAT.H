// vfeat.h : main header file for the VFEAT application
//

#if !defined(AFX_VFEAT_H__A427EE26_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
#define AFX_VFEAT_H__A427EE26_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CVfeatApp:
// See vfeat.cpp for the implementation of this class
//

class CVfeatApp : public CWinApp
{
public:
	CVfeatApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CVfeatApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CVfeatApp)
	afx_msg void OnAppAbout();
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_VFEAT_H__A427EE26_E6E3_11D1_A4D1_E898DCB8424A__INCLUDED_)
