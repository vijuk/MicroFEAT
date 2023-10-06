// vfeatDoc.cpp : implementation of the CVfeatDoc class
//

#include <stdafx.h>
#include <afxtempl.h>
#include "vfeat.h"

#include "vfeatDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CVfeatDoc

IMPLEMENT_DYNCREATE(CVfeatDoc, CDocument)

BEGIN_MESSAGE_MAP(CVfeatDoc, CDocument)
	//{{AFX_MSG_MAP(CVfeatDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CVfeatDoc construction/destruction

CVfeatDoc::CVfeatDoc()
{
	// TODO: add one-time construction code here

}

CVfeatDoc::~CVfeatDoc()
{
}

BOOL CVfeatDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}



/////////////////////////////////////////////////////////////////////////////
// CVfeatDoc serialization

void CVfeatDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CVfeatDoc diagnostics

#ifdef _DEBUG
void CVfeatDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CVfeatDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CVfeatDoc commands

void CVfeatDoc::AddLine(CArray <VmPoint3D, VmPoint3D &> & pt_array)
{
	CLine	line;
	line.m_pt_array.Copy(pt_array);
    m_lines.Add(line);
}


void CVfeatDoc::OnDraw(CDC* pDC)
{
	for (int i=0; i< m_lines.GetSize(); i++)
		m_lines[i].OnDraw(pDC);
}
