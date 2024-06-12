// pdb.h: interface for the Cpdb class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PDB_H__AFA10698_D877_4941_821A_83A6AF3B7221__INCLUDED_)
#define AFX_PDB_H__AFA10698_D877_4941_821A_83A6AF3B7221__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>

#define MAXRES 200

class CpdbRec
{
public:
	char Name[4];
	int num;
	char Conform;
};


class CResidue
{
public:
	int GetNum(char **pStr = NULL);
	void ModifyName(int nSh);
	char Name[8], NumStr[6], ChainName[4];
	int m_nNumInSite;
	int num;
};

class CAtom : public CpdbRec
{
public:
	bool IsBackboneAtom();
	CResidue *m_Res;
	double xyz[3];
	void ParseAtomStr(char* str, char conform);
	bool operator>(CAtom &At2);
	bool operator==(CAtom &At2);
};

class CpdbSite  
{
public:

	int m_nResidues;
	CResidue* Residues;
	int m_nAtoms;
	CAtom *Atoms;
	char *FName;
	int m_nSiteNum;
	char SiteName[10];
	char Conform;

	static bool m_bOnlyBackbone;

	static CResidue* FindResidue(CResidue &Res, CResidue ResArr[], int nRes, int nShift, char *ChName);
	static void SortResidues(CResidue ResArr[], int nRes);

	int HasResidue(CResidue &Res);
	int FindItself(CResidue Res[], int nRes, CpdbSite Sites[]);

	void LeaveOnlyBackbone();
	void AllocResidues();
	void SortAtoms();
	int ReadFromFile(char *FName);
	void AllocAtoms();

	void GetSiteName();

	CpdbSite();
	CpdbSite(int nAt, int nRes);
	//CpdbSite(const CpdbSite &st);
	void operator=(CpdbSite &st);
	virtual ~CpdbSite();

};

class CTmpAtom : public CAtom
{
public:
	int m_nSiteNum, m_nAtNum;
};

class CPdbReader
{
protected:
	int ReadResidues(FILE *in, char conform = ' ');
	int m_nTmpAtoms;
	CTmpAtom TmpAtoms[8000], CurAtom;
	CResidue TmpResidues[6000], CurResidue;
	char buff[256];
	bool ParseStr(char conform = '\0');

public:
	int ReadTemplate(char *FName);
	int m_nTmpSites;
	CpdbSite TmpSites[50];
	CpdbSite Template;

	int ReadPDBFile(char *FName);
};


#endif // !defined(AFX_PDB_H__AFA10698_D877_4941_821A_83A6AF3B7221__INCLUDED_)
