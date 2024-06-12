// Compare.h: interface for the CCompare class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COMPARE_H__E8BDEAAE_138E_4408_8622_7E13D8C2E10A__INCLUDED_)
#define AFX_COMPARE_H__E8BDEAAE_138E_4408_8622_7E13D8C2E10A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include"pdb.h"


#define float double
#define MAXATOMS 1000


class CTransform
{
public:
	void Init(double par[]);
	void Move(double xyz[]);
	void SimpleRotate(double xyz[]);
	void Rotate(double xyz[], double der[][3]);
	double sz, cz, sx, cx, sy, cy, dxyz[3];
};

class CAtomPair
{
public:
	double CalcPairMove(CTransform &tra, double dyda[], bool bSave);
	double xyz1[3], xyz2[3], diff;
	char Name[4];
	int res_num;
	int nAltPair;
};


class CCompare  
{
protected:
	double param[7], best_param[7], chisq, alamda, y[MAXATOMS], sig[MAXATOMS];
	double *covar[7], *alpha[7], fcov[49], falf[49];
	int ia[7], x[MAXATOMS];
	double RealOptim2();
	void CenterBoth();
	double SobOptimReal(int nSteps);

public:
	int m_nAtomPairs, m_nAltPairs;
	int m_nResidues;
	CAtomPair Pairs[MAXATOMS], AltPairs[(MAXATOMS / 8)][2];
	int PI_steps;
	bool bUsePair[MAXATOMS];
	bool bCA;
	double ResidueImpact[MAXRES];
	int ResCnt[MAXRES];

	CTransform tra;
	//CTransform dtra[6];
	
	double CalcRMSD(double par[], double der[], bool bFinal = false);
	void PrepareToOptim();
	double Optim();
	double SobOptim();
	int GetAtomPairs(CpdbSite& s1, CpdbSite& s2);
	double CalcMove(int n, double par[], double dyda[], bool bSave = false);
	void AddToGlobalImpact(double GlobResImp[], double MaxImpact[], double MinImpact[], int GlobResCnt[], int MinCase[][2], int MaxCase[][2], int n1, int n2);
	void SaveToFile(int n1, int n2, CpdbSite &s1, CpdbSite &s2);
	void PrepareSymmetric(CpdbSite &s1);
	CCompare();
	virtual ~CCompare();

};

#endif // !defined(AFX_COMPARE_H__E8BDEAAE_138E_4408_8622_7E13D8C2E10A__INCLUDED_)
