/* 
 new file for multithreaded version 
 Author: Vsevolod Tanchuk

*/


#pragma once
#include "OptimBase.h"

class CPowellOptim : public COptimBase
{
public:
	CPowellOptim(void *pobj1, double (*nrfunc1)(void *pobj, double [])) : COptimBase(pobj1, nrfunc1){};
	void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret);
};
