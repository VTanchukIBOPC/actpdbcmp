
/*
 parent clas foroptimization routines from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

#pragma once
#define NRANSI
#include "nrutil.h"


class COptimBase
{
protected:
	int ncom;
	double pcom[MAX_VARS], xicom[MAX_VARS];

	double f1dim(double x);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);
	void linmin(double p[], double xi[], int n, double *fret);
	double brent(double ax, double bx, double cx, double tol, double *xmin);

public:
	void *pobj;
	double (*nrfunc)(void *pobj, double []);	

	double golden(double ax, double bx, double cx, double tol, double *xmin);

	COptimBase(void *pobj1, double (*nrfunc1)(void *pobj, double []));
	~COptimBase(void);
};
