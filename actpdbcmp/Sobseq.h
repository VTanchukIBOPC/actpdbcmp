/*
 based on sobseq routine from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

#pragma once

#define MAXBIT 30
#define MAXDIM 6

class CSobseq
{
protected:
	double fac;
	unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
	unsigned long mdeg[MAXDIM+1];
	unsigned long ip[MAXDIM+1];
	unsigned long iv[MAXDIM*MAXBIT+1];

	unsigned long iv2[MAXDIM*MAXBIT+1];
public:
	void sobseq(int *n, double x[]);

	CSobseq(void);
};

#undef MAXBIT
