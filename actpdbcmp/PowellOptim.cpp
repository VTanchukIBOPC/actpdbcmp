
/*
 based on powell routine from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

#include <math.h>
#include "PowellOptim.h"
#define NRANSI
#include "nrutil.h"
#define ITMAX 500


void CPowellOptim::powell(double p[], double **xi, int n, double ftol, int *iter, double *fret)
{
	int i,ibig,j;
	//double del,fp,fptt,t,*pt,*ptt,*xit;
	double del,fp,fptt,t,pt[MAX_VARS],ptt[MAX_VARS],xit[MAX_VARS];

	*fret=(*nrfunc)(pobj, p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			return;
		}
		if (*iter == ITMAX){
			nrerror("powell exceeding maximum iterations.");
			return;
			}
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*nrfunc)(pobj, ptt);
		if(fptt > 0 && *iter > n){
			nrerror("wrong direction"); 
			*iter = -10; //to sign direction failure
		}
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
#undef NRANSI
