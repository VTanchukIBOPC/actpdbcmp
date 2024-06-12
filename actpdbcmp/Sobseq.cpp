/*
 based on sobseq routine from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

#define NRANSI

#include <string.h>
#include "nrutil.h"
#include "Sobseq.h"

#define MAXBIT 30
#define MAXDIM 6

void CSobseq::sobseq(int *n, double x[])
{
	int j,k,l;
	unsigned long i,im,ipp;

	if (*n < 0) {
		memmove(iv, iv2, sizeof(iv2));
		memset(ix, 0, sizeof(ix));
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
		for (k=1;k<=MAXDIM;k++) {
			for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
			for (j=mdeg[k]+1;j<=MAXBIT;j++) {
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) {
					if (ipp & 1) i ^= iu[j-l][k];
					ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
		fac=1.0/(1L << MAXBIT);
		in=0;
	} else {
		im=in;
		for (j=1;j<=MAXBIT;j++) {
			if (!(im & 1)) break;
			im >>= 1;
		}
		//if (j > MAXBIT) nrerror("MAXBIT too small in sobseq");
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(*n,MAXDIM);k++) {
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
		}
		in++;
	}
}

/* (C) Copr. 1986-92 Numerical Recipes Software Y"1i41)jD`521R5ks)15. */

CSobseq::CSobseq(void)
{
	static unsigned long smdeg[MAXDIM+1]={0,1,2,3,3,4,4};
	static unsigned long sip[MAXDIM+1]={0,0,1,1,2,1,4};
	static unsigned long siv[MAXDIM*MAXBIT+1];

	static unsigned long siv2[MAXDIM*MAXBIT+1]={
		0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

	memmove(mdeg, smdeg, sizeof(smdeg));
	memmove(ip, sip, sizeof(sip));;
	memmove(iv, siv, sizeof(siv));
	memmove(iv2, siv2, sizeof(siv2));   

	fac = 0;
	in = 0;
	memset(ix, 0, sizeof(ix));
	memset(iu, 0, sizeof(iu));
}

#undef MAXBIT
#undef MAXDIM
