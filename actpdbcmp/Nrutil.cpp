/*
 based on nrutil.c from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the same subdirectory or archive
   as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */


#define NR_END 1
#define FREE_ARG char*

#define NRANSI
#include "nrutil.h"

//extern void NRError(char *error_text);

void nrerror(char error_text[])
{
 	//NRError(error_text);
}

/* (C) Copr. 1986-92 Numerical Recipes Software Y"1i41)jD`521R5ks)15. */
