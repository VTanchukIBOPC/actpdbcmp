/*
 based on nrutil.h routine from 
 Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P. (1995)
 Numerical recipes in C: the art of scientific computing. 
 Cambridge University Press, New York, USA.

 modified by Vsevolod Tanchuk
*/

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#define float double
//#include <math.h>

#define MAX_VARS 10
//added


inline double sqr_nr(float x){ return x * x;}
#define SQR(a) (sqr_nr(a))

#define DSQR(a) (sqr_nr(a))

inline double dmax_nr(double a, double b){ return (a > b) ? a : b;}
#define DMAX(a,b) (dmax_nr(a, b)))

inline double dmin_nr(double a, double b){ return (a < b) ? a : b;}
#define DMIN(a,b) (dmin_nr(a, b))

#define FMAX(a,b) (dmax_nr(a, b))

#define FMIN(a,b) (dmin_nr(a, b))

inline long lmax_nr(long a, long b){ return (a > b) ? a : b;}
#define LMAX(a,b) (lmax_nr(a, b))

inline long lmin_nr(long a, long b){ return (a < b) ? a : b;}
#define LMIN(a,b) (lmin_nr(a, b))

inline long imax_nr(int a, int b){ return (a > b) ? a : b;}
#define IMAX(a,b) (lmax_nr(a, b))

inline long lmin_nr(int a, int b){ return (a < b) ? a : b;}
#define IMIN(a,b) (lmin_nr(a, b))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#ifdef __cplusplus
extern "C" {
#endif

void nrerror(char error_text[]);

#ifdef __cplusplus
}
#endif

#endif /* _NR_UTILS_H_ */
