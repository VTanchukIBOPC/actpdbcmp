
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the CLASTER_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// CLASTER_API functions as being imported from a DLL, wheras this DLL sees symbols
// defined with this macro as being exported.
/*
#ifdef CLASTER_EXPORTS
#define CLASTER_API __declspec(dllexport)
#else
#define CLASTER_API __declspec(dllimport)
#endif
*/
// This class is exported from the claster.dll


extern void clust(int nd,int*iopt,double*xsim,double*clevel,int*iclson,int*icrson,int*iptr);
//extern "C" CLASTER_API void clast2(int nd,int*iopt,double*xsim,double*clevel,int*iclson,int*icrson,int*iptr);

