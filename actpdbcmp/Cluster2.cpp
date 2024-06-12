//version 17/09/95
//#include <stdio.h>
#include <math.h>
#include "claster.h"


/*   IMSL ROUTINE NAME   - OCLINK                                                
 *
 *-----------------------------------------------------------------------        
 *
 *   PURPOSE             - PERFORM A SINGLE-LINKAGE OR COMPLETE-LINKAGE          
 *                           HIERARCHICAL CLUSTER ANALYSIS GIVEN A               
 *                           SIMILARITY MATRIX                                   
 *
 *   USAGE               - CALL OCLINK (ND,IOPT,XSIM,CLEVEL,ICLSON,              
 *                                     ICRSON,IPTR,IER)                          
 *
 *   ARGUMENTS    ND     - INPUT NUMBER OF DATA POINTS TO BE CLUSTERED.          
 *                           ND-1 CLUSTERS ARE FORMED NUMBERED CONSEC-           
 *                           UTIVELY FROM ND+1 TO ND+(ND-1).  ND MUST            
 *                           BE GREATER THAN 2. (SEE REMARKS)                    
 *                IOPT   - INPUT OPTIONS VECTOR OF LENGTH 2.                     
 *                           IOPT(1)=0 IMPLIES SINGLE-LINKAGE DESIRED.           
 *                             OTHERWISE, COMPLETE-LINKAGE IS PERFORMED.         
 *                           IOPT(2)=0 IMPLIES THE SIMILARITIES ARE              
 *                             DISTANCE-LIKE (I.E.  SMALLER IMPLIES              
 *                             CLOSER).  OTHERWISE, THE SIMILARITIES ARE         
 *                             ASSUMED TO BE POSITIVE CORRELATION-LIKE           
 *                             SIMILARITIES. (SEE REMARKS).                      
 *                XSIM   - INPUT/OUTPUT VECTOR OF LENGTH ((ND+1)*ND)/2           
 *                           CONTAINING THE SIMILARITY MATRIX IN SYM-            
 *                           METRIC STORAGE MODE.  FOR I GREATER THAN J,         
 *                           XSIM(((I-1)*I)/2 + J) CONTAINS THE SIMILAR-         
 *                           ITY OF THE I-TH AND J-TH DATA POINTS.               
 *                           ON INPUT, THE DIAGONAL ELEMENTS (SIMILARITY         
 *                           TO ITSELF) ARE ARBITRARY AND DO NOT NEED TO         
 *                           BE DEFINED.  XSIM IS DESTROYED ON OUTPUT.           
 *                CLEVEL - OUTPUT VECTOR OF LENGTH ND-1.  CLEVEL(K)              
 *                           CONTAINS THE SIMILARITY LEVEL AT WHICH              
 *                           CLUSTER ND+K WAS FORMED.                            
 *                ICLSON - OUTPUT VECTOR OF LEFTSONS OF LENGTH ND-1.             
 *                           CLUSTER NUMBER ND+K WAS FORMED BY MERGING           
 *                           CLUSTER ICLSON(K) WITH CLUSTER ICRSON(K).           
 *                ICRSON - OUTPUT VECTOR OF RIGHTSONS OF LENGTH ND-1.            
 *                           THE RIGHTSON OF CLUSTER ND+K IS CONTAINED           
 *                           IN ICRSON(K).                                       
 *                IPTR   - WORK VECTOR OF LENGTH ND.                             
 *
 *   PRECISION/HARDWARE  - SINGLE/ALL                                            
 *
 *   REQD. IMSL ROUTINES - UERTST,UGETIO                                         
 *
 *   NOTATION            - INFORMATION ON SPECIAL NOTATION AND                   
 *                           CONVENTIONS IS AVAILABLE IN THE MANUAL              
 *                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP          
 *
 *   REMARKS  1.  THE DATA CLUSTERS ARE NUMBERED 1 TO ND. THE ND-1               
 *                CLUSTERS FORMED BY MERGING ARE NUMBERED ND+1 TO                
 *                ND+(ND-1) AND DECREASE IN SIMILARITY, MAKING IT EASY           
 *                TO IDENTIFY THE MOST SIMILAR CLUSTERS.                         
 *            2.  SIMILARITIES GENERALLY SHOULD BE NONNEGATIVE. RAW              
 *                CORRELATIONS TAKE ON VALUES R WHERE R LIES IN THE              
 *                CLOSED INTERVAL (-1, 1) AND SHOULD BE MADE POSITIVE.           
 *                IF R = 1 AND R = -1 BOTH MEAN HIGH SIMILARITY, THEN THE        
 *                TRANSFORMATIONS, RR EQUAL TO THE ABSOLUTE VALUE OF R,          
 *                OR, RR EQUAL TO R SQUARED, ARE APPROPRIATE. IF R = -1          
 *                REPRESENTS VERY LOW SIMILARITY, THEN THE TRANSFORMA-           
 *                TION, RR = 1-R, BECOMES A DISTANCE-LIKE SIMILARITY.            
 *            3.  NOTE THAT THE ORIGINAL DATA MATRIX OF THE USER DOES            
 *                NOT ENTER OCLINK, ALLOWING THE USER TO DEFINE, VIA             
 *                XSIM, WHATEVER MEASURE OF SIMILARITY SEEMS MOST                
 *                APPROPRIATE.                                                   
 *            4.  A COMPUTER PRINTING OF A BINARY TREE IS AVAILABLE              
 *                THROUGH IMSL ROUTINE USTREE.                                   
 *
 *   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.               
 *
 *   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN         
 *                           APPLIED TO THIS CODE. NO OTHER WARRANTY,            
 *                           EXPRESSED OR IMPLIED, IS APPLICABLE.                
 *
 *-----------------------------------------------------------------------        
 */
void clust(int nd,int*iopt,double*xsim,double*clevel,int*iclson,int*icrson,int*iptr)
{
int  i, icl, ilo, im1, inum, iopt1, ip, ir, irr, is, isgn, ism1, 
	 isp1, istar, itype, iup, j, jp, jsave, jsm1, jsp1, jstar, 
	 jtype, kstar;
double temp, xmin;
static double zero = 0.0, tenth = 0.1, onep1 = 1.1, xbig = 1.e9;

//if(nd<3)error("Error in Cluster, number of point =%i<3",nd);

iopt1 = iopt[0];
isgn=(iopt[1])?-1:1;

xsim[0] = onep1;
ir = 1;
im1 = 1;
for( i = 2; i <= nd; i++ ){
		if( isgn == 1 )
			goto L_10;
		for( j = 0; j < im1; j++ )
			xsim[ir + j] = -xsim[ir + j];
L_10:
		ir = ir + i;
		im1 = i;
		xsim[ir - 1] = (double)( i ) + tenth;
		}
	/*                                  IPTR(I) CONTAINS THE COLUMN WITH             
	 *                                    THE SMALLEST SIMILARITY IN ROW I           
	 *                                    (IPTR(1) UNDEFINED-ROW 1 IS NULL)           */
	iptr[1] = 1;
	/*                                  INUM POINTS TO ROW WITH SMALLEST             
	 *                                    SIMILARITY.                                 */
	inum = 2;
	xmin = xsim[1];
	ir = 3;
	im1 = 2;
	for( i = 3; i <= nd; i++ ){
		temp = xsim[ir];
		jsave = 1;
		for( j = 1; j < im1; j++ ){
			if( temp < xsim[ir + j] )
				goto L_20;
			temp = xsim[ir + j];
			jsave = j + 1;
L_20:
			;
			}
		im1 = i;
		ir = ir + i;
		iptr[i - 1] = jsave;
		if( xmin < temp )
			goto L_25;
		inum = i;
		xmin = temp;
L_25:
		;
		}
	icl = 0;
	/*                                  MERGE THE TWO CLUSTERS CLOSEST               
	 *                                    TOGETHER.  THEY ARE IN ROWS                
	 *                                    ISTAR AND JSTAR.                            */
L_30:
	icl = icl + 1;
	istar = inum;
	jstar = iptr[inum - 1];
	ism1 = ((istar - 1)*istar)/2;
	jsm1 = ((jstar - 1)*jstar)/2;
	isp1 = ((istar + 1)*istar)/2;
	jsp1 = ((jstar + 1)*jstar)/2;
	/*                                  UPDATE LEFT AND RIGHT SON ARRAYS             
	 *                                  RECOVER CURRENT CLUSTER NUMBER               
	 *                                    ON DIAGONAL.                                */
	itype = int(xsim[isp1 - 1]);
	jtype = int(xsim[jsp1 - 1]);
	if( itype > nd && jtype > nd )
		goto L_35;
	/*                                  AT LEAST ONE SON IS A PURE DATA              
	 *                                    CLUSTER.  MAKE RIGHTSON DATA ALWAYS         */
	iclson[icl - 1] = (itype>jtype)?itype:jtype;
	icrson[icl - 1] = (itype<jtype)?itype:jtype;
	goto L_45;
	/*                                  BOTH SONS ARE HIGH LEVEL CLUSTERS.           
	 *                                  MAKE LEFTSON MORE SIMILAR OF TWO              */
L_35:
	if( clevel[jtype - nd - 1] < clevel[itype - nd - 1] )
		goto L_40;
	iclson[icl - 1] = itype;
	icrson[icl - 1] = jtype;
	goto L_45;
L_40:
	iclson[icl - 1] = jtype;
	icrson[icl - 1] = itype;
	/*                                  RECORD LEVEL OF MERGING                       */
L_45:
	clevel[icl - 1] = xsim[ism1 + jstar - 1];
	if( icl >= nd - 1 )
		goto L_125;
	/*                                  MARK MERGED CLUSTERS-DELETE ROW ISTAR         */
	xsim[jsp1 - 1] = (double)( nd + icl ) + tenth;
	xsim[isp1 - 1] = -tenth;
	/*                                  UPDATE SIMILARITY MATRIX                      */
	ir = 0;
	for( i = 1; i <= nd; i++ ){
		irr = ir;
		ir = ir + i;
		/*                                  SKIP CURRENT AND DELETED CLUSTERS             */
		if( xsim[ir - 1] < zero )
			goto L_75;
		if( i == jstar )
			goto L_75;
		/*                                  FIND XSIM(ISTAR,I) AND XSIM(JSTAR,I)          */
		if( i > istar )
			goto L_50;
		ip = ism1 + i;
		goto L_55;
L_50:
		ip = irr + istar;
		goto L_60;
L_55:
		if( i > jstar )
			goto L_60;
		jp = jsm1 + i;
		goto L_65;
L_60:
		jp = irr + jstar;
		/*                                  SINGLE/COMPLETE LINKAGE OPTION                */
L_65:
		if( iopt1 == 0 )
			goto L_70;
		xsim[jp - 1] =(xsim[ip-1]>xsim[jp-1])?xsim[ip-1]:xsim[jp-1];
		goto L_75;
L_70:
		xsim[jp - 1] =(xsim[ip-1]<xsim[jp-1])?xsim[ip-1]:xsim[jp-1];
L_75:
		;
		}
	/*                                  UPDATE IPTR ARRAY, FINDING NEW ROW           
	 *                                    MINIMUMS AND OVERALL MINIMUM ALSO           */
	xmin = xbig;
	/*                                  FIRST DO THE NEW CLUSTER ROW, JSTAR.         
	 *                                  IF JSTAR IS LESS THAN 3, THEN                
	 *                                    IPTR ARRAY IS OK.                           */
	if( jstar < 3 )
		goto L_85;
	temp = xbig;
	iup = jstar - 1;
	ir = 0;
	for( j = 1; j <= iup; j++ ){
		ir = ir + j;
		if( xsim[ir - 1] < zero )
			goto L_80;
		jp = jsm1 + j;
		if( temp < xsim[jp - 1] )
			goto L_80;
		temp = xsim[jp - 1];
		iptr[jstar - 1] = j;
L_80:
		;
		}
	/*                                  FIND OVERALL MIN OF 1ST JSTAR ROWS            */
L_85:
	if( jstar == 1 )
		goto L_95;
	ir = 1;
	for( i = 2; i <= jstar; i++ ){
		irr = ir;
		ir = ir + i;
		if( xsim[ir - 1] < zero )
			goto L_90;
		ip = irr + iptr[i - 1];
		if( xmin < xsim[ip - 1] )
			goto L_90;
		xmin = xsim[ip - 1];
		inum = i;
L_90:
		;
		}
	/*                                  UPDATE SIMILARITIES IN ROWS AFTER            
	 *                                    JSTAR.                                      */
L_95:
	ilo = jstar + 1;
	ir = ((ilo - 1)*ilo)/2;
	/*                                  SINCE JSTAR .LT. ISTAR, ILO .LE. ND          
	 *                                    IS ASSURED                                  */
	im1 = ilo - 1;
	for( i = ilo; i <= nd; i++ ){
		is = ir;
		ir = ir + i;
		if( xsim[ir - 1] < zero )
			goto L_120;
		kstar = iptr[i - 1];
		/*                                  IF KSTAR=ISTAR, THE LAST MINIMUM             
		 *                                    WAS DELETED-SEARCH WHOLE ROW                */
		if( kstar == istar )
			goto L_105;
		if( iopt1 == 0 )
			goto L_100;
		/*                                  COMPLETE LINKAGE OPTION.  VALUE IN           
		 *                                    COLUMN JSTAR WAS INCREASED.  IF            
		 *                                    KSTAR=JSTAR, SEARCH WHOLE ROW.             
		 *                                    OTHERWISE IPTR AND MIN ARE OK.              */
		if( kstar == jstar )
			goto L_105;
		goto L_115;
		/*                                  SINGLE LINKAGE OPTION.  VALUE IN             
		 *                                    COLUMN JSTAR WAS DECREASED.  IF            
		 *                                    KDTAR=JSTAR, IPTR IS OK. OTHERWISE,        
		 *                                    COMPARE OLD MIN KSTAR WITH JSTAR.           */
L_100:
		if( kstar == jstar )
			goto L_115;
		ip = is + kstar;
		jp = is + jstar;
		if( xsim[ip - 1] < xsim[jp - 1] )
			goto L_115;
		iptr[i - 1] = jstar;
		goto L_115;
		/* SEARCH ENTIRE ROW FOR MINIMUM */
L_105:
		temp = xbig;
		irr = 0;
		for( j = 1; j <= im1; j++ ){
			irr = irr + j;
			if( xsim[irr - 1] < zero )
				goto L_110;
			jp = is + j;
			if( temp < xsim[jp - 1] )
				goto L_110;
			temp = xsim[jp - 1];
			iptr[i - 1] = j;
L_110:
			;
			}
		/*   LOOK FOR OVERALL MINIMUM */
L_115:
		ip = is + iptr[i - 1];
		if( xmin < xsim[ip - 1] )
			goto L_120;
		xmin = xsim[ip - 1];
		inum = i;
L_120:
		im1 = i;
		}
	goto L_30;
	/* CHECK CLEVEL USING IOPT AND ISGN */
L_125:
	if( isgn == 1 ) return;
	for( i = 0; i < icl; i++ )
		clevel[i] = -clevel[i];
}

