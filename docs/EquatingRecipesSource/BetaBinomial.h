/* 	
  BetaBinomial.h 
*/ 

#include "ERutilities.h"

#ifndef BETABINOMIAL_H
#define BETABINOMIAL_H

#define RMAXIT 20	  /* maximum number of iterations for computing upper limit */
#define KACC .001 	    /* accuracy for computing lower limit in CalcBetaParaLS */

void Wrapper_Smooth_BB(struct USTATS *x, int nparm, double rel,
                      struct BB_SMOOTH *s);
void Smooth_BB(int n, int nitems, double *fd, double *mts, 
               int nparm, double rel, struct BB_SMOOTH *s);
void Print_BB(FILE *fp, char tt[], struct USTATS *x, struct BB_SMOOTH *s);
void Wrapper_RB(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct BB_SMOOTH *bbx, struct BB_SMOOTH *bby, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);
void Print_RB(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

short Beta4Smooth(double *rfreq, struct BB_SMOOTH *betafit);
short Beta2Smooth(double *rfreq, struct BB_SMOOTH *betafit);

short CalcBetaPara(int nitems,double *moment,double *nctmoment,double *para);
short EstNegHypGeo(int nitems,double *moment,double *para);
short CalcBetaParaMM(int nitems,double *moment,double *para);
short CalcBetaParaLS(int nitems,double *moment,double *nctmoment,double *para);
short FindUpper(double *para,double *tmoment);
short FindLower(double *para,double *tmoment);
double CalcKurt(double *para);
short Kurtfuncd(double x,int nitems,double *para,double *tmoment,double *nctmoment,
		double *kurt);
short KurtfuncUpper(double x,int nitems,double *para,double *tmoment,double *nctmoment,
		double *kurt);
		
void BetaMoments(int n,double k, const double *rmoment,double *tmoment,double *nctmoment);

short ObsDensity(int n, int nsamp, double *beta, double *scounts);
void ObsDenK(int n,double k,double *scounts);
double dotprod(register int n, register double *v1, register double *v2);
void CalcM3(int n, double *beta, register double **m3);
short CalcM24(int n, double para, double *m);
void CalcM15(int n, double para, double *m);

double LRChiSqr(int n, register const double *rawc, register const double *fitc,
	register int *ncat);
double PChiSqr(int n, double minexpcount, register const double *rawc,
	register const double *fitc,register int *ncat);
	
double CalcLordk(double kr20, int nitems, double *rmoment);

#endif 
