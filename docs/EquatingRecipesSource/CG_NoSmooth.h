/*
  CG_NoSmooth.h                                            
*/

#include "ERutilities.h"

#ifndef CG_NOSMOOTH_H
#define CG_NOSMOOTH_H

void Wrapper_CN(char design, char method, char smoothing, 
               double w1, int anchor, double rv1, double rv2,
               struct BSTATS *xv, struct BSTATS *yv, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);  
void CI_LinEq(double mnx1, double sdx1, double mnv1, double sdv1, double covxv1,
              double mny2, double sdy2, double mnv2, double sdv2, double covyv2,
              double w1, int anchor, char method, double min, double max,
              double inc, int nm, double *msx ,double *msy, double *ssx, 
              double *ssy, double *gamma1, double *gamma2, double *a, double *b, 
              double **eraw);
void CI_LinObsEq(double mnx1, double sdx1, double mnv1, double sdv1, 
                 double mny2, double sdy2, double mnv2, double sdv2,
                 double w1, char method, double gamma1, double gamma2,
                 double *msx, double *msy, double *ssx, double *ssy,
                 double *a, double *b);
void Print_CN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

#endif 