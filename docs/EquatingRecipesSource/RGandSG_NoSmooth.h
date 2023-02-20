/*
  RGandSG_NoSmooth.h                                        
*/

#include "ERutilities.h"

#ifndef RGANDSG_NOSMOOTH_H
#define RGANDSG_NOSMOOTH_H

void Wrapper_RN(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);
void Wrapper_SN(char design, char method, char smoothing,  struct BSTATS *xy, 
               int rep, struct PDATA *inall, struct ERAW_RESULTS *r);
void RGandSG_LinEq(double mnx, double sdx, double mny, double sdy,
                   char method, double min, double max, double inc, 
                   double *a, double *b, double *eraw);
void Print_RN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
void Print_SN(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

#endif 