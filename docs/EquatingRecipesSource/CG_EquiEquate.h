/* 
  CG_EquiEquate.h  
*/

#include "ERutilities.h"

#ifndef CG_EQUIEQUATE_H
#define CG_EQUIEQUATE_H

void FEorMFE_EE(double w1, int internal, int nsv, 
                int nsx, double minx, double maxx, 
                int nsy, double miny, double maxy, double inc,
                double **bxvin, double **byvin, double rv1, double rv2,
                double *fxs, double *gys, double *eraw,
                double *a, double *b, double *erawBH);
void SyntheticDensities(double w1, int internal, int nsv, int nsx, double **bxv, 
                        int nsy, double **byv, double rv1, double rv2,
                        double *fxs, double *gys);
void MixSmooth(int nsv, int nsx, double unimix, double **braw, 
               double *fx, double *hv);
void CondBivDist(int nsx, int nsv, double **bxv, double *hv);
void BH_LinEq(double minx, double maxx, double miny, double maxy, double inc, 
              double *fxs, double *gys, double *a, double *b);
void ModCondBivDist(int internal, int nsv, int nsx, double rv1, double rv2, 
                    double muv1, double muv2, double **bxv);
void Chained_EE(int nsx, double *prx1, double minv, double maxv,
                double incv, int nsv, double *Gv1, double miny, 
                double incy, int nsy, double *Gy2, double *Gv2,
                double *eraw);
void Print_SynDens(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

#endif 