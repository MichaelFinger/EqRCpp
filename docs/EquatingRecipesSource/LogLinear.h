/*
  LogLinear.h
*/

#include "ERutilities.h"

#ifndef LOGLINEAR_H
#define LOGLINEAR_H

void mtranspose(double **a, int nr, int nc, double **at);
void mmult_a_b(double **a, double **b, 
               int nra, int nca, int nrb, int ncb,
               double **s);
void mmult_a_v(double **a, double *v, int nr, int nc, double *s);
void design_matrix(int nsu, double minu, double incu, 
                   int nsv, double minv, double incv,  
                   int cu, int cv, int cuv, int cpm[][2], int scale, 
                   double **B_raw, double **B);
void get_nct_bfd(int anchor, int nsx, int nsv, double **bfd, 
                 double *nct);
void get_bfd_mct(int anchor, int nsx, int nsv, double *mct, 
                 double **bfd);
void get_BtSmB(double **B, double *m, int ns, int nc, double N,
               double **BtSmB);
void get_Btnm(double **B, double *n, double *m, int ns, int nc, 
              double *Btnm);
void get_Beta0(double **B, double *n, double N, int ns, int nc,
               double *Beta0, FILE *fp);
double get_mct(double **B, double *Beta, double *uin, double N, 
               int ns, int nc, double *m, FILE *fp);
int iteration(FILE *fp, double **B, double **B_raw, double *nct, double N, 
              double *uin, int ns, int cu, int cv, int cuv, int cpm[][2], 
              int max_nit, int ctype, int Btype, double crit,
              double *Beta, double *mct, double *n_mts, double *m_mts, 
              double *n_mts_raw, double *m_mts_raw, double *lrc, 
			  int *nzero, double *ap);
void get_LLmoments(double **B, double **B_raw, double *f, double N,  
                   int ns, int cu, int cv, int cuv, int cpm[][2], 
                   double *mts, double *mts_raw);
int crit_mts(int nc, int cu, int ctype, int Btype,
             double *n_mts, double *m_mts, double crit);
void Print_iteration_heading(FILE *fp, int ns, int nc, double *nct,
                             double *n_mts, double *n_mts_raw,
                             int ctype, int Btype, double crit);
void Wrapper_Smooth_ULL(struct USTATS *x, int c,
						int scale, int Btype, int ctype, double crit,
						FILE *fp, struct ULL_SMOOTH *s);
void Smooth_ULL(int n, int ns, double min, double inc,
                double *fd, int c, 
				int scale, int Btype, int ctype, double crit,
				FILE *fp, struct ULL_SMOOTH *s);
void Print_ULL(FILE *fp, char tt[], struct USTATS *x,
               struct ULL_SMOOTH *s, int print_dm, int print_mts);
void Wrapper_Smooth_BLL(struct BSTATS *xv, int anchor,
                        int cu, int cv, int cuv, int cpm[][2],
						int scale, int Btype, int ctype, double crit,
                        FILE *fp, struct BLL_SMOOTH *s);
void Smooth_BLL(int n, int nsu, double minu, double incu,
                int nsv, double minv, double incv, double *nct, 
                int anchor, int cu, int cv, int cuv, int cpm[][2],
				int scale, int Btype, int ctype, double crit,
                FILE *fp, struct BLL_SMOOTH *s);
void Print_BLL(FILE *fp, char tt[], struct BSTATS *xv, struct BLL_SMOOTH *s, 
			   int print_dm, int print_mts, int print_freq, int print_bfd);
void Wrapper_RL(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);
void Print_RL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
void Wrapper_SL(char design, char method, char smoothing,  struct BSTATS *xy, 
               struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall, 
			   struct ERAW_RESULTS *r);
void Print_SL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
void Wrapper_CL(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r);
void Print_CL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

#endif 