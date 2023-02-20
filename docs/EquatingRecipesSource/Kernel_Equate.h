/*
  Kernel_Equate.h
*/

#ifndef KERNEL_EQUATE_H
#define KERNEL_EQUATE_H

typedef double (*PTR_FTN1)(int, double*, double*, double, double );

double KernelContinuPdf(int ncat, double *scores, double *fd, double hx, double x);
double KernelContinuCdf(int ncat, double *scores, double *fd, double hx, double x);
double StdNormalCdf(double x) ;
double StdNormalPdf(double x);
double NormalPdf(int npara, double *para, double x);
void CalcKernelMoments(PTR_FTN1 pdf, double a, double b, int nDistCat, double *score, 
					   double *fd, double hx, double *moments);

double Pen1(int nDistCat, double *scores, double *fd, double hx);

double KernelPdfDerivative(int nDistCat, double *scores, double *fd, double hx, double x);

double Pen2(int nDistCat, double *scores, double *fd, double hx);

double Pen(int nDistCat, double *scores, double *fd, double hx, double K);

double Optimalh(int nDistCat, double *scores, double *fd, double K);

double KernelInverseCdf(int nDistCat, double *scores, double *fd, double h, double cdf);

void KernelEquate(int nDistCatx, double *scoresx, double *fdx, double hx, 
					int nDistCaty, double *scoresy, double *fdy, double hy,
					double *Equatedx);

void ComputeCmatrix(int ncat, int degree, long np, double **B, double *fd, double **Cr);

void ComputeCmatrixGen(int ncat, int degree, long np, double **B, double *fd, double **Cr);

void PartialFPartialr(int nDistCat, double *scores, double *fd, double hx, double *Fr, double x);

double FrCrSqNorm(int nDistCat, int degree, double *Fr, double **Cr);

void vPMN(int ncatx, int ncaty, double **bdist, double *vP, double **M, double **N);

void vPT(int ncatx, int ncaty, double **bdist, double *vPP);

void MatrixMultiVector(int nrow, int ncol, double **m, double *v, double *r);

void MatrixMultiMatrix(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r);

void VectorMultiMatrix(int nrow, int ncol, double *v, double **m, double *r);

double VectorNormSq(int ncat, double *v);

void KernelEquateSEERG(int nDistCatx, int degreex, long npx, double *scoresx, double *fdx, 
					   int nDistCaty, int degreey, long npy, double *scoresy, 
					   double *fdy, double *Equatedx, double *SEE);

void KernelEquateSEESG(struct BLL_SMOOTH *bivar, double *Equatedx, double *SEE);

void KernelEquateSG(struct BLL_SMOOTH *bivar, double *Equatedx);

void KernelEquateSEECB(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *wts, double *Equatedx, 
					   double *SEE);

void KernelEquateSEECB2(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *wts, double *Equatedx, 
					   double *SEE);

void KernelEquateNEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx);

void KernelBootStrapNEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						    double *BootStrapSEE, double *MeanEq);

void KernelBootStrapNEATPS2(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						    double *BootStrapSEE, double *MeanEq);

void KernelEquateSEENEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx, double *SEE);

void KernelEquateSEENEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2,  
							double *Equatedx, double *SEE);

void KernelEquateNEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *Equatedx);

void KernelBootStrapNEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						    double *BootStrapSEE, double *MeanEq);

void KernelBootStrapNEATChn2(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						    double *BootStrapSEE, double *MeanEq);

void Wrapper_RK(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);

void Print_RK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

void Wrapper_SK(char design, char method, char smoothing,  struct BSTATS *xy, 
               struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall, 
			   struct ERAW_RESULTS *r);

void Print_SK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

void Wrapper_CK(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r);           

void Print_CK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

double var(double *data, unsigned long n);

double mean(double *data, unsigned long n);

#endif 

