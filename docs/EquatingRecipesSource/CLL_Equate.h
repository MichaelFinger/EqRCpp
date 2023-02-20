/*
  CLL_Equate.h                                            
*/

#ifndef CLL_EQUATE_H
#define CLL_EQUATE_H

typedef double (*PTR_FTN2)(int, double*, double );
typedef double (*PTR_FTN3)(int, double*, double, int );
typedef double (*PTR_FTN4)(double, double, int, double *, double);
typedef double (*PTR_FTN5)(struct BLL_SMOOTH *, double, double);


double CLLEGPdf(double min, double max, int npara, double *para,  
					double x, double nc);

double CLLEGCdf(double min, double max, int npara, double *para,  
					double x, double nc);

double GaussianQuadrature16(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature32(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature64(PTR_FTN2 func, double a, double b, int npara,
						  double *para);

double GaussianQuadrature64i(PTR_FTN3 func, double a, double b, int npara,
						  double *para,  int i);

double ExpPolynomial(int npara, double *para, double x);

double ExpPolynomialxi(int npara, double *para, double x, int i);

void CalcLLContinuMoments(PTR_FTN4 pdf, double a, double b, int npara,
						  double *para, double *moments);
void CLLEquateEG(double minx, double maxx, int nparax, double *parax,  
			   double miny, double maxy, int nparay, double *paray,	
			   int nDistCatx, double *scoresx, double *Equatedx);

double CLLInverseCdf(double min, double max, int npara, double *para, 
							  double cdf, double nc);

double BivExpPolynomial(struct BLL_SMOOTH *bivar, double x, double y);

double BivGaussianQuadrature64(PTR_FTN5 func, struct BLL_SMOOTH *bivar, double ax, 
							   double bx, double ay, double by);

double BivGaussianQuadrature32(PTR_FTN5 func, struct BLL_SMOOTH *bivar, double ax, double bx,
							   double ay, double by);

double CLLMargYPdf(struct BLL_SMOOTH *bivar, double y, double nc);

double CLLMargXPdf(struct BLL_SMOOTH *bivar, double x, double nc);

		/* nc was the last parameter in CLLBivPdf(), but not in function;
		   nc eliminated by rlb on 8-18-09=8 */

double CLLBivPdf(struct BLL_SMOOTH *bivar, double x, double y);

double CLLBivCdf(struct BLL_SMOOTH *bivar, double x, double y, double nc);

double CLLNEATPSMargPdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *fv, 
						double wts, double x);

double CLLNEATPSMargCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
								 double x);

double CLLNEATPSInverseCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, 
										double wts, double cdf);

int CLLEquateNEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx);

double CLLMargInverseXCdf(struct BLL_SMOOTH *bivar1, double xcdf, double nc);

double CLLMargInverseYCdf(struct BLL_SMOOTH *bivar1, double ycdf, double nc);

double CLLMargInverseCBYCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wtsy, double ycdf, double nc1, double nc2);

int CLLEquateNEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *Equatedx);

int CLLEquateSG(struct BLL_SMOOTH *bivar, double *Equatedx);

void CLLSEEEG(long npx, int nparax, double *parax, int nCatx, double *scoresx,
			  long npy, int nparay, double *paray,	int nCaty, double *scoresy, double *SEE);

void Wrapper_RC(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r);

void Print_RC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

void Wrapper_SC(char design, char method, char smoothing,  struct BSTATS *xy, 
               struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall, 
			   struct ERAW_RESULTS *r);

void Print_SC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

void Wrapper_CC(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r);            

void Print_CC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);

#endif 
