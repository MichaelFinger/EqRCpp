#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#define FALSE	0
#define TRUE	1
#define PI  3.1415926535897932384626433832795 

typedef int BOOL;

/* functions related to wrapper and print */
void Wrapper_Smooth_CubSpl(char design, 
                           struct PDATA *xtoy, struct ERAW_RESULTS *r_xtoy, 
						   double *se_xtoy, struct CS_SMOOTH *cs_xtoy,
                           struct PDATA *ytox, struct ERAW_RESULTS *r_ytox, 
						   double *se_ytox, struct CS_SMOOTH *cs_ytox, 
                           double prlow, double prhigh, double s, int rep, 
                           struct PDATA *inall, struct ERAW_RESULTS *r);
void Smooth_CubSpl(char design, struct PDATA *z, struct ERAW_RESULTS *r, 
				   double *se, struct CS_SMOOTH *cs_xtoy, 
				   double prlow, double prhigh, double s, 
			       int rep, int ns1, int ns2, int inverse);
void Print_CubSpl(FILE *fp, char tt[], struct PDATA *xtoy, struct PDATA *ytox,  
                  struct PDATA *inall, struct ERAW_RESULTS *r, int parm);
void Print_Coefs_CubSpl(FILE *fp, struct PDATA *z, char c);

/* functions related to numerical linear algebra */ 
void	dpdch(double * mx,int dim);
void	dtdmt(double *tdm2,double *tdm1, int rdim1, int cdim1);
void	dtdsm(double *tdm2,double c, double *tdm1, int rdim1, int cdim1);
void	dtdmm(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1, int cdim2);
void	dtdma(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1);
void	subbk(double* b, double *utm, int dim);
void	subfd(double* b, double *utm, int dim);
void	chsol(double *vx, double* vb, double *mpd, int dim); 
void	dtdmv(double *v2, double * mtd, double * v1, int rdim, int cdim); 

/* functions related to numerical linear/cubic	*/
/* polynomial			*						*/ 
double	linearPoly(double x0, double y0, double x1, double y1, double xvalue);
double	cubicPoly(double xleft, double xright, double ai, double bi, double ci, 
		double di, double xvalue);

/* functions related to post-smoothing method	*/
/* using cubic splines							*/
void	setupT(BOOL edist, double * mt, int n, double * h);
void	setupQt(BOOL edist, double * mqt,int n, double * h); 
int		sspline(double *x, double *y,double *dyi, int n, double s, double *cmx);  
void	postSmooth(double *x, double *y, double *dyi, int num1, double s, 
		int xlow, int xhigh, double ky,double *vectX, int num2, double *vectY, double *cmat);

/* functions related to the inverse of the      */
/* post-smoothing method using cubic splines	*/
double	inverseCubicPoly(double left, double right, double ai, double bi, 
		double ci, double di, double yvalue);
void	inverseSSpline(double *ynodes, double * cmat, int n2,double *vectX, int n1,double *vectY); 
void	inversePostSmooth(double *yvalues, double *xvalues, double *dxi, int num1, double s, 
		int ylow, int yhigh, double kx,double *vectX, int num2, double *vectY, double *cmat);

#endif 