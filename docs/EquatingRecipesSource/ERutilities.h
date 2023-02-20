 /*
ERutilities.h 

  Includes:  (a) all structure declarations, and 
             (b) prototypes for functions in ERutilities.c                                        
*/


#ifndef ERUTILITIES_H
#define ERUTILITIES_H

#pragma warning(disable:4996)
#pragma warning(disable:4706)
#pragma warning(disable:4715)
#pragma warning(disable:4701)
#pragma warning(disable:4127)
#pragma warning(disable:4244)
#pragma warning(disable:4100)
#pragma warning(disable:4305) 

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "IRTst.h"
#include "IRTeq.h"
/************************ structure definitions **************************/

struct USTATS{
  /* 
    raw-score statistics for a univariate distribution 
  */
  char fname[100];                                 /* name of input file */
  char id;                                        /* single character id */
  int n;                                          /* number of examinees */
  double mind;                                      /* min score in data */
  double maxd;                                      /* max score in data */
  double min;                                      /* min score for fd[] */
  double max;                                      /* max score for fd[] */
  double inc;                       /* increment between adjacent scores */
  int ns;                            /* number of scores (or categories) */
  int *fd;                                 /* freq dist fd[0]...fd[ns-1] */
  double *dbl_fd;                              /* double version of fd[] */
  int *cfd;                                             /* cum freq dist */
  double *rfd;                                     /* relative freq dist */
  double *crfd;                                /* cum relative freq dist */
  double *prd;                                   /* percentile rank dist */
  double mts[4];                        /* moments: mean, sd, skew, kurt */
};

struct BSTATS{
  /* 
    raw-score statistics for a bivariate distribution 
  */
  char fname[100];                                 /* name of input file */
  int n;                                          /* number of examinees */
                 /* var 1 -- rows of bfd[][] */
  char id1;                                /* Var 1: single character id */
  double mind1;                              /* var 1: min score in data */
  double maxd1;                              /* var 1: max score in data */
  double min1;                             /* var 1: min score for fd1[] */
  double max1;                             /* var 1: max score for fd1[] */
  double inc1;               /* var 1: increment between adjacent scores */
  int ns1;                                    /* Var 1: number of scores */
  int *fd1;                      /* var 1: freq dist fd1[0]...fd1[ns1-1] */
  double *dbl_fd1;                            /* double version of fd1[] */
  int *cfd1;                                     /* Var 1: cum freq dist */
  double *rfd1;                             /* Var 1: relative freq dist */
  double *crfd1;                        /* Var 1: cum relative freq dist */
  double *prd1;                           /* Var 1: percentile rank dist */
  double mts1[4];                /* var 1: moments: mean, sd, skew, kurt */
                 /* var 2 -- columns of bfd[][] */ 
  char id2;                                /* Var 2: single character id */
  double mind2;                              /* var 2: min score in data */
  double maxd2;                              /* var 2: max score in data */
  double min2;                             /* var 2: min score for fd2[] */
  double max2;                             /* var 2: max score for fd2[] */
  double inc2;               /* var 2: increment between adjacent scores */
  int ns2;                                    /* Var 2: number of scores */
  int *fd2;                      /* var 2: freq dist fd2[0]...fd2[ns2-1] */
  double *dbl_fd2;                            /* double version of fd2[] */
  int *cfd2;                                     /* Var 2: cum freq dist */
  double *rfd2;                             /* Var 2: relative freq dist */
  double *crfd2;                        /* Var 2: cum relative freq dist */
  double *prd2;                           /* Var 2: percentile rank dist */
  double mts2[4];                /* var 2: moments: mean, sd, skew, kurt */
  
  int **bfd;        /* bivariate freq dist bfd[0][0]...bfd[ns1-1][ns2-1] */
  double **dbl_bfd;                         /* double version of bfd[][] */
  double cov;                                              /* covariance */
  double corr;                                            /* correlation */             
  double **bp12;                          /* biv proportions[var1][var2] */                                                                 
};

struct PDATA{
  /* 
    structure that contains all input for a particular 
    design/method/smoothing.  Used extensively in the argument lists for
    for Wrapper and Print functions 
  */
  char xfname[100];                        /* name of x or xv input file */
  char yfname[100];                        /* name of y or yv input file */
  char xyfname[100];                            /* name of xy input file */
  struct USTATS *x;        /* structure for summary raw data for x and v */
  struct USTATS *y;        /* structure for summary raw data for y and v */
  struct BSTATS *xv;       /* structure for summary raw data for x and v */
  struct BSTATS *yv;       /* structure for summary raw data for y and v */
  struct BSTATS *xy;       /* structure for summary raw data for x and y */
  char design;                       /* 'R' = RG; 'S' = SG; 'C' = CINEG  */
  char method; /* 'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs); 
			      'For CG design,
			      'E' = Freq esti FE) with Braun-Holland (BH-FE) 
                  'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE) 
                  'G' = FE + BH-FE + MFE + BH-MFE
                  'C' = Chained
			      'H' = FE + BH-FE + Chained
                  'A' = FE + BH-FE + MFE + BH-MFE + Chained */
  char smoothing; /* 'N'= no; 'L' = loglin; 'B' = betabin; 'S' = cub spl;
				                                 'K' = kernel; 'Z' = CLL */
  double w1;                               /* weight for synthetic pop 1 */
  int anchor;                          /* = 0 (external); = 1 (internal) */
  double rv1;                    /* reliability of common items in pop 1 */
  double rv2;                    /* reliability of common items in pop 2 */
  int nm;                                           /* number of methods */
  char **names;                                      /* names of methods */
  double min;                                         /* min score for x */
  double max;                                         /* max score for x */
  double inc;                 /* increment between adjacent scores for x */
  int *fdx;                                         /* fd for new form x */
  int n;                           /* number of examinees for new form x */
  double minp;      /* min raw score for yct-- see ReadSSConvTableForY() */
  double maxp;      /* max raw score for yct-- see ReadSSConvTableForY() */
  double incp;      /* raw score inc for yct-- see ReadSSConvTableForY() */
  char nameyct[100];                  /* name of file containing yct[][] */
  double **yct;                                /* conversion table for Y */
  int round;                   /* if round = i, then round to i-th place */
  int lprss;                      /* lowest possible rounded scale score */
  int hprss;                     /* highest possible rounded scale score */
  int nrep;                      /* number of replications for bootstrap */
  int rep;     /* rep number for bootstrap; set to 0 for actual equating */
  struct BB_SMOOTH *bbx;  /* structure for beta binomial smoothing for x */
  struct BB_SMOOTH *bby;  /* structure for beta binomial smoothing for y */ 
  struct ULL_SMOOTH *ullx; /* structure for univ log-lin smoothing for x */ 
  struct ULL_SMOOTH *ully; /* structure for univ log-lin smoothing for y */
  struct BLL_SMOOTH *bllxv; /* struc for biv log-lin smoothing for x & v */ 
  struct BLL_SMOOTH *bllyv; /* struc for biv log-lin smoothing for y & v */
  struct BLL_SMOOTH *bllxy; /* struc for biv log-lin smoothing for x & y */
  struct CS_SMOOTH *cs;      /* structure for cubic-spline postsmoothing */
  struct IRT_INPUT *IRT_Input;                /* structure for IRT input */
};

struct ERAW_RESULTS{
  /* 
    equated raw-score results 
  */
  double msx[4];                         /* mean for x for synthetic pop */
  double msy[4];                         /* mean for y for synthetic pop */
  double ssx[4];                           /* sd for x for synthetic pop */
  double ssy[4];                           /* sd for y for synthetic pop */
  double gamma1[4];                                   /* gamma for pop 1 */
  double gamma2[4];                                   /* gamma for pop 2 */
  double a[4];                                                  /* slope */
  double b[4];                                              /* intercept */
  double **eraw;                                   /* equated raw scores */
  double **mts;                        /* moments for equated raw scores */
  double **fxs;     /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
  double **gys;     /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */ 
};

struct ESS_RESULTS{
  /* 
    equated scale-score results 
  */
  double **essu;                       /* unrounded equated scale scores */
  double **essr;                         /* rounded equated scale scores */
  double **mtsu;           /* moments for equated unrounded scale scores */
  double **mtsr;             /* moments for equated rounded scale scores */
};
                                                       
struct BOOT_ERAW_RESULTS{
  /* 
    equated raw-score results for bootstrap
  */  
  double **mn;                  /* for matrix of sum and mean for scores */
  double **sd;                   /* for matrix of sum2 and sd for scores */
  double *bse;                                 /* overall bootstrap se's */
};

struct BOOT_ESS_RESULTS{
  /* 
    equated scale-score results for bootstrap
  */
  double **mnu;   /* for matrix of sum and mn for unrounded scale scores */
  double **sdu;  /* for matrix of sum2 and sd for unrounded scale scores */
  double *bseu;     /* overall bootstrap se's for unrounded scale scales */ 
  double **mnr;   /* for matrix of sum and mean for rounded scale scores */
  double **sdr;    /* for matrix of sum2 and sd for rounded scale scores */
  double *bser;       /* overall bootstrap se's for rounded scale scales */
};

struct BB_SMOOTH{
  /* 
    beta-binomial smoothing: 2 or 4 parameter beta with 
    binomial or compound binomial error
  */
  int num_items;		                      /* number of items on test */
  int num_persons;	                                      /* sample size */
  int nparm;                            /* number of parameters (2 or 4) */
  double rel;                       /* reliability -- almost always Kr20 */
  double lordk;		  /* Lord's k for approximation of compound binomial */
  double beta[4];		        /* parameters of true score distribution */
  double rmoments[4];	                            /* Raw score moments */
  double fmoments[4];	                     /* Fitted raw score moments */
  double tmoments[4];	                          /* True score  moments */
  double lrchisq;		  /* likelihood ratio chi-square for fitted dist */
  double pchisq;		           /* Pearson chi-square for fitted dist */
  short  momentsfit;	                        /* number of moments fit */
  double *density;	   /* pointer to fitted raw score dist (proportions) */
  double *crfd;          /* pointer to cum rel freq dist for fitted dist */
  double *prd;         /*pointer to percentile rank dist for fitted dist */
};

struct ULL_SMOOTH{
  /* 
    univariate log-linear smoothing
  */
  int num_persons;                                  /* number of persons */
  int ns;                                  /* number of score categories */
  double min;                                       /* minimum raw score */
  double inc;                                 /* increment in raw scores */
  int c;                               /* number of degrees of smoothing */
  double **B_raw;                        /* design matrix for raw scores */
  double **B;                         /* design matrix used for solution */
  double *nct;                                      /* actual frqeuncies */
  double *mct;                                     /* fitted frequencies */
  double *Beta;                                          /* coefficients */
  double ap;                         /* normalizing constant used in CLL */
  double *n_mts;                                 /* actual moments for B */
  double *m_mts;                                 /* fitted moments for B */
  double *n_mts_raw;                   /* central moments for actual and */
  double *m_mts_raw;                   /* central moments for fitted and */
  int nit;                        /* number of iterations to convergence */
  int max_nit;                   /* maximum number of iterations allowed */
  int ctype;    /* comparison for criterion (0-->absolute; 1-->relative) */
  int Btype;     /* 0--> use B mts for crit; 1--> use B_raw and cent mts */
  int scale;      /* 0 --> no scaling for B design matrix; 1--> scaling  */
  double crit;                                  /* convergence criterion */
  double lrchisq;                         /* likelihood-ratio chi-square */
  int nzero;                                  /* # 0's (see iteration()) */
  double *density;	   /* pointer to fitted raw score dist (proportions) */
  double *crfd;          /* pointer to cum rel freq dist for fitted dist */
  double *prd;         /*pointer to percentile rank dist for fitted dist */
};

struct BLL_SMOOTH{
  /* 
    Bivariate log-linear smoothing, where bivariate distribution that is
    smoothed is (rows = non-common items) x (cols = common items);
	i.e. u by v.  nct[] is the row-major version of this matrix for 
	actual frequencies; mct[] is the corresponding matrix for fitted
	frequencies.  bfd[][]is for the x by v matrix of 
	fitted frequencies; subsequent elements are marginals for this 
	matrix.  For an external anchor, u by v is identical to x by v.

	Note that notation is in terms of x/u/v (implicityly for pop 1); 
	same logic applies to y/u/v (implcitly for pop 2). 
	
	For a single-group design, the "anchor" is effectively external 
	in the sense that x=u, v plays the role of y, and Wrapper_SL 
	will put x (rows) on scale of y (cols).   
  */
  int anchor;                          /* = 0 (external); = 1 (internal) */ 
  int num_persons;                                  /* number of persons */
  int nsu;        /* number of score categories for u (non-common items) */
  double minu;                                /* minimum raw score for u */
  double incu;                          /* increment in raw scores for u */
  int nsv;            /* number of score categories for v (common items) */
  double minv;                                /* minimum raw score for v */
  double incv;                          /* increment in raw scores for v */
  int ns;                            /* total number of score categories */
  int cu;                        /* number of degrees of smoothing for u */
  int cv;                        /* number of degrees of smoothing for v */
  int cuv;        /* number of (u,v) cross product moments for smoothing */
  int **cpm;                      /* cross-product moments for smoothing */
  int nc;                          /* number of columns in design matrix */
  double **B_raw;                        /* design matrix for raw scores */
  double **B;                         /* design matrix used for solution */
  double *nct;                                      /* actual frqeuncies */
  double *mct;                                     /* fitted frequencies */
  double *Beta;                                          /* coefficients */
  double ap;                         /* normalizing constant used in CLL */
  double *n_mts;                                 /* actual moments for B */
  double *m_mts;                                 /* fitted moments for B */
  double *n_mts_raw;           /* central moments for actual using B_raw */
  double *m_mts_raw;           /* central moments for fitted using B_raw */
  int nit;                        /* number of iterations to convergence */
  int max_nit;                   /* maximum number of iterations allowed */
  int ctype;    /* comparison for criterion (0-->absolute; 1-->relative) */
  int Btype;     /* 0--> use B mts for crit; 1--> use B_raw and cent mts */
  int scale;      /* 0 --> no scaling for B design matrix; 1--> scaling  */
  double crit;                                  /* convergence criterion */
  double lrchisq;                         /* likelihood-ratio chi-square */
  int nzero;                                  /* # 0's (see iteration()) */
  int nsx;    /* = nsu + nsv -1 if internal anchor; = nsu if ext anchor  */
  double minx;                                /* minimum raw score for x */
  double incx;                          /* increment in raw scores for x */
  double **bfd;      /* fitted biv freq dist for x by v (see note below) */
  double *fd_x;                   /* pointer to fitted frequencies for x */
  double *density_x;          /* pointer to fitted raw score props for x */
  double *crfd_x;    /* pointer to cum rel freq dist for fit. dist for x */
  double *prd_x;             /* pointer to PR dist for fitted dist for x */
  double *fd_v;                   /* pointer to fitted frequencies for v */
  double *density_v;    /* pointer to fitted raw score proportions for v */
  double *crfd_v;  /* pointer to cum rel freq dist for fitted dist for v */
  double *prd_v; /*pointer to percentile rank dist for fitted dist for v */

  /* Note: if external anchor, bfd[][] has dimensions nsu x nsv; 
           if internal anchor, bfd[][] has dimensions nsx x nsv, where
           nsx = nsu + nsv - 1 (bfd[][] includes structural zeros);
		   density[], crfd[], and prd[] are associated with rows of 
		   bfd[][]; density_v[], crfd_v[], and prd_v[] are associated 
		   with columns of bfd[][]; */ 
  
  double **brfd;                  /* fitted biv rel freq dist for x by v */
  double *crfd_vector_bfd;           /* cum rel fd as a row-major vector 
								        from bfd[][]; used for bootstrap */
};

struct CS_SMOOTH{
  /* 
    cubic-spline postsmoothing; need one of these structures in
	PDATA structure for putting x on scale of y, and one for
	PDATA structure for putting y on scale of x.
  */
  int ns;                                  /* number of score categories */
  double s;                      /* smoothing or "flexibility" parameter */
  double prlow;     /* percentile rank for lowest score that is smoothed */
  double prhigh;   /* percentile rank for highest score that is smoothed */
  int low;                 /* lowest (pseudo raw) score that is smoothed */
  int high;               /* highest (pseudo raw) score that is smoothed */
  int nsb;      /* high-low+1 = number of score categories in [low,high] */
  double *eeq;                    /* equipercentile equivalents: eeq[ns] */
  double *se;                                 /* standard errors: se[ns] */
  double *cmat;      /* vector containing a, b, c, d coeffs: cmat[4*ns];
					                       Note--dimensioned for maximum */
  double *eeqs;           /* cubic-spline smoothed equivalents including
			                               interpolated values: eeqs[ns] */
  double *inv; /* inverse of eeqs[]; computed only for Y to X; number of
			   score categories is ns for X, NOT ns for Y: inv[ns for X] */  
};

struct IRT_TRUE{
  /*
    IRT true score equating: 
	item parameters and other statitistics for a test form.
  */
  int num_items;		                      /* number of items on test */
  int num_scores;                     /* number of test score categories */
  int *num_icat;          /* pointer to numbers of item score categories */
  int *num_choice;         /* pointer to numbers of choices for MC items */
  double *right;                          /* poiner to right item scores */
  double *wrong;                          /* poiner to wrong item scores */
  double *a_par;                              /* pointer to a parameters */
  double *b_par;                              /* pointer to b parameters */
  double *c_par;                              /* pointer to c parameters */
  double **d_par;                 /* pointer to pointers to d parameters */ 
  double *scores;                              /* pointer to test scores */
  double *theta;                 /* theta values corresponding to scores */
  double min;                             /* minimal possible test score */
            /* min may not equal to the actually minimal score scores[0] */
  double max;                             /* maximum possible test score */
  double inc;                                         /* score increment */
  double chance;                          /* the chance level test score */
};

/************************ function prototypes ****************************/

int loc(double x, double min, double inc);
int nscores(double max, double min, double inc);
double score(int loc, double min, double inc);
int skipcols(FILE *fp, int k);
int atodouble(FILE *fp, char *str, int begin, int end, double *x);
int atointeger(FILE *fp, char *str, int begin, int end, int *x);
void flushline(FILE *fp);
void runerror(char error_text[]);
void convertFtoW(char fname[], int nv, int fields[][3], char tname[]);

void cum_rel_freqs(double min, double max, double inc, 
                   double *rfd, double *crfd);
double perc_rank(double min, double max, double inc, double *crfd, double x);
double interpolate(double x, int ns, double *f);                 

int ReadRawGet_moments(char fname[], int scol,
                       double *moments, double *mind, double *maxd);
int ReadRawGet_mind_maxd(char fname[], int scol, double *mind, double *maxd);
void ReadFdGet_USTATS(char fname[], int scol, int fcol, double min, 
					  double max, double inc, char id, struct USTATS *s);
void ReadRawGet_USTATS(char fname[], int scol, double min, 
                       double max, double inc, char id, struct USTATS *s);
void ReadRawGet_BSTATS(char fname[], int rows, int cols, 
                double rmin, double rmax, double rinc, 
                double cmin, double cmax, double cinc, 
                char rid, char cid, struct BSTATS *s);
void ReadSSConvTableForY(char nameyct[], double minp, double maxp, double inc,
                         double **yct);
int MomentsFromFD(double min, double max, double inc, double *scores,
                  int *fd, double *moments);
void MomentsFromRFD(double min, double max, double inc, double *scores,
                    double *fd, double *moments);
                     
void EquiEquate(int nsy, double miny, double incy, 
               double *crfdy, int nsx, double *prdx, double *eraw);
double perc_point(int ns, double min, double inc, 
                  double *crfd, double pr);
 
void Wrapper_ESS(struct PDATA *inall, struct ERAW_RESULTS *r, double minp, 
                double maxp, double incp, char *nameyct, int round, 
                int lprss, int hprss, struct ESS_RESULTS *s);
void Equated_ss(double min, double max, double inc, double minp, double maxp,
                double incp, double *eraw, double **yct, int round, int lprss,
                int hprss, double *essu, double *essr);
                
void Print_ESS(FILE *fp, char tt[], struct PDATA *inall, struct ESS_RESULTS *s);
void Print_USTATS(FILE *fp, char tt[], struct USTATS *s);
void Print_BSTATS(FILE *fp, char tt[], struct BSTATS *s, int pbfd);
void Print_vector(FILE *fp, char label[], double *vector, int nrows, 
                  char rowhead[], char colhead[]);
void Print_matrix(FILE *fp, char label[], double **matrix, int nrows, 
                  int ncols, char rowhead[], char colheads[]);
void Print_file(char FileName[], FILE *fp);

/* The next four lines of code provide an enumeration constant and 
   prototypes for three functions that have the same functionality as some 
   of the public-domain Numerical Recipes functions (see NRutilities.h and
   NRutilities.c) */

enum	datatype {CHARACTER, INTEGER, FLOAT,DOUBLE}  ;  
void	er_error2(char err_source[], char err_message[]); 
void	release_matrix(enum datatype mode, void **matrix, long rid_low, long cid_low);
void **	allocate_matrix(enum datatype mode,long rid_low, long rid_high, long cid_low, long cid_high); 

void    er_scale(double vect1[], double scale,double vect2[],int n,int offset);
double  er_dot(double vect1[], double vect2[], int n, int offset);
void    er_daxpy(double vectY[], double scale, double vectX[], int n, int offset);
void    er_r1update(double ** matx,double scale,double vect[],int n,int offset);
void    er_mvmult(double vect2[], double **matx,double vect1[], int n, int offset);

/* float ran2(long *idum);                                                    * 
 * ================= functions for random number generation ================  */  
float er_random(long *seed);

/*  void ludcmp(double **a, int n, int *indx, double *d);      
 *  void lubksb(double **a, int n, int *indx, double b[]);     
 * ================= functions for gaussj ==================================  */ 
void	er_ludcmp(double **a, int n, int *pivot, double *det); 
void    er_lubksb(double **a, int n, int pivot[], double b[]);  

/* void sort(unsigned long n, float arr[]);                    
 * ================= functions for gaussj ==================================  */  
void	er_sort(float* vector, int left, int right);   

/* void gaussj(double **a, int n, double **b, int m);          
 * er_matrix_inverse() replaces gaussj(), MatrixInverse0(),    
 * and MatrixInverse1().                                       
 * ================= functions for gaussj ==================================  */  
void    er_matrix_inverse(int n, double **a); 

/* void qrdcmp(double **a, int n, double *c, double *d,		  
 *		int *sing);										       
 * void genqrdcmp(double **a, int n, int m, double *c,	       
 *		double *d, int *sing); 
 * er_qrdcmp() replaces both qrdcmp() and genqrdcmp() 
 * ================= functions for gaussj ==================================  */  
void	er_qrdcmp(double **a, int nrow, int ncol, double *coeff, double *diag); 

/* void		dfpmin(double p[], int n, double gtol, int *iter, double *fret,
 *  		double(*func)(double []), void (*dfunc)(double [], double []));
 * void		lnsrch(int n, double xold[], double fold, double g[], double p[], 
 *			double x[], double *f, double stpmax, int *check, 
 *			double (*func)(double [])); 
 * double	rtsafe(void (*funcd)(double, double *, double *), double x1, 
 *			double x2, double xacc); 
 * ================= functions for dfpmin ==================================  */  
int     er_lnsrch(double xold[], int n, double gk[], double sk[],  double maxstep, 
			     double (*ftn_ptr)(double []), double xnew[]) ;
void    er_dfpmin(double xold[], int n, double error, int * numiter, double *fvalue, 
				  double(*ftn_ptr)(double []), void (*dftn_ptr)(double [], double [])); 
double  er_rtsafe(void (*funcd)(double, double *, double *), double x0, double x1, double error); 

#endif 

