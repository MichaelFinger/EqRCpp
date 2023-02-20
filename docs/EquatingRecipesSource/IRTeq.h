/*    
  Header file for IRT True Score Equating for Mixed-Format Tests

  Author: Seonghoon Kim, with contribtuions by
          Robert L. Brennan and Tianyou Wang

  Date of last revision: September 24, 2008   
*/

#include <stdio.h>
#include <stdlib.h>
#include "IRTst.h"

#ifndef _IRTEQUATE_H_
#define _IRTEQUATE_H_


struct RawFitDist {
    /*
      Structure that stores IRT fitted distributions * 
    */
	int nRaws;                                          /* Number of raw score categories */
	double *rawScrs;                                           /* Raw scores; zero-offset */
	double *newFits;                /* Fitted distribution for the new group; zero-offset */
	double *oldFits;                /* Fitted distribution for the old group; zero-offset */
	double *synFits;          /* Fitted distribution for the synthetic group; zero-offset */
	double mtsng[4];                      /* moments for fitted distribtion for new group */
	double mtsog[4];                      /* moments for fitted distribtion for old group */
	double mtssg[4];                /* moments for fitted distribtion for synthetic group */
};

struct RawTruObsEquiv {
    /*
      Structure that stores IRT true-score and observed-score equating results * 
    */
	int nRaws;                                 /* Number of new form raw score categories */
	double TCCnewMin;                   /* lower limit of true test score on the new form */
	double TCColdMin;                   /* lower limit of true test score on the old form */
	double *thetaTru;                                 /* Theta-equivalent of Form X score */
	double *unroundedEqTru;/* Unrounded raw-to-raw conversion for IRT true score equating */
	double *roundedEqTru;    /* Rounded raw-to-raw conversion for IRT true score equating */
	double *unroundedEqObs; /* Unrounded raw-to-raw conversion for IRT obs score equating */
	double *roundedEqObs;     /* Rounded raw-to-raw conversion for IRT obs score equating */
	double mtsTru[4];          /* new form equated (true score) score (unrounded) moments */
	double mtsObs[4];          /* new form equated (obs. score) score (unrounded) moments */
};

struct IRT_INPUT {
  /* 
    Structure that contains all input for IRT equating conducted using Wrapper_IRTeq(). 
    This struct is a member of the PDATA structure
  */
  char method;                /* 'T' for true score; 'O' for observed score; 'A' for both */
  int *NewFD;                               /* Actual frequency distribution for new form */
  char ItemNewFile[100];        /* Name of file containing item parameter estimates for X */
  char ItemOldFile[100];        /* Name of file containing item parameter estimates for Y */
  char DistNewFile[100];              /* Name of file containing theta distribution for X */
  char DistOldFile[100];              /* Name of file containing theta distribution for Y */
  struct ItemSpec *NewItems;
  struct ItemSpec *OldItems;
  struct IRTstControl *stControl;
  struct RawFitDist *NewForm;
  struct RawFitDist *OldForm;
  struct RawTruObsEquiv *RawEq;
};

/* memory Allocation functions */

void RawFitMem(struct ItemSpec *Items, const char *oldOrnew, 
		struct IRTstControl *Handle, struct RawFitDist *Form);
void RawEqResultsMem(struct ItemSpec *NewItems, struct IRTstControl *Handle,
		struct RawTruObsEquiv *RawEq);

/* memory Deallocation functions */

void RawFitDeAlloc(struct RawFitDist *Form);
void RawEqResultsDeAlloc(struct RawTruObsEquiv *RawEq);

/* IRT observed score equating functions */

short IRTmixObsEq(struct IRTstControl *Handle, struct ItemSpec *NewItems,
		struct ItemSpec *OldItems, double wNew, double wOld,
		struct RawFitDist *newForm, struct RawFitDist *oldForm,
		struct RawTruObsEquiv *RawEq);
short IRTmixObsDist(struct ItemSpec *Items, int n, int MaxScrP, int nq, double xqpts[],
		double xqwts[], int *nscr, double xscr[], double xmarg[]);
short ObsDistGivenTheta(double theta, struct ItemSpec *Items, int n,
		int MaxCat, int MaxScrP, int *nscr, double xscr[], double xnew[]);
short recurs(int mino, int maxo, double xold[],
            int mitem, double iitem[], double xitem[],
            int * minn, int * maxn, double xnew[]);

/* IRT true score equating functions */

short trueScoreEq(struct IRTstControl *Handle, struct ItemSpec *NewItems,
		struct ItemSpec *OldItems, int nScores, const double *newScores, 
		double *eqvOld, double *theta, double *newMin, double *OldMin);
double f_mix(double theta);
double f_mixDer(double theta);
void funcd(double x, double *f, double *fd);
double trueScore(struct ItemSpec *Items, int n, double theta);


/* functions associated with different IRT models */

double ProbCCC(struct ItemSpec *Item, int CatID, double theta);
double PdCCCoverTheta(struct ItemSpec *Item, int CatID, double theta);
double Pd3PLoverTheta(int CatID, double theta, double D, double a, double b, double c);
double PdLGRoverTheta(int CatNum, int CatID, double theta, double D, double a, double b[]);
double PdGPCoverTheta(int CatNum, int CatID, double theta, double D, double a, double b[]);
double PdNRMoverTheta(int CatNum, int CatID, double theta, double a[], double c[]);


/* Wrapper and associated print functions */

 void Wrapper_IRTeq(char design, char method, double w1, 
     char ItemNewFile[], char ItemOldFile[], char DistNewFile[], char DistOldFile[], 
     struct ItemSpec *NewItems, struct ItemSpec *OldItems, 
     struct RawFitDist *NewForm, struct RawFitDist *OldForm,
     struct RawTruObsEquiv *RawEq, struct IRTstControl *StInfo, int *NewFD, 
     struct IRT_INPUT *irtall, struct PDATA *pinall, struct ERAW_RESULTS *r);  

void Print_IRTeq(FILE *fp, char tt[], struct PDATA *pinall, 
                 struct ERAW_RESULTS *r, int PrintFiles);
void Print_ESS_QD(FILE *fp, char tt[], struct PDATA *pinall, struct ESS_RESULTS *s,
                  struct RawFitDist *NewForm, struct RawFitDist *OldForm);                  

#endif
