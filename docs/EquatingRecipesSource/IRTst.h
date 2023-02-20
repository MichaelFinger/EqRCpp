/* September 5, 2005
  
   Seonghoon Kim
   ACT, Inc.

*/

#include "ERutilities.h"
#ifndef _IRTST_H_
#define _IRTST_H_

enum OnOff {off=0, on};
enum symmetry {old_scale, new_scale, symmetric};
enum ModelSpec {l3=1, gr, pc, nr};
enum LossSpec {mix_ha, mix_sl};

struct ItemSpec {
	int ItemID;
	int CatNum;
	double ScaleConst;
	double *ScoreFunc;
	double *a;
	double *b;
	double *c;
	double *d;
	enum ModelSpec model;
};

struct CommonItemSpec { /* structure for specification of common items */
	int NewID;
	int OldID;
	int CatNum;
	double ScaleConst;
	double *ScoreFunc;
	double *Na;
	double *Nb;
	double *Nc;
	double *Nd;	
	double *Oa;
	double *Ob;
	double *Oc;
	double *Od;
	enum ModelSpec model;
};
   
struct IRTstControl {
	int NewItemNum; 
	int OldItemNum;  /* Numbers of items on the new (X) and old (Y) forms */
	int ComItemNum;  /* Number of the common items */
	int NewThetaNum; 
	int OldThetaNum; /* Numbers of ability points on the new and old scales */
	double NewRawMin;
	double NewRawMax;
	double NewRawInc;
	double OldRawMin;
	double OldRawMax;
	double OldRawInc;
	double *NewThetaValues;
	double *NewThetaWeights;
	double *OldThetaValues;
	double *OldThetaWeights;
};


/* Functions for IRT scale transformation */

/* ItemInfoR.c */
struct ItemSpec *ItemInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle);


/* ItemPairsR.c */
struct CommonItemSpec *ItemPairsRead(FILE *inf, struct ItemSpec *NewItem,
		struct ItemSpec *OldItem, struct IRTstControl *Handle);


/* ProbDeriv.c */
double ProbOld(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I);
double ProbNew(struct CommonItemSpec *Item, int CatID, double theta, enum OnOff original,
		double S, double I);
double PdOldOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I);
double PdOldOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I);
double PdNewOverS(struct CommonItemSpec *Item, int CatID, double theta, double S, double I);
double PdNewOverI(struct CommonItemSpec *Item, int CatID, double theta, double S, double I);
double Prob3PL(int CatID, double theta, double D, double a, double b, double c,
		char scale[], double S, double I);
double Pd3PLOldOverS(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I);
double Pd3PLOldOverI(int CatID, double theta, double D, double na, double nb, double nc,
		double S, double I);
double Pd3PLNewOverS(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I);
double Pd3PLNewOverI(int CatID, double theta, double D, double oa, double ob, double oc,
		double S, double I);
double CumProbLGR(int CatID, double theta, double D, double a, double b,
		char scale[], double S, double I);
double ProbLGR(int CatNum, int CatID, double theta, double D, double a, double b[],
               char scale[], double S, double I);
double PdLGROldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I);
double PdLGROldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I);
double PdLGRNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I);
double PdLGRNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I);
double ProbGPC(int CatNum, int CatID, double theta, double D, double a, double b[],
		char scale[], double S, double I);
double PdGPCOldOverS(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I);
double PdGPCOldOverI(int CatNum, int CatID, double theta, double D, double na, double nb[],
		double S, double I);
double PdGPCNewOverS(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I);
double PdGPCNewOverI(int CatNum, int CatID, double theta, double D, double oa, double ob[],
		double S, double I);
double ProbNRM(int CatNum, int CatID, double theta, double a[], double c[],
		char scale[], double S, double I);
double PdNRMOldOverS(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I);
double PdNRMOldOverI(int CatNum, int CatID, double theta, double na[], double nc[],
		double S, double I);
double PdNRMNewOverS(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I);
double PdNRMNewOverI(int CatNum, int CatID, double theta, double oa[], double oc[],
		double S, double I);


/* ScaleTransform.c */
void ScaleTransform(const char *ItemOutF, const char *DistOutF, double slope,
		double intercept, struct ItemSpec *NewItem, struct IRTstControl *Handle);


/* StHaebara.c */
void StHaebara(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept);
double FuncHaebara(double x[]);
void GradHaebara(double x[], double grad[]);


/* StMeanMean.c */
void StMeanMean(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept);


/* StMeanSigma.c */
void StMeanSigma(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		double *slope, double *intercept);


/* StMemDeAlloc.c */
void StItemDeAlloc(struct ItemSpec *Item, const char *oldOrnew, struct IRTstControl *Handle);
void StComItemDeAlloc(struct CommonItemSpec *ComItem, struct IRTstControl *Handle);
void StContDeAlloc(struct IRTstControl *Handle);


/* StStockingLord.c */
void StStockingLord(struct IRTstControl *Handle, struct CommonItemSpec *ComItem,
		enum symmetry SYM, enum OnOff FuncStd, double S0, double I0,
		double *slope, double *intercept);
double FuncStockingLord(double x[]);
void GradStockingLord(double x[], double grad[]);


/* ThetaInfoR.c */
void ThetaInfoRead(FILE *inf, const char *oldOrnew, struct IRTstControl *Handle);

/* Wrapper */

 void Wrapper_IRTst(FILE *outf, char tt[],
     char ItemNewFile[], char ItemOldFile[], char ItemCommonFile[], 
     char DistNewFile[], char DistOldFile[], 
     int HA, enum symmetry HAsym, enum OnOff HAfs, double HAs, double HAi,
     int SL, enum symmetry SLsym, enum OnOff SLfs, double SLs, double SLi,
     char ST[], int PrintFiles);

#endif
