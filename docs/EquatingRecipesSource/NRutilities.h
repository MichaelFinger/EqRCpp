/*
  NRutilities.h
  
  Header for selected utility functions from Numerical Recipes.
  These are in the public domain
*/ 

#ifndef NRUTILITIES_H
#define NRUTILITIES_H

#pragma warning(disable:4996)
#pragma warning(disable:4706)
#pragma warning(disable:4715)
#pragma warning(disable:4701)
#pragma warning(disable:4127)
#pragma warning(disable:4244)
#pragma warning(disable:4100)
#pragma warning(disable:4305)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(char error_text[]);

unsigned char *cvector(long nl, long nh);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
char **cmatrix(long nrl, long nrh, long ncl, long nch);

void free_cvector(unsigned char *v, long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_cmatrix(int **m, long nrl, long nrh, long ncl, long nch);

#endif 
