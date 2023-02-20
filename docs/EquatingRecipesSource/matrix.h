/*
  matrix.h
*/

#ifndef MATRIX_H
#define MATRIX_H

void MatrixMultiVector0(int nrow, int ncol, double **m, double *v, double *r);
void MatrixMultiVector1(int nrow, int ncol, double **m, double *v, double *r);


void VectorMultiMatrix0(int nrow, int ncol, double *v, double **m, double *r);
void VectorMultiMatrix1(int nrow, int ncol, double *v, double **m, double *r);

void MatrixMultiMatrix0(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r);
void MatrixMultiMatrix1(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r);

double VectorNorm0(int ncat, double *v);

void PTransposeMultiP1(int nr, int nc, double **P, double **PTP);

void PTransposeMultiP0(int nr, int nc, double **P, double **PTP);

double Determinant(double **a, int n);
/* ===== Before version 1.0 update ===== 
void MatrixInverse1(int n, double **m); 
void MatrixInverse0(int n, double **m);
*/ 
void PtAP0(int nr, int nc, double **P, double **A, double **PtAP);

double vtAv0(int n, double *v, double **A);

#endif 

