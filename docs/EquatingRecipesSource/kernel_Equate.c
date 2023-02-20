/*
kernel_Equate.c   code for kernel equating

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

*/

#include <math.h>
#include <stdlib.h> 
#include "NRutilities.h"
#include "LogLinear.h"
#include "CG_EquiEquate.h"
#include "kernel_EQuate.h"
#include "ERutilities.h"

/*--------------------------------------------------------------------------
	KernelContinuPdf
	
	functionality

	Computes kernel smoohxing function from von Davier, Holland, & Thayer 
	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 10/29/2004.
	
	input
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the relative frequency distribution
		hx          bandwidhx for the kernel smoohxing
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		The function returns the kernel pdf
--------------------------------------------------------------------------*/
double KernelContinuPdf(int ncat, double *scores, double *fd, double hx, double x)
{
	int i;
	double zscore;                 /*the standardized normal random variable */
	double mu=0.0, sigma=0.0;      /*the mean and standard deviation of the
							                         unsmoothed distribution */
	double ax, sumsq=0.0;                            /*ax, and second moment */
	double temp1;                                       /*temporary variable */
	double pdf=0;


	for (i=0; i<ncat; i++) {
		mu += scores[i] * fd[i];
		sumsq += scores[i] * scores[i] * fd[i];
	}

	sigma = sumsq - mu * mu;
	ax = sqrt(sigma / (sigma + hx * hx));

	for (i=0; i<ncat; i++) {
		zscore = (x - ax * scores[i] - (1 - ax) * mu) / ax / hx;
		temp1 = StdNormalPdf(zscore);
		pdf += fd[i] * temp1 / ax / hx;
	}

	return pdf;

}

/*--------------------------------------------------------------------------
	KernelContinuCdf
	
	functionality:

	Computes kernel smoohxing function from von Davier, Holland, & Thayer 
	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 10/29/2004.
	
	input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the relative frequency distribution
		hx          bandwidhx for the kernel smoohxing
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output:
 		The function returns the kernel cdf
--------------------------------------------------------------------------*/
double KernelContinuCdf(int ncat, double *scores, double *fd, double hx, double x)
{
	int i;
	double zscore;               /*the standardized normal random variable */
	double mu=0.0, sigma=0.0;     /*the mean and standard deviation of the
							                     unsmoothed distribution   */
	double ax, sumsq=0.0;                          /*ax, and second moment */
	double temp2;                                     /*temporary variable */
	double cdf=0;

	for (i=0; i<ncat; i++) {
		mu += scores[i] * fd[i];
		sumsq += scores[i] * scores[i] * fd[i];
	}

	sigma = sumsq - mu * mu;
	ax = sqrt(sigma / (sigma + hx * hx));

	for (i=0; i<ncat; i++) {
		zscore = (x - ax * scores[i] - (1 - ax) * mu) / ax / hx;
		temp2 = StdNormalCdf(zscore);
		cdf += fd[i] * temp2;
	}

	return cdf;

}

/*--------------------------------------------------------------------------
	StdNormalCdf
	
	functionality

	Computes the cdf for a z score from a standard normal distribution
	
	author: Tianyou Wang 10/29/2004.
	
	input
		x           a z score        

 
	output
 		the function returns the normal cdf for x.
--------------------------------------------------------------------------*/

double StdNormalCdf(double x) 
{
	double sum = 0.0;
	double t, z;
     
	t = 1 / (1 + .2316419 * fabs(x));
	z = 0.398942280401433 * exp(-.5 * x * x) ;
	sum = 1 - z * t* (.319381530  + t * (-.356563782  
		  + t * (1.781477937 + t * (-1.821255978 
		  + t * 1.330274429))));
	
	if (x>=0) return sum;
	else return 1-sum;

}

/*--------------------------------------------------------------------------
	StdNormalPdf
	
	functionality

	Computes the pdf for a z score from a standard normal distribution
	
	author: Tianyou Wang 10/29/2004.
	
	input
		x           a z score        

 
	output
 		the function returns the normal cdf for x.
--------------------------------------------------------------------------*/
double StdNormalPdf(double x)
{

	return 0.398942280401433 * exp(-.5 * x * x);
	
}

/*--------------------------------------------------------------------------
	NormalPdf
	
	functionality

	Computes the pdf for a z score from a normal distribution
	
	author: Tianyou Wang 10/29/2004.
	
	input
		x           a z score        

 
	output
 		the function returns the normal cdf for x.
--------------------------------------------------------------------------*/
double NormalPdf(int npara, double *para, double x)
{
	double mu, sigma, pdf;
	if (npara==2){
		mu =  para[0];
		sigma = para[1];
	}

	pdf = 0.398942280401433 * exp(-.5 * pow((x - mu) / sigma, 2)) / sigma;
	
	return pdf;
}

/*------------------------------------------------------------------------------
	CalcKernelMoments	
	
	functionality

	calculates mean, sd, skewness and kurtosis for a continuous distribution.

	author: Tianyou Wang 10/29/2004.

	input -
        (*pdf)()    pointer to a function which is the pdf of the continuous 
		            distribution
        a   		lower limit of distribution
        b   		upper limit of distribution
		npara       number of parameters for the distribution
		para        vector of parameters for the distribution

	output -
		moments - mean, sd, skewness and kurtosis of distribution

------------------------------------------------------------------------------*/

void CalcKernelMoments(PTR_FTN1 pdf, double a, double b, int ncat, double *score, 
					   double *fd, double hx, double *moments)
{
	int j;
	double xr, xm, dx, ss, s1, s2, s3, s4;
	static double x[]={0.02435029266342, 0.07299312178780, 0.12146281929612,
		0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
		0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705, 
		0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466, 
		0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334, 
		0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511, 
		0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205, 
		0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196, 
		0.99930504173577};
	static double w[]={0.04869095700914, 0.04857546744150, 0.04834476223480, 
		0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131, 
        0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365, 
        0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024, 
        0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240, 
		0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
        0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
        0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056, 
        0.00178328072170};

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);

	ss = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		ss += w[j] * (pdf(ncat, score, fd, hx, (xm + dx))  + 
			pdf(ncat, score, fd, hx, (xm - dx)));
	}
	ss *= xr;

	s1 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s1 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * (xm + dx) + 
			pdf(ncat, score, fd, hx, (xm - dx)) / ss * (xm-dx));
	}
	s1 *= xr;

	s2 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s2 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 2) + 
			pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 2));
	}
	s2 *= xr;

	s3 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s3 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 3) + 
			pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 3));
	}
	s3 *= xr / pow(s2, 1.5);

	s4 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s4 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 4) + 
			pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 4));
	}
	s4 *= xr / pow(s2, 2);

	moments[0] = s1;
	moments[1] = sqrt(s2);
	moments[2] = s3;
	moments[3] = s4;
}	

/*****************************************************************************
  Pen1

  Functionality:

  Compute the Penality function PEN1 of von Davier, Holland, & Thayer 
  (2004, p. 62).

  author: Tianyou Wang 1/5/2005.

  Input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
		hx          bandwidhx for the kernel smoohxing

  Output:
        Pen1        The panelty function PEN1 value
*****************************************************************************/
double Pen1(int ncat, double *scores, double *fd, double hx) 
{
	int i;
	double pdf, sum=0;

	for (i=0; i<ncat; i++) {
		pdf = KernelContinuPdf(ncat, scores, fd, hx, scores[i]);
		sum += (pdf - fd[i]) * (pdf - fd[i]);
	}

	return sum;
}

/*****************************************************************************
  Pen2

  Functionality:

  Compute the Penality function PEN2 of von Davier, Holland, & Thayer 
  (2004, p. 63).

  author: Tianyou Wang 1/5/2005.

  Input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
		hx          bandwidhx for the kernel smoohxing

  Output:
        Pen1        The panelty function PEN1 value
*****************************************************************************/

double Pen2(int ncat, double *scores, double *fd, double hx) 
{
	int i;
	double pdfPL, pdfPR, sum=0, left, right, A, B;

	for (i=0; i<ncat; i++) {
        left = scores[i] - .25;
		right = scores[i] + .25;
		pdfPL = KernelPdfDerivative(ncat, scores, fd, hx, left);
		pdfPR = KernelPdfDerivative(ncat, scores, fd, hx, right);
		A = 0; B = 1;
		if (pdfPL < 0) A = 1;
		if (pdfPR > 0) B = 0;
		sum += A * (1 - B);
	}

	return sum;
}

/*--------------------------------------------------------------------------
	KernelPdfDerivative
	
	functionality

	Computes the derivative of kernel pdf function from von Davier, Holland, & 
	Thayer (2004, p. 63). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 1/5/2005.
	
	input
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		inc         Increment between consecutive raw scores for old form.
		fd          vector containing the relative frequency distribution
		hx          bandwidhx for the kernel smoohxing
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		The function returns the kernel pdf
--------------------------------------------------------------------------*/
double KernelPdfDerivative(int ncat, double *scores, double *fd, double hx, double x)
{
	int i;
	double zscore;              /*the standardized normal random variable */
	double mu=0.0, sigma=0.0; /*the mean and standard deviation of the
							    unsmoothed distribution                 */
	double ax, sumsq=0.0;         /*ax, and second moment */
	double temp1;      /*temporary variable */
	double pdfDer=0;

	for (i=0; i<ncat; i++) {
		mu += scores[i] * fd[i];
		sumsq += scores[i] * scores[i] * fd[i];
	}

	sigma = sumsq - mu * mu;
	ax = sqrt(sigma / (sigma + hx * hx));

	for (i=0; i<ncat; i++) {
		zscore = (x - ax * scores[i] - (1 - ax) * mu) / ax / hx;
		temp1 = StdNormalPdf(zscore);
		pdfDer += fd[i] * temp1 / ax / ax / hx / hx * zscore;
	}

	return -pdfDer;

}

/*****************************************************************************
  Pen

  Functionality:

  Compute the combined Penality function of von Davier, Holland, & Thayer 
  (2004, p. 63).

  author: Tianyou Wang 1/5/2005.

  Input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
		hx          bandwidhx for the kernel smoohxing
		K           The weight for Pen2

  Output:
        return the combined panelty function value
*****************************************************************************/
double Pen(int ncat, double *scores, double *fd, double hx, double K) 
{
	double PEN, PEN1, PEN2;

	PEN1 = Pen1(ncat, scores, fd, hx);
	PEN2 = Pen2(ncat, scores, fd, hx);
	PEN = PEN1 + K * PEN2;

	return PEN;
}

/*****************************************************************************
  Optimalh

  Functionality:
        Find the optimal bandwidhx parameter for kernel continuization hx based on 
        von Davier, Holland, & Thayer (2004, p. 63). The Kernel Method of Test Equating.

  author: Tianyou Wang 1/5/2005.

  Input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		inc         Increment between consecutive raw scores for old form.
		fd          vector containing the frequency distribution
		K           The weight for PEN2
  
  Output:
        return the optimal hx
*****************************************************************************/
double Optimalh(int ncat, double *scores, double *fd, double K)
{
    int niter = 0;
    double eps = .0001, hxl = .0001, hxu = 3, hxlplus = .0002;
	double hxuminxs = 2.9, hxb =1.1, hxx;
    double pl, pu, plplus, puminxs, pb, px, optimhx, absdif;
     
    pl = Pen(ncat, scores, fd, hxl, K); 
    pu = Pen(ncat, scores, fd, hxu, K); 
    plplus = Pen(ncat, scores, fd, hxlplus, K); 
    puminxs = Pen(ncat, scores, fd, hxuminxs, K); 
    pb = Pen(ncat, scores, fd, hxb, K); 
    if (pl<pb && pb<pu && pl<plplus)
        optimhx = hxl;
    else if (pu<pb && pb<pl && pu<puminxs)
        optimhx = hxu;
    else
    {
        do
        {
            niter++;
            hxb = .38197 * hxu + .61803 * hxl;
            hxx = .61803 * hxu + .38197 * hxl;
            absdif = fabs(hxu - hxl);
            if (absdif < eps)
                break;
            pb = Pen(ncat, scores, fd, hxb, K); 
            px = Pen(ncat, scores, fd, hxx, K); 

			if (px<=pb)
                hxl = hxb;
            if (px>pb) 
                hxu = hxx;
        } while (niter <=200);
        optimhx = .5 * (hxb + hxx);        
    }

	return optimhx;
}

/*--------------------------------------------------------------------------
	KernelEquate
	
	functionality:

	Computes kernel equating function based on continuized cdf in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 1/5/2005.
	
	input:
		ncatx   Number of discrete score categories for the new form
		scoresx     vector containing the discrete scores for the new form
		fdx         vector containing the relative frequency distribution  for the new form
		hx          bandwidhx for the kernel smoohxing  for the new form
		ncaty   Number of discrete score categories for the old form
		scoresy     vector containing the discrete scores for the old form
		fdy         vector containing the relative frequency distribution  for the old form
		hy          bandwidhx for the kernel smoohxing  for the old form

 
	output:
 		Equatedx   a vector containing the equated score 
--------------------------------------------------------------------------*/
void KernelEquate(int ncatx, double *scoresx, double *fdx, double hx, 
					int ncaty, double *scoresy, double *fdy, double hy,
					double *Equatedx)
{
	int i;
	double cdfx;

	for (i=0; i<ncatx; i++) {
		cdfx = KernelContinuCdf(ncatx, scoresx, fdx, hx, scoresx[i]);
//		cdfy = KernelContinuCdf(ncaty, scoresy, fdy, hy, scoresy[i]);
		Equatedx[i] = KernelInverseCdf(ncaty, scoresy, fdy, hy, cdfx);
	}

}

/*--------------------------------------------------------------------------
	KernelInverseCdf
	
	functionality:

	Computes the inverse of the cdf in von Davier, Holland, & Thayer 
	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 1/5/2005.
	
	input:
		ncat    Number of discrete score categories
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
		h           bandwidhx for the kernel smoohxing
		cdf         a particular cdf for which the score is found
 
	output:
 		The function returns the inverse of cdf
--------------------------------------------------------------------------*/
double KernelInverseCdf(int ncat, double *scores, double *fd, double h, double cdf)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .000001;

	lb = scores[0] - 5.0;
	ub = scores[ncat-1] + 5.0;
	cdfl = KernelContinuCdf(ncat, scores, fd, h, lb);
	cdfu = KernelContinuCdf(ncat, scores, fd, h, ub);
	if (cdf<cdfl) 
		return scores[0];
	else if (cdf>cdfu)
		return scores[ncat-1];
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = KernelContinuCdf(ncat, scores, fd, h, half);  
            absdif = fabs(cdf - cdfhalf);
            if (absdif < eps)
                return half;
			else if (cdfhalf<cdf) 
				lb = half;
			else
				ub = half;
        } while (niter <=200);
	}

	return half;
}

/*--------------------------------------------------------------------------
	ComputeCmatrix
	
	functionality: 

	Computes the C matrix in von Davier, Holland, & Thayer 
	(2004, Equation 3.10). The Kernel Method of Test Equating.
	call the numerical recipe function qrdcmp to do QR decomposition.
	
	author: Tianyou Wang 3/11/2005.
	
	input:
		ncat    Number of discrete score categories
		degree      The highest polynomial degree of the log-linear model
		np          sample size
		B           B matrix (design matrix)
		fd          vector containing the frequency distribution
 
	output:
 		Cr          C matrix in von Davier et. al. Equation 3.10
--------------------------------------------------------------------------*/

void ComputeCmatrix(int ncat, int degree, long np, double **B, double *fd, double **Cr)
{
	int i, j, k, l;
	double *D, **A, **a, **qt, **q, **r, rrij;
	double *c, *d; 
	double con;
    FILE *outf;
    char outfname[]="Bmatrix.out";


	D = dvector(0, ncat-1);
	A = dmatrix(0, ncat-1, 0, ncat-1);
	a = dmatrix(1, ncat, 1, ncat);
	qt = dmatrix(0, ncat-1, 0, ncat-1);
	q = dmatrix(0, ncat-1, 0, ncat-1);
	r = dmatrix(0, ncat-1, 0, ncat-1);
	c = dvector(1, ncat);
	d = dvector(1, ncat);

	/* initialize all the matrices */
	for (i=0; i<ncat; i++)	{
		D[i] = 0;
		for (j=0; j<ncat; j++)	{
			A[i][j] = 0;
			a[i+1][j+1] = 0;
			q[i][j] = 0;
			qt[i][j] = 0;
			r[i][j] = 0;
		}
	}

	for (i=0; i<ncat; i++)	
		for (j=0; j<degree; j++) Cr[i][j] = 0;

	for (i=1; i<=ncat; i++)	{
		c[i] = 0;
		d[i] = 0;
	}

	for (i=0; i<ncat; i++)	{
		for (j=0; j<ncat; j++)	{
			D[i] = sqrt(fd[i]);
			rrij = sqrt(fd[i]) * fd[j];
			if (i==j) 
				A[i][j] = D[i] - rrij;
			else
				A[i][j] = - rrij;

		}
	}

	outf = fopen(outfname, "w");
	fprintf(outf, "B matrix \n ");
	for(i=0; i<degree; i++) {
		for(j=0; j<ncat; j++) {
			fprintf(outf, "%10.6e ", B[i][j]);
		}
		fprintf(outf, "\n ");
	}
	fprintf(outf, "\n ");


	for (i=0; i<ncat; i++)	{
		for (j=0; j<ncat; j++)	{
			for (k=0; k<ncat; k++)	{
				a[i+1][j+1] += A[i][k] * B[j][k];
			}
		}
	}

	
	fprintf(outf, "a matrix before decomposition 1\n ");
	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, " %13.10e ", a[i+1][j+1]);
		}
		fprintf(outf, "\n\n ");
	}
	fprintf(outf, "\n\n");

	/* ===== Before version 1.0 update =====
	qrdcmp(a, ncat, c, d, &sing);
	*/
	er_qrdcmp(a, ncat, ncat, c, d);
	/* ===== End of version 1.0 update ===== */ 
	
	
	fprintf(outf, "a matrix after decomposition 1\n ");
	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, " %13.10e ", a[i+1][j+1]);
		}
		fprintf(outf, "\n\n ");
	}
	fprintf(outf, "\n\n");

    /* compute the Q and R matrices */
    for (k=0;k<ncat;k++) {
		for (l=0;l<ncat;l++) {
			if (l > k) {
				r[k][l]=a[k+1][l+1];
                q[k][l]=0.0;
			} 
			else if (l < k) {
				r[k][l]=q[k][l]=0.0;
			} 
			else {
				r[k][l]=d[k+1];
                q[k][l]=1.0;
			} 
        } 
	 } 
		
	for (i=ncat-2;i>=0;i--) {
		for (con=0.0,k=i;k<ncat;k++) 
			con += a[k+1][i+1]*a[k+1][i+1];
			con /= 2.0;
            for (k=i;k<ncat;k++) {
				for (l=i;l<degree;l++) {
					qt[k][l]=0.0;
					for (j=i;j<ncat;j++) {
						qt[k][l] += q[j][l]*a[k+1][i+1]*a[j+1][i+1]/con;
					}
				}
            }
		for (k=i;k<ncat;k++) 
			for (l=i;l<degree;l++) {
				q[k][l] -= qt[k][l];
			}
	}

	fprintf(outf, "q matrix \n ");
	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, " %13.10e  ", q[i][j]);
		}
		fprintf(outf, "\n ");
	}
	fprintf(outf, "\n ");

	/* compute the Cr matrix */
	for(i=0; i<ncat; i++) 	
		for(j=0; j<degree; j++) 
			Cr[i][j] += (1.0/sqrt((double)np))*D[i]*q[i][j];

	/* output Cr matrix */
	fprintf(outf, "Cr matrix \n ");
	for(i=0; i<ncat; i++) 	{
		for(j=0; j<degree; j++) {
				fprintf(outf, " %13.10e ", Cr[i][j]*1000);
		}
		fprintf(outf, "\n ");

	}

	fclose(outf);
	free_dvector(D, 0, ncat-1);
	free_dmatrix(A, 0, ncat-1, 0, ncat-1);
	free_dmatrix(a, 1, ncat, 1, ncat);
	free_dmatrix(qt, 0, ncat-1, 0, ncat-1);
	free_dmatrix(q, 0, ncat-1, 0, ncat-1);
	free_dmatrix(r, 0, ncat-1, 0, ncat-1);
	free_dvector(c, 1, ncat);
	free_dvector(d, 1, ncat);

}

/*--------------------------------------------------------------------------
	ComputeCmatrixGen
	
	functionality: 

	Computes the C matrix in von Davier, Holland, & Thayer 
	(2004, Equation 3.10) for a general case where non-square matrixes 
	are involved. The Kernel Method of Test Equating.
	call the numerical recipe function qrdcmp to do QR decomposition.
	
	author: Tianyou Wang 3/11/2005.
	
	input:
		ncat    Number of discrete score categories
		degree      The highest polynomial degree of the log-linear model
		np          sample size
		B           B matrix (design matrix)
		fd          vector containing the frequency distribution
 
	output:
 		Cr          C matrix in von Davier et. al. Equation 3.10
--------------------------------------------------------------------------*/

void ComputeCmatrixGen(int ncat, int degree, long np, double **B, double *fd, double **Cr)
{
	int i, j, k, l;
	double *D, **A, **a, **b, **qt, **q, **r, rrij, sum;
	double *c, *d; 
	double con;
    FILE *outf;
    char outfname[]="Bmatrix.out";
	static int ii=0;


	D = dvector(0, ncat-1);
	A = dmatrix(0, ncat-1, 0, ncat-1);
	a = dmatrix(1, ncat, 1, degree);
	b = dmatrix(1, ncat, 1, degree);
	qt = dmatrix(0, ncat-1, 0, degree-1);
	q = dmatrix(0, ncat-1, 0, degree-1);
	r = dmatrix(0, ncat-1, 0, degree-1);
	c = dvector(1, ncat);
	d = dvector(1, ncat);

	/* initialize all the matrices */
	for (i=0; i<ncat; i++)	{
		D[i] = 0;
		c[i+1] = 0;
		d[i+1] = 0;
		for (j=0; j<ncat; j++)	A[i][j] = 0;
	}

	for (i=0; i<ncat; i++)	{
		for (j=0; j<degree; j++)	{
			a[i+1][j+1] = 0;
			q[i][j] = 0;
			qt[i][j] = 0;
			r[i][j] = 0;
			Cr[i][j] = 0;
		}
	}


	for (i=0; i<ncat; i++)	{
		D[i] = sqrt(fd[i]);
		for (j=0; j<ncat; j++)	{
			rrij = sqrt(fd[i]) * fd[j];
			if (i==j) 
				A[i][j] = D[i] - rrij;
			else
				A[i][j] = - rrij;

		}
	}

	outf = fopen(outfname, "w");
	fprintf(outf, "B matrix \n ");
	for(i=0; i<degree; i++) {
		for(j=0; j<ncat; j++) {
			fprintf(outf, "%1.0lf ", B[i][j]);
		}
		fprintf(outf, "\n\n ");
	}
	fprintf(outf, "\n\n ");


	for (i=0; i<ncat; i++)	{
		for (j=0; j<degree; j++) {
			for (k=0; k<ncat; k++)	{
				a[i+1][j+1] += A[i][k] * B[j][k];
			}
			b[i+1][j+1] = a[i+1][j+1];
		}
	} 

	/* ===== Before version 1.0 update =====
    genqrdcmp(a, ncat, degree, c, d, &sing);  
	*/  

	er_qrdcmp(a, ncat, degree, c, d);
	/* ===== End of version 1.0 update ===== */ 

	fprintf(outf, "a matrix after decomposition 2\n ");
	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, "%10.8lf ", a[i+1][j+1]);
		}
		fprintf(outf, "\n ");
	}
	fprintf(outf, "\n ");


	 /* compute the Q and R matrices */
     for (k=0;k<ncat;k++) {
		for (l=0;l<degree;l++) {
			if (l > k) {
				r[k][l]=a[k+1][l+1];
                q[k][l]=0.0;
			} 
			else if (l < k) {
				r[k][l]=q[k][l]=0.0;
			} 
			else {
				r[k][l]=d[k+1];
                q[k][l]=1.0;
			} 
        } 
	 } 

	/* in the next line, original i=ncat-2, now i=degree-1. i=degree-2 does not work */	
	for (i=degree-1;i>=0;i--) {
		for (con=0.0,k=i;k<ncat;k++) 
			con += a[k+1][i+1]*a[k+1][i+1];
			con /= 2.0;
            for (k=i;k<ncat;k++) {
				for (l=i;l<degree;l++) {
					qt[k][l]=0.0;
					for (j=i;j<ncat;j++) {
						qt[k][l] += q[j][l]*a[k+1][i+1]*a[j+1][i+1]/con;
					}
				}
            }
		for (k=i;k<ncat;k++) 
			for (l=i;l<degree;l++) {
				q[k][l] -= qt[k][l];
			}
	}

	fprintf(outf, "q matrix \n ");
	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, "%10.8lf  ", q[i][j]);
		}
		fprintf(outf, "\n ");
	}
	fprintf(outf, "\n ");

	fprintf(outf, "r matrix \n ");
	for(i=0; i<degree; i++) {
		for(j=0; j<degree; j++) {
			fprintf(outf, "%10.8lf ", r[i][j]);
		}
		fprintf(outf, "\n ");
	}
	fprintf(outf, "\n ");

	for(i=0; i<ncat; i++) {
		for(j=0; j<degree; j++) {
			a[i+1][j+1] = 0;
			for (k=0; k<degree; k++) {
				a[i+1][j+1] += q[i][k] * r[k][j];
			}
			if (fabs((a[i+1][j+1]-b[i+1][j+1])/a[i+1][j+1]) > 0.000001)
				fprintf(outf, "not match at %d %d \n", i, j );			
		}
	}

	for(i=0; i<degree; i++) {
		for(j=i; j<degree; j++) {
			sum = 0;
			for (k=0; k<ncat; k++) {
				sum += q[k][i] * q[k][j];
			}
			if (i!=j && fabs(sum) > 0.000000001)
				fprintf(outf, "not orthogonal for %d %d \n", i, j );			
		}
	}

	/* compute the Cr matrix */
	for(i=0; i<ncat; i++) 	
		for(j=0; j<degree; j++) 
			Cr[i][j] += (1.0/sqrt((double)np))*D[i]*q[i][j];

	/* output Cr matrix */
	fprintf(outf, "Cr matrix \n ");
	for(i=0; i<ncat; i++) 	{
		for(j=0; j<degree; j++) {
				fprintf(outf, " %13.10e ", Cr[i][j]);
		}
		fprintf(outf, "\n ");

	}

	fclose(outf);
	free_dvector(D, 0, ncat-1);
	free_dmatrix(A, 0, ncat-1, 0, ncat-1);
	free_dmatrix(a, 1, ncat, 1, degree);
	free_dmatrix(qt, 0, ncat-1, 0, degree-1);
	free_dmatrix(q, 0, ncat-1, 0, degree-1);
	free_dmatrix(r, 0, degree-1, 0, degree-1);
	free_dvector(c, 1, ncat);
	free_dvector(d, 1, ncat);

}
 
/*--------------------------------------------------------------------------
	PartialFPartialr
	
	functionality: 

	Computes the partial derivative of F to r in von Davier, Holland, & Thayer 
	(2004, Equation 5.21). The Kernel Method of Test Equating.
	call the numerical recipe function qrdcmp to do QR decomposition.
	
	author: Tianyou Wang 3/15/2005.
	
	input:
		ncat    Number of discrete score categories
		degree      The highest polynomial degree of the log-linear model
		np          sample size
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
		hx          kernel continuizing parameter
		x           X score values at which the partial derivatives are evaluated.
		            x can be non-integer values.
 
	output:
 		Fr          vector containing partial derivatives of F with respect to r
		            in von Davier et. al. Equation 5.12
--------------------------------------------------------------------------*/
void PartialFPartialr(int ncat, double *scores, double *fd, double hx, double *Fr, double x)
{
	int j;
	double Rjx, zjx;              /*the standardized normal random variable */
	double mu=0.0, sigma=0.0;         /*the mean and standard deviation of the
							                        unsmoothed distribution */
	double ax, sumsq=0.0;                           /*ax, and second moment */
	double temp2, Mjx, pdf;                            /*temporary variable */

	for (j=0; j<ncat; j++) {
		mu += scores[j] * fd[j];
		sumsq += scores[j] * scores[j] * fd[j];
	}

	sigma = sumsq - mu * mu;
	ax = sqrt(sigma / (sigma + hx * hx));
    
	for (j=0; j<ncat; j++) {
		Rjx = (x - ax * scores[j] - (1 - ax) * mu) / ax / hx;
		temp2 = StdNormalCdf(Rjx);
		zjx = (scores[j] - mu) / sqrt(sigma);
		Mjx = .5 * (x - mu) * (1 - ax * ax) * zjx * zjx + (1 - ax) * scores[j];
		pdf = KernelContinuPdf(ncat, scores, fd, hx, x);
		Fr[j] =  temp2 - Mjx * pdf;
	}
}

/*--------------------------------------------------------------------------
	FrCrSqNorm
	
	functionality: 

	Computes the sqaured norm of Fr multiplied by C matrix in von Davier, Holland, & Thayer 
	(2004, Equation 7.5). The Kernel Method of Test Equating.
	call the numerical recipe function qrdcmp to do QR decomposition.
	
	author: Tianyou Wang 3/15/2005.
	
	input:
		ncat    Number of discrete score categories
		degree      The highest polynomial degree of the log-linear model
		np          sample size
		scores      vector containing the discrete scores
		fd          vector containing the frequency distribution
 		Fr          vector containing partial derivatives of F with respect to r
 
	output:
 		return the square of the norm in von Davier, Holland, & Thayer 
	    (2004, Equation 7.5)
--------------------------------------------------------------------------*/
double FrCrSqNorm(int ncat, int degree, double *Fr, double **Cr)
{
	int i, j;
	double temp, norm=0;

	for(j=0; j<degree; j++) {
		temp = 0;
		for(i=0; i<ncat; i++) 
			temp += Fr[i] * Cr[i][j];
		norm += temp * temp;
	}

	return norm;

}

/*--------------------------------------------------------------------------
	vPMN
	
	functionality:

	Computes vectorized P and M and N matrix from the bivariate distribution 
	defined in von Davier, 	Holland, & Thayer (2004, Equations 2.8, 2.9, 2.10). 
	The Kernel Method of Test Equating.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		bdist       bivariate fitted distribution
		ncatx       Number of discrete score categories for the new form
		ncaty       Number of discrete score categories for the old form

 
	output:
	    vP          vectorized P 
 		M           M matrix
		N           N matrix
--------------------------------------------------------------------------*/

void vPMN(int ncatx, int ncaty, double **bdist, double *vP, double **M, double **N)
{
	int i, j, ncat;

	ncat = ncatx * ncaty;

	for (i=0; i<ncaty; i++) 
		for (j=0; j<ncatx; j++) 
			vP[i*ncatx+j] = bdist[j][i];

	for (i=0; i<ncatx; i++) 
		for (j=0; j<ncat; j++) M[i][j] = 0;

	for (i=0; i<ncatx; i++) 
		for (j=0; j<ncaty; j++) M[i][j*ncatx+i] = 1;

	for (i=0; i<ncaty; i++) 
		for (j=0; j<ncat; j++) N[i][j] = 0;

	for (i=0; i<ncaty; i++) 
		for (j=0; j<ncatx; j++) N[i][i*ncatx+j] = 1;

}

/*--------------------------------------------------------------------------
	vPP
	
	functionality:

	Computes vectorized P and M and N matrix from the bivariate distribution 
	defined in von Davier, 	Holland, & Thayer (2004, Equations 2.8, 2.9, 2.10). 
	The Kernel Method of Test Equating.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		bdist       bivariate fitted distribution
		ncatx       Number of discrete score categories for the new form
		ncaty       Number of discrete score categories for the old form

 
	output:
	    vPP         vectorized P 
--------------------------------------------------------------------------*/

void vPT(int ncatx, int ncaty, double **bdist, double *vPP)
{
	int i, j, ncat;

	ncat = ncatx * ncaty;

	for (i=0; i<ncatx; i++) 
		for (j=0; j<ncaty; j++) 
			vPP[i*ncaty+j] = bdist[i][j];


}

/*--------------------------------------------------------------------------
	MatrixMultiVector
	
	functionality:

	Computes a nrow x ncol matrix multipled by a ncol vector
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           ncol vector

 
	output:
	    r           resulting vector of nrow elements
--------------------------------------------------------------------------*/
void MatrixMultiVector(int nrow, int ncol, double **m, double *v, double *r)
{
	int i, j;

	for (i=0; i<nrow; i++) {
		r[i] = 0;
		for (j=0; j<ncol; j++) 
			r[i] += m[i][j] * v[j];
	}
}

/*--------------------------------------------------------------------------
	VectorMultiMatrix
	
	functionality:

	Computes a nrow vector  multipled by a nrow x ncol matrix
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow        number of rows
		ncol        number of columns
		m           nrow x ncol matrix
		v           nrow vector

 
	output:
	    r           resulting vector of ncol elements
--------------------------------------------------------------------------*/
void VectorMultiMatrix(int nrow, int ncol, double *v, double **m, double *r)
{
	int i, j;

	for (i=0; i<ncol; i++) {
		r[i] = 0;
		for (j=0; j<nrow; j++) r[i] += v[j] * m[j][i];
	}
}

/*--------------------------------------------------------------------------
	MatrixMultiMatrix
	
	functionality:

	Computes m (nrow1 x ncol1row2 matrix) multipled by n ( ncol1row2 x ncol2 matrix)
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		nrow1       number of rows for m
		ncol1row2   number of columns for m and number of rows for n
		ncol2       number of columnss for n
		m           nrow1 x ncol matrix
		n           ncol x ncol2 matrix

 
	output:
	    r           resulting matrix of nrow1 x ncol2
--------------------------------------------------------------------------*/
void MatrixMultiMatrix(int nrow1, int ncol1row2, int ncol2, double **m, double **n, double **r)
{
	int i, j, k;

	for (i=0; i<nrow1; i++) 
		for (j=0; j<ncol2; j++) r[i][j] = 0;
	
	for (i=0; i<nrow1; i++) 
		for (j=0; j<ncol2; j++) 
			for (k=0; k<ncol1row2; k++)  
				r[i][j] += m[i][k] * n[k][j];
	
}

/*--------------------------------------------------------------------------
	VectorNormSq
	
	functionality:

	Computes the square of the norm of a vector.
	
	author: Tianyou Wang 1/16/2005.
	
	input:
		ncat       Number of discrete score categories 
		v          vector with ncat elements
 
 
	output:
 		return the norm of the vector
--------------------------------------------------------------------------*/
double VectorNormSq(int ncat, double *v)
 {
	 int i;
	 double sum=0;

	 for (i=0; i<ncat; i++)	
		 sum += v[i] * v[i];

	 return sum;
 }


/*-----------------------------------------------------------------------------------------
	KernelEquateSEERG
	
	functionality: 

	Computes the standard error of equating (SEE) for the Random Groups Design
	(called EG design in von Davier, Holland, & Thayer 
	(2004, Equation 7.5)). The Kernel Method of Test Equating.
	call the numerical recipe function qrdcmp to do QR decomposition.
	
	author: Tianyou Wang 3/15/2005.
	
	input:
		ncatx        Number of discrete score categories for the new form
		degreex      The highest polynomial degree of the log-linear model for the new form
		npx          sample size for the new form
		scoresx      vector containing the discrete scores for the new form
		fdx          vector containing the relative frequency distribution for the new form  
		hx           kernel continuizing parameter for the new form
		Equatedx     the equated X scores
		ncaty        Number of discrete score categories for the old form
		degreey      The highest polynomial degree of the log-linear model for the old form
		npy          sample size for the old form
		scoresy      vector containing the discrete scores for the old form
		fdy          vector containing the relative frequency distribution for the old form  
		hy           kernel continuizing parameter for the old form
 
	output:
	    Equatedx     vector containing equated scores
 		SEE          vector containing the standard error of equating 
-----------------------------------------------------------------------------------------*/

void KernelEquateSEERG(int ncatx, int degreex, long npx, double *scoresx, double *fdx, 
					   int ncaty, int degreey, long npy, double *scoresy, 
					   double *fdy, double *Equatedx, double *SEE) 
{
	int i, j, nparax, nparay;
	double **Crx, **Cry, **Bx, **By;
	double *Fr, *Gs, Gp;
	double SqNormx, SqNormy;
	double hx, hy;
	

	nparax = degreex ;
	nparay = degreey ;
	Crx = dmatrix(0, ncatx-1, 0, ncatx-1);
	Cry = dmatrix(0, ncaty-1, 0, ncaty-1);
	Bx = dmatrix(0, nparax-1, 0, ncatx-1);
	By = dmatrix(0, nparay-1, 0, ncaty-1);
	Fr = dvector(0, ncatx-1);
	Gs = dvector(0, ncaty-1);

	for (i=0; i<nparax; i++)	{
		for (j=0; j<ncatx; j++)	{
			Bx[i][j] = 2 * pow(scoresx[j], (i+1));
		}
	}

	for (i=0; i<nparay; i++)	{
		for (j=0; j<ncaty; j++)	{
			By[i][j] = 2 * pow(scoresy[j], (i+1));
		}
	}

	ComputeCmatrixGen(ncatx, nparax, npx, Bx, fdx, Crx);
	ComputeCmatrixGen(ncaty, nparay, npy, By, fdy, Cry);

	hx = Optimalh(ncatx, scoresx, fdx, 1); 
	hy = Optimalh(ncaty, scoresy, fdy, 1); 

	KernelEquate(ncatx, scoresx, fdx, hx, ncaty, scoresy, fdy, hy, 
		Equatedx); 


	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatx, scoresx, fdx, hx, Fr, scoresx[i]);
        PartialFPartialr(ncaty, scoresy, fdy, hy, Gs, Equatedx[i]);
        Gp = KernelContinuPdf(ncaty, scoresy, fdy, hy, Equatedx[i]);
		SqNormx = FrCrSqNorm(ncatx, nparax, Fr, Crx); 
		SqNormy = FrCrSqNorm(ncaty, nparay, Gs, Cry);
		SEE[i] = sqrt(SqNormx + SqNormy) / Gp;
	}
 
	free_dmatrix(Crx, 0, ncatx-1, 0, ncatx-1);
	free_dmatrix(Cry, 0, ncaty-1, 0, ncaty-1);
	free_dmatrix(Bx, 0, ncatx-1, 0, ncatx-1);
	free_dmatrix(By, 0, ncaty-1, 0, ncaty-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);

}

/*--------------------------------------------------------------------------
	KernelEquateSEESG
	
	functionality:

	Computes kernel equating function and SEE for the Single Group (SG) design
	based on continuized cdf in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 3/26/2005.
	
	input:
		bivar       information about bivariate distribution of X and Y
 
	output:
 		Equatedx    a vector containing the equated score 
		SEE         a vector containing standard error of equating 
--------------------------------------------------------------------------*/
void KernelEquateSEESG(struct BLL_SMOOTH *bivar, double *Equatedx, double *SEE)
{
	int i, j, k, ncat, ncatx, ncaty, npara, cu, cv, cuv;
	long np;
	double **fitbdist;
	double *vP, *vPP, **M, **N, *r, *s, **U, **V;
	double hx, hy;
	double **Cr, **B;
	double *Fr, *Gs, Gp, *FrU, *GsV;
	double *scoresx, *scoresy;
	int *interx, *intery;
	int **cpm;

	cuv = bivar->cuv;
	cu = bivar->cu;
	cv = bivar->cv;
	cpm = imatrix(0, 1, 0, cuv);
    for(i=0;i<cuv;i++) for(j=0;j<2;j++) cpm[j][i] = bivar->cpm[i][j];
	
	ncatx = bivar->nsx;
	ncaty = bivar->nsv;
	npara = bivar->cu + bivar->cv + bivar->cuv;
	interx = cpm[0];
	intery = cpm[1];
	np = bivar->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar->minx + i * bivar->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar->minv + i * bivar->incv;
	ncat = ncatx * ncaty;

	Cr = dmatrix(0, ncat-1, 0, npara-1);
	B = dmatrix(0, npara-1, 0, ncat-1); /*be aware B transpose(!!) as same dimension as Cr.  */
	U = dmatrix(0, ncatx-1, 0, npara-1);
	V = dmatrix(0, ncaty-1, 0, npara-1);
	Fr = dvector(0, ncatx-1);
	Gs = dvector(0, ncaty-1);
	FrU = dvector(0, npara-1);
	GsV = dvector(0, npara-1);
 	fitbdist = dmatrix(0, ncatx-1, 0, ncaty-1);

	for(i = 0; i <ncatx; i++) 
		for(j = 0; j <ncaty; j++) fitbdist[i][j] = bivar->bfd[i][j] / np;

	for (i=0; i<npara; i++)
		for (j=0; j<ncat; j++) B[i][j] = 0;

    /* The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. The B matrix is a transpose of the design matrix in loglinear model
	   (without the first column)
       First, assign natrual design matrix first for rows corresponding to x */	

	for (i=0; i<cu; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatx)+j] = pow((double) j, (double) (i+1));
			}
		}
	}
    /* then assign values to rows corresponding to y */
	for (i=cu; i<cu+cv; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatx)+j] = pow((double) k,(double) (i-cu+1));
			}
		}
	}
/* assign value to the last rows corresponding to the interaction terms */
	for (i=0; i<cuv; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i+cu+cv][k*(ncatx)+j] = pow((double) j, (double) interx[i]) * 
                                          pow((double) k, (double) intery[i]);
			}
		}
	}

	vP = dvector(0, ncat-1);
	vPP = dvector(0, ncat-1);
	M = dmatrix(0, ncatx-1, 0, ncat-1);
	N = dmatrix(0, ncaty-1, 0, ncat-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);

	vPMN(ncatx, ncaty, fitbdist, vP, M, N);
	vPT(ncatx, ncaty, fitbdist, vPP);
	MatrixMultiVector(ncatx, ncat, M, vP, r);
	MatrixMultiVector(ncaty, ncat, N, vP, s);

	ComputeCmatrixGen(ncat, npara, np, B, vP, Cr);
	MatrixMultiMatrix(ncatx, ncat, npara, M, Cr, U);
	MatrixMultiMatrix(ncaty, ncat, npara, N, Cr, V);

	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 

	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, 
		Equatedx); 

	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatx, scoresx, r, hx, Fr, scoresx[i]);
        PartialFPartialr(ncaty, scoresy, s, hy, Gs, Equatedx[i]);
        Gp = KernelContinuPdf(ncaty, scoresy, s, hy, Equatedx[i]);
		VectorMultiMatrix(ncatx, npara, Fr, U, FrU);
		VectorMultiMatrix(ncaty, npara, Gs, V, GsV);
		for (j=0; j<npara; j++) FrU[j] -= GsV[j];
		SEE[i] = sqrt(VectorNormSq(npara, FrU)) / Gp;
	}
 
	free_dvector(vP, 0, ncat-1);
	free_dvector(vPP, 0, ncat-1);
	free_dmatrix(M, 0, ncatx-1, 0, ncat-1);
	free_dmatrix(N, 0, ncaty-1, 0, ncat-1);
	free_dmatrix(Cr, 0, ncat-1, 0, npara-1);
	free_dmatrix(B, 0, ncat-1, 0, ncat-1);
	free_dmatrix(U, 0, ncatx-1, 0, npara-1);
	free_dmatrix(V, 0, ncaty-1, 0, npara-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);
	free_dvector(FrU, 0, npara-1);
	free_dvector(GsV, 0, npara-1);
 	free_dmatrix(fitbdist, 0, ncatx-1, 0, ncaty-1);

}

/*--------------------------------------------------------------------------
	KernelEquateSG
	
	functionality:

	Computes kernel equating function and SEE for the Single Group (SG) design
	based on continuized cdf in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 3/26/2005.
	
	input:
		bivar       information about bivariate distribution of X and Y
 
	output:
 		Equatedx    a vector containing the equated score 
--------------------------------------------------------------------------*/
void KernelEquateSG(struct BLL_SMOOTH *bivar, double *Equatedx)
{
	int i, j, ncat, ncatx, ncaty;
	long np;
	double **fitbdist;
	double *vP, **M, **N, *r, *s;
	double hx, hy;
	double *scoresx, *scoresy;
	
	ncatx = bivar->nsx;
	ncaty = bivar->nsv;
	np = bivar->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar->minx + i * bivar->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar->minv + i * bivar->incv;
	ncat = ncatx * ncaty;

 	fitbdist = dmatrix(0, ncatx-1, 0, ncaty-1);

	for(i = 0; i <ncatx; i++) 
		for(j = 0; j <ncaty; j++) fitbdist[i][j] = bivar->bfd[i][j] / np;


	vP = dvector(0, ncat-1);
	M = dmatrix(0, ncatx-1, 0, ncat-1);
	N = dmatrix(0, ncaty-1, 0, ncat-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);

	vPMN(ncatx, ncaty, fitbdist, vP, M, N);
	MatrixMultiVector(ncatx, ncat, M, vP, r);
	MatrixMultiVector(ncaty, ncat, N, vP, s);


	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 

	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, 
		Equatedx); 

 	free_dmatrix(fitbdist, 0, ncatx-1, 0, ncaty-1);
	free_dvector(vP, 0, ncat-1);
	free_dmatrix(M, 0, ncatx-1, 0, ncat-1);
	free_dmatrix(N, 0, ncaty-1, 0, ncat-1);
	free_dvector(r, 0, ncatx-1);
	free_dvector(s, 0, ncaty-1);

}

/*--------------------------------------------------------------------------
	KernelEquateSEECB
	
	functionality:

	Computes kernel equating function and SEE for the counter balance (CB) 
	design based on continuized cdf in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 3/27/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         vector containing weights for the form taken first.
		            wts[0] is for form X, wts[1] is for form Y

 
	output:
 		Equatedx    a vector containing the equated score 
		SEE         a vector containing standard error of equating 
--------------------------------------------------------------------------*/
void KernelEquateSEECB(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *wts, double *Equatedx, 
					double *SEE)
{
	int i, j, k, ncat, ncatx, ncaty, npara, cu, cv, cuv;
	long np1, np2;
	double *vP12, *vP21, **M, **N, *r1, *s1, *r2, *s2, *r, *s;
	/* r1 is for X taken first(first group), r2 is for X taken second (second group),
	   s1 is for Y taken first(second group), s2 is for Y taken second (first group)
	   vP12 is for first group, vP21 is for second group */
	double **U12, **U21, **V12, **V21;
	double hx, hy;
	double **Cr12, **Cr21, **B;
	double *Fr, *Gs, Gp, *FrU12, *GsV12, *FrU21, *GsV21;
	double *scoresx, *scoresy;
	double **fitbdist1, **fitbdist2;
	int *interx, *intery;
	int **cpm1;

	cuv = bivar1->cuv;
	cu = bivar1->cu;
	cv = bivar1->cv;
	cpm1 = imatrix(0, 1, 0, cuv);
    for(i=0;i<cuv;i++) for(j=0;j<2;j++) cpm1[j][i] = bivar1->cpm[i][j];

	ncatx = bivar1->nsx;
	ncaty = bivar1->nsv;
	npara = bivar1->cu + bivar1->cv + bivar1->cuv ;
	interx = cpm1[0];
	intery = cpm1[1];
	np1 = bivar1->num_persons;
	np2 = bivar2->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar1->minx + i * bivar1->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar1->minv + i * bivar1->incv;
	ncat = ncatx * ncaty;

	vP12 = dvector(0, ncat-1);
	vP21 = dvector(0, ncat-1);
	M = dmatrix(0, ncatx-1, 0, ncat-1);
	N = dmatrix(0, ncaty-1, 0, ncat-1);
	r1 = dvector(0, ncatx-1);
	s1 = dvector(0, ncaty-1);
	r2 = dvector(0, ncatx-1);
	s2 = dvector(0, ncaty-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);
	Cr12 = dmatrix(0, ncat-1, 0, npara-1);
	Cr21 = dmatrix(0, ncat-1, 0, npara-1);
	B = dmatrix(0, npara-1, 0, ncat-1);
	U12 = dmatrix(0, ncatx-1, 0, npara-1);
	V12 = dmatrix(0, ncaty-1, 0, npara-1);
	U21 = dmatrix(0, ncatx-1, 0, npara-1);
	V21 = dmatrix(0, ncaty-1, 0, npara-1);
	Fr = dvector(0, ncatx-1);
	Gs = dvector(0, ncaty-1);
	FrU12 = dvector(0, npara-1);
	GsV12 = dvector(0, npara-1);
	FrU21 = dvector(0, npara-1);
	GsV21 = dvector(0, npara-1);
	fitbdist1 = dmatrix(0, ncatx-1, 0, ncaty-1);
	fitbdist2 = dmatrix(0, ncatx-1, 0, ncaty-1);



    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. note that this B is the transpose of the design matrix in the loglinear 
	   model. */

	for (i=0; i<npara; i++)
		for (j=0; j<ncat; j++) B[i][j] = 0;

	/*First, assign natrual design matrix first for rows corresponding to x */	
	for (i=0; i<cu; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatx)+j] = pow((double)j,(double)(i+1));
			}
		}
	}
    /* then assign values to rows corresponding to y */
	for (i=cu; i<cu+cv; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatx)+j] = pow((double)k,(double)(i-cu+1));
			}
		}
	}
/* assign value to the last rows corresponding to the interaction terms */
	for (i=0; i<cuv; i++)
	{
		for (j=0; j<ncatx; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i+cu+cv][k*(ncatx)+j] = pow((double)j, (double)interx[i]) *
                                          pow((double)k, (double)intery[i]);
			}
		}
	}

	for(i = 0; i <ncatx; i++) {
		for(j = 0; j <ncaty; j++) {
			fitbdist1[i][j] = bivar1->bfd[i][j] / np1;
			fitbdist2[i][j] = bivar2->bfd[i][j] / np2;
		}
	}

	vPMN(ncatx, ncaty, fitbdist1, vP12, M, N);
	vPMN(ncatx, ncaty, fitbdist2, vP21, M, N);
	MatrixMultiVector(ncatx, ncat, M, vP12, r1);
	MatrixMultiVector(ncaty, ncat, N, vP12, s2);
	MatrixMultiVector(ncatx, ncat, M, vP21, r2);
	MatrixMultiVector(ncaty, ncat, N, vP21, s1);

	for (i=0; i<ncatx; i++) 
		r[i] = wts[0] * r1[i] + (1 - wts[0]) * r2[i];
	for (i=0; i<ncaty; i++) 
		s[i] = wts[1] * s1[i] + (1 - wts[1]) * s2[i];


	ComputeCmatrixGen(ncat, npara, np1, B, vP12, Cr12);
	ComputeCmatrixGen(ncat, npara, np2, B, vP21, Cr21);
	MatrixMultiMatrix(ncatx, ncat, npara, M, Cr12, U12);
	MatrixMultiMatrix(ncaty, ncat, npara, N, Cr12, V12);
	MatrixMultiMatrix(ncatx, ncat, npara, M, Cr21, U21);
	MatrixMultiMatrix(ncaty, ncat, npara, N, Cr21, V21);

	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 

 	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, 
		Equatedx); 

	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatx, scoresx, r, hx, Fr, scoresx[i]);
        PartialFPartialr(ncaty, scoresy, s, hy, Gs, Equatedx[i]);
        Gp = KernelContinuPdf(ncaty, scoresy, s, hy, Equatedx[i]);
		VectorMultiMatrix(ncatx, npara, Fr, U12, FrU12);
		VectorMultiMatrix(ncaty, npara, Gs, V12, GsV12);
		VectorMultiMatrix(ncatx, npara, Fr, U21, FrU21);
		VectorMultiMatrix(ncaty, npara, Gs, V21, GsV21);
		for (j=0; j<npara; j++) {
			FrU12[j] *= wts[0];
			FrU12[j] -= (1 - wts[1]) * GsV12[j];
			FrU21[j] *= (1 - wts[0]);
			FrU21[j] -= wts[1] * GsV21[j];
		}
		SEE[i] = sqrt(VectorNormSq(npara, FrU12) + VectorNormSq(npara, FrU21)) / Gp;
	}
    	
	free_dvector(vP12, 0, ncat-1);
	free_dvector(vP21, 0, ncat-1);
	free_dmatrix(M, 0, ncatx-1, 0, ncat-1);
	free_dmatrix(N, 0, ncaty-1, 0, ncat-1);
	free_dvector(r1, 0, ncatx-1);
	free_dvector(s1, 0, ncaty-1);
	free_dvector(r2, 0, ncatx-1);
	free_dvector(s2, 0, ncaty-1);
	free_dvector(r, 0, ncatx-1);
	free_dvector(s, 0, ncaty-1);
	free_dmatrix(Cr12, 0, ncat-1, 0, npara-1);
	free_dmatrix(Cr21, 0, ncat-1, 0, npara-1);
	free_dmatrix(B, 0, ncat-1, 0, ncat-1);
	free_dmatrix(U12, 0, ncatx-1, 0, npara-1);
	free_dmatrix(V12, 0, ncaty-1, 0, npara-1);
	free_dmatrix(U21, 0, ncatx-1, 0, npara-1);
	free_dmatrix(V21, 0, ncaty-1, 0, npara-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);
	free_dvector(FrU12, 0, npara-1);
	free_dvector(GsV12, 0, npara-1);
	free_dvector(FrU21, 0, npara-1);
	free_dvector(GsV21, 0, npara-1);
	free_dmatrix(fitbdist1, 0, ncatx-1, 0, ncaty-1);
	free_dmatrix(fitbdist2, 0, ncatx-1, 0, ncaty-1);
}

/*--------------------------------------------------------------------------
	KernelEquateSEENEATPS
	
	functionality:

	Computes kernel equating function and SEE for the NEAT design 
	with post stratification method based on procedures in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	Author: Tianyou Wang 3/29/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         wts for group taking form X.
 
	output:
 		Equatedx    a vector containing the equated score 
		SEE         a vector containing standard error of equating 
--------------------------------------------------------------------------*/
void KernelEquateSEENEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx, double *SEE)
{
	int i, j, k, ncat1, ncat2, ncatx, ncaty, ncatv, npara1, npara2;
	int cu1, cv1, cuv1, cu2, cv2, cuv2;
	long np1, np2;
	double *vP1, *vP2, *vPP1, *vPP2, **M1, **N1, **M2, **N2, *r1, *s1, *r2, *s2, *r, *s,  *v1, *v2;
	/* r1 is for X taken first(first group), r2 is for X taken second (second group),
	   s1 is for Y taken first(second group), s2 is for Y taken second (first group)
	   vP1 is for first group, vP2 is for second group */
	double *vP121, *vP212;
	double **UR, **US, **VR, **VS, **UP, **UPStar, **UQ, **UQStar, **Utemp, **Vtemp;
	double *pl, *vPl, *ql, *vQl;
	double hx, hy;
	double **Cr1, **Cr1Star, **Cr2, **Cr2Star, **B1, **B2;
	double *Fr, *Gs, Gp, *FrUR, *GsVR, *FrUS, *GsVS;
	double *scoresx, *scoresy;
	int *interx, *interv1, *interv2, *intery;
	double **fitbdist1, **fitbdist2, sum1=0, sum2=0;
	int currentrow1, currentrow2;
	int **cpm1, **cpm2;

	cuv1 = bivar1->cuv;
	cu1 = bivar1->cu;
	cv1 = bivar1->cv;
	cuv2 = bivar2->cuv;
	cu2 = bivar2->cu;
	cv2 = bivar2->cv;
	cpm1 = imatrix(0, 1, 0, cuv1);
    for(i=0;i<cuv1;i++) for(j=0;j<2;j++) cpm1[j][i] = bivar1->cpm[i][j];
	cpm2 = imatrix(0, 1, 0, cuv2);
    for(i=0;i<cuv2;i++) for(j=0;j<2;j++) cpm2[j][i] = bivar2->cpm[i][j];


	ncatx = bivar1->nsx;
	ncaty = bivar2->nsx;
	ncatv = bivar1->nsv;
	npara1 = bivar1->cu + bivar1->cv + bivar1->cuv + 10 ;
	npara2 = bivar2->cu + bivar2->cv + bivar2->cuv + 10 ;
	interx = cpm1[0];
	interv1 = cpm1[1];
	intery = cpm2[0];
	interv2 = cpm2[1];
	np1 = bivar1->num_persons;
	np2 = bivar2->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar1->minx + i * bivar1->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar2->minx + i * bivar2->incx;
	ncat1 = ncatx * ncatv;
	ncat2 = ncaty * ncatv;

	vP1 = dvector(0, ncat1-1);
	vP2 = dvector(0, ncat2-1);
	vPP1 = dvector(0, ncat1-1);
	vPP2 = dvector(0, ncat2-1);
	vP121 = dvector(0, ncat1-1);
	vP212 = dvector(0, ncat2-1);
	M1 = dmatrix(0, ncatx-1, 0, ncat1-1);
	N1 = dmatrix(0, ncatv-1, 0, ncat1-1);
	M2 = dmatrix(0, ncaty-1, 0, ncat2-1);
	N2 = dmatrix(0, ncatv-1, 0, ncat2-1);
	r1 = dvector(0, ncatx-1);
	s1 = dvector(0, ncaty-1);
	r2 = dvector(0, ncatx-1);
	s2 = dvector(0, ncaty-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);
	v1 = dvector(0, ncatv-1);
	v2 = dvector(0, ncatv-1);
	Cr1 = dmatrix(0, ncat1-1, 0, npara1-1);
	Cr2 = dmatrix(0, ncat2-1, 0, npara2-1);
	Cr1Star = dmatrix(0, ncat1-1, 0, npara1-1);
	Cr2Star = dmatrix(0, ncat2-1, 0, npara2-1);
	B1 = dmatrix(0, npara1-1, 0, ncat1-1);
	B2 = dmatrix(0, npara2-1, 0, ncat2-1);
	UR = dmatrix(0, ncatx-1, 0, npara1-1);
	UP = dmatrix(0, ncatx-1, 0, npara1-1);
	UPStar = dmatrix(0, ncatx-1, 0, npara1-1);
	Utemp = dmatrix(0, ncatx-1, 0, npara1-1);
	VR = dmatrix(0, ncaty-1, 0, npara1-1);
	US = dmatrix(0, ncatx-1, 0, npara2-1);
	VS = dmatrix(0, ncaty-1, 0, npara2-1);
	UQ = dmatrix(0, ncaty-1, 0, npara2-1);
	UQStar = dmatrix(0, ncaty-1, 0, npara2-1);
	Vtemp = dmatrix(0, ncaty-1, 0, npara2-1);
	pl = dvector(0, ncatx-1);
	vPl = dvector(0, npara1-1);
	ql = dvector(0, ncaty-1);
	vQl = dvector(0, npara2-1);
	Fr = dvector(0, ncatx-1);
	Gs = dvector(0, ncaty-1);
	FrUR = dvector(0, npara1-1);
	GsVR = dvector(0, npara1-1);
	FrUS = dvector(0, npara2-1);
	GsVS = dvector(0, npara2-1);
    fitbdist1 = dmatrix(0, ncatx-1, 0, ncatv-1);
    fitbdist2 = dmatrix(0, ncaty-1, 0, ncatv-1);

    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. note that this B is the transpose of the design matrix in the loglinear 
	   model. */

	/*First, assign natrual design matrix first for rows corresponding to x */	
	for (i=0; i<cu1; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[i][k*(ncatx)+j] = pow((double)j,(double)(i+1));
			}
		}
	}

    /* then assign values to rows corresponding to v */
	currentrow1 = cu1 ;
	for (i=0; i<cv1; i++)	{
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[currentrow1+i][k*(ncatx)+j] = pow((double)k,(double)(i+1));
			}
		}
	}
/* assign value to the last rows corresponding to the interaction terms */
	currentrow1 += cv1;
	for (i=0; i<cuv1; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[currentrow1+i][k*(ncatx)+j] = pow((double)j, (double)interx[i]) * 
                                                 pow((double)k, (double)interv1[i]);
			}
		}
	}

/* assign values to the rows for the 0 score on form X */
	currentrow1 += cuv1 ;
	for (j=0; j<ncatx; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (j==0) 
				B1[currentrow1][k*(ncatx)+j] = 1;
			else
				B1[currentrow1][k*(ncatx)+j] = 0;
		}
	}

/* assign values to the rows for the 0 score on common set v */
	currentrow1++;
	for (j=0; j<ncatx; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (k==0) 
				B1[currentrow1][k*(ncatx)+j] = 1;
			else
				B1[currentrow1][k*(ncatx)+j] = 0;
		}
	}

	/*assign value to the rows for score gaps in x 	*/
	currentrow1++;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				if ((j%5)==0 && j!=0) 
					B1[currentrow1+i][k*(ncatx)+j] = pow((double)j,(double)i);
				else
					B1[currentrow1+i][k*(ncatx)+j] = 0;
			}
		}
	}

	/*assign value to the rows for score gaps in v */
	currentrow1 += 4;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				if (((k)%5)==0 && k!=0) 
					B1[currentrow1+i][k*(ncatx)+j] = pow((double)k,(double)i);
				else
					B1[currentrow1+i][k*(ncatx)+j] = 0;
			}
		}
	}

    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. note that this B2 is the transpose of the design matrix in the loglinear 
	   model. */

	/*First, assign natrual design matrix first for rows corresponding to y */	
	for (i=0; i<cu2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[i][k*(ncaty)+j] = pow((double)j,(double)(i+1));
			}
		}
	}
    /* then assign values to rows corresponding to v */
	currentrow2 = cu2 ;
	for (i=0; i<cv2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[currentrow2+i][k*(ncaty)+j] = pow((double)k,(double)(i+1));
			}
		}
	}

/* assign value to the last rows corresponding to the interaction terms */
	currentrow2 += cv2;
	for (i=0; i<cuv2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[currentrow2+i][k*(ncaty)+j] = pow((double)j, (double)intery[i]) * 
                                                 pow((double)k, (double)interv2[i]);
			}
		}
	}


    /* assign values to the rows for the 0 score on form Y  */
	currentrow2 += cuv2 ;
	for (j=0; j<ncaty; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (j==0) 
				B2[currentrow2][k*(ncaty)+j] = 1;
			else
				B2[currentrow2][k*(ncaty)+j] = 0;
		}
	}

    /* assign values to the rows for the 0 score on common set V  */
	currentrow2++;
	for (j=0; j<ncaty; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (k==0) 
				B2[currentrow2][k*(ncaty)+j] = 1;
			else
				B2[currentrow2][k*(ncaty)+j] = 0;
		}
	}

	/*assign value to the rows for score gaps in y 	*/
	currentrow2++;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				if ((j%5)==0 && j!=0) 
					B2[currentrow2+i][k*(ncaty)+j] = pow((double)j,(double)i);
				else
					B2[currentrow2+i][k*(ncaty)+j] = 0;
			}
		}
	}

	/*assign value to the rows for score gaps in v */
	currentrow2 += 4;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				if (((k)%5)==0 && k!=0) 
					B2[currentrow2+i][k*(ncaty)+j] = pow((double)k,(double)i);
				else
					B2[currentrow2+i][k*(ncaty)+j] = 0;
			}
		}
	}

	sum1 = 0;
	sum2 = 0;
	for(j = 0; j <ncatv; j++) {
		for(i = 0; i <ncatx; i++) {
			fitbdist1[i][j] = bivar1->bfd[i][j] / np1;
			sum1 += fitbdist1[i][j];
		}
		for(i = 0; i <ncaty; i++) {
			fitbdist2[i][j] = bivar2->bfd[i][j] / np2;
			sum2 += fitbdist2[i][j];
		}
	}


	vPMN(ncatx, ncatv, fitbdist1, vP1, M1, N1);
	vPMN(ncaty, ncatv, fitbdist2, vP2, M2, N2);
	vPT(ncatx, ncatv, fitbdist1, vPP1);
	vPT(ncaty, ncatv, fitbdist2, vPP2);
	MatrixMultiVector(ncatx, ncat1, M1, vP1, r1);
	MatrixMultiVector(ncatv, ncat1, N1, vP1, v1);
	MatrixMultiVector(ncaty, ncat2, M2, vP2, s2);
	MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);

	ComputeCmatrixGen(ncat1, npara1, np1, B1, vP1, Cr1);
	ComputeCmatrixGen(ncat2, npara2, np2, B2, vP2, Cr2);


	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncatx; i++) {
			vP121[i+j*ncatx] = vP1[i+j*ncatx] * v2[j] / v1[j];
			for (k=0; k<npara1; k++) 
				Cr1Star[i+j*ncatx][k] = Cr1[i+j*ncatx][k] * v2[j] / v1[j];
		}
	}

	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncaty; i++) {
			vP212[i+j*ncaty] = vP2[i+j*ncaty] * v1[j] / v2[j];
			for (k=0; k<npara2; k++) 
				Cr2Star[i+j*ncatx][k] =  Cr2[i+j*ncatx][k] * v1[j] / v2[j];
		}
	}
	for (i=0; i<ncatx; i++) {
		for (k=0; k<npara1; k++) Utemp[i][k] = 0;
		for (k=0; k<npara2; k++) US[i][k] = 0;
	}

	for (i=0; i<ncaty; i++) {	
		for (k=0; k<npara2; k++) Vtemp[i][k] = 0;
		for (k=0; k<npara1; k++) VR[i][k] = 0;
	}

	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncatx; i++) pl[i] = fitbdist1[i][j];
		for (k=0; k<npara1; k++) {
			for (vPl[k]=0, i=0; i<ncatx; i++) 
				vPl[k] += Cr1[j*ncatx+i][k];
		}
		for (i=0; i<ncatx; i++) 
			for (k=0; k<npara1; k++) 
				Utemp[i][k] += v2[j] / v1[j] / v1[j] * pl[i] * vPl[k];

		for (i=0; i<ncaty; i++) ql[i] = fitbdist2[i][j];
		for (k=0; k<npara2; k++) {
			for (vQl[k]=0, i=0; i<ncaty; i++) 
				vQl[k] += Cr2[j*ncaty+i][k];
		}
		for (i=0; i<ncaty; i++) 
			for (k=0; k<npara2; k++) 
				Vtemp[i][k] += v1[j] / v2[j] /  v2[j]* ql[i] * vQl[k];

		for (i=0; i<ncaty; i++) 
			for (k=0; k<npara1; k++) VR[i][k] += wts / v2[j] * ql[i] * vPl[k];

		for (i=0; i<ncatx; i++) 
			for (k=0; k<npara2; k++) US[i][k] += (1 - wts) / v1[j] * pl[i] * vQl[k];
	}

	MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1, UP);
	MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2, UQ);
	MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1Star, UPStar);
	MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2Star, UQStar);
	MatrixMultiVector(ncatx, ncat1, M1, vP121, r2);
	MatrixMultiVector(ncatx, ncat2, M2, vP212, s1);


	for (i=0; i<ncatx; i++) 
		r[i] = wts * r1[i] + (1 - wts) * r2[i];
	for (i=0; i<ncaty; i++) 
		s[i] = wts * s1[i] + (1 - wts) * s2[i];

	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 
	
	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, Equatedx); 

	for (i=0; i<ncatx; i++) {
		for (k=0; k<npara1; k++) 
			UR[i][k] = wts * UP[i][k] + (1 - wts) * (UPStar[i][k] - Utemp[i][k]);
	}

	for (i=0; i<ncaty; i++) {
		for (k=0; k<npara2; k++) 
			VS[i][k] = (1 - wts) * UQ[i][k] + wts * (UQStar[i][k] - Vtemp[i][k]);
	}
  	
	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatx, scoresx, r, hx, Fr, scoresx[i]);
        PartialFPartialr(ncaty, scoresy, s, hy, Gs, Equatedx[i]);
        Gp = KernelContinuPdf(ncaty, scoresy, s, hy, Equatedx[i]);
		VectorMultiMatrix(ncatx, npara1, Fr, UR, FrUR);
		VectorMultiMatrix(ncaty, npara1, Gs, VR, GsVR);
		VectorMultiMatrix(ncatx, npara2, Fr, US, FrUS);
		VectorMultiMatrix(ncaty, npara2, Gs, VS, GsVS);
		for (j=0; j<npara1; j++) 
			FrUR[j] -= GsVR[j];
		for (j=0; j<npara2; j++) 
			FrUS[j] -= GsVS[j];
		SEE[i] = sqrt(VectorNormSq(npara1, FrUR) + VectorNormSq(npara2, FrUS)) / Gp;
	}
    	
	free_dvector(vP1, 0, ncat1-1);
	free_dvector(vP2, 0, ncat2-1);
	free_dvector(vP121, 0, ncat1-1);
	free_dvector(vP212, 0, ncat2-1);
	free_dmatrix(M1, 0, ncatx-1, 0, ncat1-1);
	free_dmatrix(N1, 0, ncatv-1, 0, ncat1-1);
	free_dmatrix(M2, 0, ncaty-1, 0, ncat2-1);
	free_dmatrix(N2, 0, ncatv-1, 0, ncat2-1);
	free_dvector(r1, 0, ncatx-1);
	free_dvector(s1, 0, ncaty-1);
	free_dvector(r2, 0, ncatx-1);
	free_dvector(s2, 0, ncaty-1);
	free_dvector(r, 0, ncatx-1);
	free_dvector(s, 0, ncaty-1);
	free_dvector(v1, 0, ncatv-1);
	free_dvector(v2, 0, ncatv-1);
	free_dmatrix(Cr1, 0, ncat1-1, 0, npara1-1);
	free_dmatrix(Cr2, 0, ncat2-1, 0, npara2-1);
	free_dmatrix(Cr1Star, 0, ncat1-1, 0, npara1-1);
	free_dmatrix(Cr2Star, 0, ncat2-1, 0, npara2-1);
	free_dmatrix(B1, 0, npara1-1, 0, ncat1-1);
	free_dmatrix(B2, 0, npara2-1, 0, ncat2-1);
	free_dmatrix(UR, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(UP, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(UPStar, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(Utemp, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(VR, 0, ncaty-1, 0, npara1-1);
	free_dmatrix(US, 0, ncatx-1, 0, npara2-1);
	free_dmatrix(VS, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(UQ, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(UQStar, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(Vtemp, 0, ncaty-1, 0, npara2-1);
	free_dvector(pl, 0, ncatx-1);
	free_dvector(vPl, 0, npara1-1);
	free_dvector(ql, 0, ncaty-1);
	free_dvector(vQl, 0, npara2-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);
	free_dvector(FrUR, 0, npara1-1);
	free_dvector(GsVR, 0, npara1-1);
	free_dvector(FrUS, 0, npara2-1);
	free_dvector(GsVS, 0, npara2-1);
	free_dmatrix(fitbdist1, 0, ncatx-1, 0, ncatv-1);
	free_dmatrix(fitbdist2, 0, ncaty-1, 0, ncatv-1);
}

/*--------------------------------------------------------------------------
	KernelEquateSEENEATPS3
	
	functionality:

	Computes kernel equating function and SEE for the NEAT design 
	with post stratification method based on procedures in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 3/29/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         wts for group taking form X.
 
	output:
 		Equatedx    a vector containing the equated score 
		SEE         a vector containing standard error of equating 
--------------------------------------------------------------------------*/
void KernelEquateSEENEATPS3(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx, double *SEE)
{
	int i, j, k, ncat1, ncat2, ncatx, ncaty, ncatv, npara1, npara2;
	int cu1, cv1, cuv1, cu2, cv2, cuv2;
	long np1, np2;
	double *vP1, *vP2, *vPP1, *vPP2, **M1, **N1, **M2, **N2, *r1, *s1, *r2, *s2, *r, *s,  *v1, *v2;
	/* r1 is for X taken first(first group), r2 is for X taken second (second group),
	   s1 is for Y taken first(second group), s2 is for Y taken second (first group)
	   vP1 is for first group, vP2 is for second group */
	double *vP121, *vP212;
	double **UR, **US, **VR, **VS, **UP, **UPStar, **UQ, **UQStar, **Utemp, **Vtemp;
	double *pl, *vPl, *ql, *vQl;
	double hx, hy;
	double **Cr1, **Cr1Star, **Cr2, **Cr2Star, **B1, **B2;
	double *Fr, *Gs, Gp, *FrUR, *GsVR, *FrUS, *GsVS;
	double *scoresx, *scoresy;
	int *interx, *interv1, *interv2, *intery;
	double **fitbdist1, **fitbdist2, sum1=0, sum2=0;
	int currentrow1, currentrow2;
	int **cpm1, **cpm2;

	cuv1 = bivar1->cuv;
	cu1 = bivar1->cu;
	cv1 = bivar1->cv;
	cuv2 = bivar2->cuv;
	cu2 = bivar2->cu;
	cv2 = bivar2->cv;
	cpm1 = imatrix(0, 1, 0, cuv1);
    for(i=0;i<cuv1;i++) for(j=0;j<2;j++) cpm1[j][i] = bivar1->cpm[i][j];
	cpm2 = imatrix(0, 1, 0, cuv2);
    for(i=0;i<cuv2;i++) for(j=0;j<2;j++) cpm2[j][i] = bivar2->cpm[i][j];

	ncatx = bivar1->nsx;
	ncaty = bivar2->nsx;
	ncatv = bivar1->nsv;
	npara1 = bivar1->cu + bivar1->cv + bivar1->cuv +  10;
	npara2 = bivar2->cu + bivar2->cv + bivar2->cuv +  10;
	interx = cpm1[0];
	interv1 = cpm1[1];
	intery = cpm2[0];
	interv2 = cpm2[1];
	np1 = bivar1->num_persons;
	np2 = bivar2->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar1->minx + i * bivar1->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar2->minx + i * bivar2->incx;
	ncat1 = ncatx * ncatv;
	ncat2 = ncaty * ncatv;

	vP1 = dvector(0, ncat1-1);
	vP2 = dvector(0, ncat2-1);
	vPP1 = dvector(0, ncat1-1);
	vPP2 = dvector(0, ncat2-1);
	vP121 = dvector(0, ncat1-1);
	vP212 = dvector(0, ncat2-1);
	M1 = dmatrix(0, ncatx-1, 0, ncat1-1);
	N1 = dmatrix(0, ncatv-1, 0, ncat1-1);
	M2 = dmatrix(0, ncaty-1, 0, ncat2-1);
	N2 = dmatrix(0, ncatv-1, 0, ncat2-1);
	r1 = dvector(0, ncatx-1);
	s1 = dvector(0, ncaty-1);
	r2 = dvector(0, ncatx-1);
	s2 = dvector(0, ncaty-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);
	v1 = dvector(0, ncatv-1);
	v2 = dvector(0, ncatv-1);
	Cr1 = dmatrix(0, ncat1-1, 0, npara1-1);
	Cr2 = dmatrix(0, ncat2-1, 0, npara2-1);
	Cr1Star = dmatrix(0, ncat1-1, 0, npara1-1);
	Cr2Star = dmatrix(0, ncat2-1, 0, npara2-1);
	B1 = dmatrix(0, npara1-1, 0, ncat1-1);
	B2 = dmatrix(0, npara2-1, 0, ncat2-1);
	UR = dmatrix(0, ncatx-1, 0, npara1-1);
	UP = dmatrix(0, ncatx-1, 0, npara1-1);
	UPStar = dmatrix(0, ncatx-1, 0, npara1-1);
	Utemp = dmatrix(0, ncatx-1, 0, npara1-1);
	VR = dmatrix(0, ncaty-1, 0, npara1-1);
	US = dmatrix(0, ncatx-1, 0, npara2-1);
	VS = dmatrix(0, ncaty-1, 0, npara2-1);
	UQ = dmatrix(0, ncaty-1, 0, npara2-1);
	UQStar = dmatrix(0, ncaty-1, 0, npara2-1);
	Vtemp = dmatrix(0, ncaty-1, 0, npara2-1);
	pl = dvector(0, ncatx-1);
	vPl = dvector(0, npara1-1);
	ql = dvector(0, ncaty-1);
	vQl = dvector(0, npara2-1);
	Fr = dvector(0, ncatx-1);
	Gs = dvector(0, ncaty-1);
	FrUR = dvector(0, npara1-1);
	GsVR = dvector(0, npara1-1);
	FrUS = dvector(0, npara2-1);
	GsVS = dvector(0, npara2-1);
    fitbdist1 = dmatrix(0, ncatx-1, 0, ncatv-1);
    fitbdist2 = dmatrix(0, ncaty-1, 0, ncatv-1);

    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. note that this B is the transpose of the design matrix in the loglinear 
	   model. */

	/*First, assign natrual design matrix first for rows corresponding to x */	
	for (i=0; i<cu1; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[i][k*(ncatx)+j] = pow((double)j,(double)(i+1));
			}
		}
	}

    /* then assign values to rows corresponding to v */
	currentrow1 = cu1 ;
	for (i=0; i<cv1; i++)	{
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[currentrow1+i][k*(ncatx)+j] = pow((double)k,(double)(i+1));
			}
		}
	}

    /* assign value to the last rows corresponding to the interaction terms */
	currentrow1 += cv1;
	for (i=0; i<cuv1; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				B1[currentrow1+i][k*(ncatx)+j] = pow((double)j, (double)interx[i]) * 
                                                 pow((double)k, (double)interv1[i]);
			}
		}
	}

    /* assign values to the rows for the 0 score on form X */
	currentrow1 += cuv1 ;
	for (j=0; j<ncatx; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (j==0) 
				B1[currentrow1][k*(ncatx)+j] = 1;
			else
				B1[currentrow1][k*(ncatx)+j] = 0;
		}
	}

    /* assign values to the rows for the 0 score on common set v */
	currentrow1++;
	for (j=0; j<ncatx; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (k==0) 
				B1[currentrow1][k*(ncatx)+j] = 1;
			else
				B1[currentrow1][k*(ncatx)+j] = 0;
		}
	}

	/*assign value to the rows for score gaps in x 	*/
	currentrow1++;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				if ((j%5)==0 && j!=0) 
					B1[currentrow1+i][k*(ncatx)+j] = pow((double)j,(double)i);
				else
					B1[currentrow1+i][k*(ncatx)+j] = 0;
			}
		}
	}

	/*assign value to the rows for score gaps in v */
	currentrow1 += 4;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncatx; j++)	{
			for (k=0; k<ncatv; k++)	{
				if (((k)%5)==0 && k!=0) 
					B1[currentrow1+i][k*(ncatx)+j] = pow((double)k,(double)i);
				else
					B1[currentrow1+i][k*(ncatx)+j] = 0;
			}
		}
	}




    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. note that this B2 is the transpose of the design matrix in the loglinear 
	   model. */

	/*First, assign natrual design matrix first for rows corresponding to y */	
	for (i=0; i<cu2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[i][k*(ncaty)+j] = pow((double)j,(double)(i+1));
			}
		}
	}

    /* then assign values to rows corresponding to v */
	currentrow2 = cu2 ;
	for (i=0; i<cv2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[currentrow2+i][k*(ncaty)+j] = pow((double)k,(double)(i+1));
			}
		}
	}

    /* assign value to the last rows corresponding to the interaction terms */
	currentrow2 += cv2;
	for (i=0; i<cuv2; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				B2[currentrow2+i][k*(ncaty)+j] = pow((double)j, (double)intery[i]) * 
                                                 pow((double)k, (double)interv2[i]);
			}
		}
	}


    /* assign values to the rows for the 0 score on form Y  */
	currentrow2 += cuv2 ;
	for (j=0; j<ncaty; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (j==0) 
				B2[currentrow2][k*(ncaty)+j] = 1;
			else
				B2[currentrow2][k*(ncaty)+j] = 0;
		}
	}

    /* assign values to the rows for the 0 score on common set V  */
	currentrow2++;
	for (j=0; j<ncaty; j++)	{
		for (k=0; k<ncatv; k++)	{
			if (k==0) 
				B2[currentrow2][k*(ncaty)+j] = 1;
			else
				B2[currentrow2][k*(ncaty)+j] = 0;
		}
	}

	/*assign value to the rows for score gaps in y 	*/
	currentrow2++;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				if ((j%5)==0 && j!=0) 
					B2[currentrow2+i][k*(ncaty)+j] = pow((double)j,(double)i);
				else
					B2[currentrow2+i][k*(ncaty)+j] = 0;
			}
		}
	}

	/*assign value to the rows for score gaps in v */
	currentrow2 += 4;
	for (i=0; i<=3; i++) {
		for (j=0; j<ncaty; j++)	{
			for (k=0; k<ncatv; k++)	{
				if (((k)%5)==0 && k!=0) 
					B2[currentrow2+i][k*(ncaty)+j] = pow((double)k,(double)i);
				else
					B2[currentrow2+i][k*(ncaty)+j] = 0;
			}
		}
	}

	sum1 = 0;
	sum2 = 0;
	for(j = 0; j <ncatv; j++) {
		for(i = 0; i <ncatx; i++) {
			fitbdist1[i][j] = bivar1->bfd[i][j] / np1;
			sum1 += fitbdist1[i][j];
		}
		for(i = 0; i <ncaty; i++) {
			fitbdist2[i][j] = bivar2->bfd[i][j] / np2;
			sum2 += fitbdist2[i][j];
		}
	}


	vPMN(ncatx, ncatv, fitbdist1, vP1, M1, N1);
	vPMN(ncaty, ncatv, fitbdist2, vP2, M2, N2);
	vPT(ncatx, ncatv, fitbdist1, vPP1);
	vPT(ncaty, ncatv, fitbdist2, vPP2);
	MatrixMultiVector(ncatx, ncat1, M1, vP1, r1);
	MatrixMultiVector(ncatv, ncat1, N1, vP1, v1);
	MatrixMultiVector(ncaty, ncat2, M2, vP2, s2);
	MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);

	ComputeCmatrixGen(ncat1, npara1, np1, B1, vP1, Cr1);
	ComputeCmatrixGen(ncat2, npara2, np2, B2, vP2, Cr2);


	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncatx; i++) {
			vP121[i+j*ncatx] = vP1[i+j*ncatx] * v2[j] / v1[j];
			for (k=0; k<npara1; k++) 
				Cr1Star[i+j*ncatx][k] = Cr1[i+j*ncatx][k] * v2[j] / v1[j];
		}
	}

	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncaty; i++) {
			vP212[i+j*ncaty] = vP2[i+j*ncaty] * v1[j] / v2[j];
			for (k=0; k<npara2; k++) 
				Cr2Star[i+j*ncatx][k] =  Cr2[i+j*ncatx][k] * v1[j] / v2[j];
		}
	}
	for (i=0; i<ncatx; i++) {
		for (k=0; k<npara1; k++) Utemp[i][k] = 0;
		for (k=0; k<npara2; k++) US[i][k] = 0;
	}

	for (i=0; i<ncaty; i++) {	
		for (k=0; k<npara2; k++) Vtemp[i][k] = 0;
		for (k=0; k<npara1; k++) VR[i][k] = 0;
	}

	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncatx; i++) pl[i] = fitbdist1[i][j];
		for (k=0; k<npara1; k++) {
			for (vPl[k]=0, i=0; i<ncatx; i++) 
				vPl[k] += Cr1[j*ncatx+i][k];
		}
		for (i=0; i<ncatx; i++) 
			for (k=0; k<npara1; k++) 
				Utemp[i][k] += v2[j] / v1[j] / v1[j] * pl[i] * vPl[k];

		for (i=0; i<ncaty; i++) ql[i] = fitbdist2[i][j];
		for (k=0; k<npara2; k++) {
			for (vQl[k]=0, i=0; i<ncaty; i++) 
				vQl[k] += Cr2[j*ncaty+i][k];
		}
		for (i=0; i<ncaty; i++) 
			for (k=0; k<npara2; k++) 
				Vtemp[i][k] += v1[j] / v2[j] /  v2[j]* ql[i] * vQl[k];

		for (i=0; i<ncaty; i++) 
			for (k=0; k<npara1; k++) VR[i][k] += wts / v2[j] * ql[i] * vPl[k];

		for (i=0; i<ncatx; i++) 
			for (k=0; k<npara2; k++) US[i][k] += (1 - wts) / v1[j] * pl[i] * vQl[k];
	}

	MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1, UP);
	MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2, UQ);
	MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1Star, UPStar);
	MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2Star, UQStar);
	MatrixMultiVector(ncatx, ncat1, M1, vP121, r2);
	MatrixMultiVector(ncatx, ncat2, M2, vP212, s1);


	for (i=0; i<ncatx; i++) 
		r[i] = wts * r1[i] + (1 - wts) * r2[i];
	for (i=0; i<ncaty; i++) 
		s[i] = wts * s1[i] + (1 - wts) * s2[i];

	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 
	hx = 1.9243;
	hy = 2.0056;
	
	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, Equatedx); 

	for (i=0; i<ncatx; i++) {
		for (k=0; k<npara1; k++) 
			UR[i][k] = wts * UP[i][k] + (1 - wts) * (UPStar[i][k] - Utemp[i][k]);
	}

	for (i=0; i<ncaty; i++) {
		for (k=0; k<npara2; k++) 
			VS[i][k] = (1 - wts) * UQ[i][k] + wts * (UQStar[i][k] - Vtemp[i][k]);
	}
  	
	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatx, scoresx, r, hx, Fr, scoresx[i]);
        PartialFPartialr(ncaty, scoresy, s, hy, Gs, Equatedx[i]);
        Gp = KernelContinuPdf(ncaty, scoresy, s, hy, Equatedx[i]);
		VectorMultiMatrix(ncatx, npara1, Fr, UR, FrUR);
		VectorMultiMatrix(ncaty, npara1, Gs, VR, GsVR);
		VectorMultiMatrix(ncatx, npara2, Fr, US, FrUS);
		VectorMultiMatrix(ncaty, npara2, Gs, VS, GsVS);
		for (j=0; j<npara1; j++) 
			FrUR[j] -= GsVR[j];
		for (j=0; j<npara2; j++) 
			FrUS[j] -= GsVS[j];
		SEE[i] = sqrt(VectorNormSq(npara1, FrUR) + VectorNormSq(npara2, FrUS)) / Gp;
	}
    	
	free_dvector(vP1, 0, ncat1-1);
	free_dvector(vP2, 0, ncat2-1);
	free_dvector(vP121, 0, ncat1-1);
	free_dvector(vP212, 0, ncat2-1);
	free_dmatrix(M1, 0, ncatx-1, 0, ncat1-1);
	free_dmatrix(N1, 0, ncatv-1, 0, ncat1-1);
	free_dmatrix(M2, 0, ncaty-1, 0, ncat2-1);
	free_dmatrix(N2, 0, ncatv-1, 0, ncat2-1);
	free_dvector(r1, 0, ncatx-1);
	free_dvector(s1, 0, ncaty-1);
	free_dvector(r2, 0, ncatx-1);
	free_dvector(s2, 0, ncaty-1);
	free_dvector(r, 0, ncatx-1);
	free_dvector(s, 0, ncaty-1);
	free_dvector(v1, 0, ncatv-1);
	free_dvector(v2, 0, ncatv-1);
	free_dmatrix(Cr1, 0, ncat1-1, 0, npara1-1);
	free_dmatrix(Cr2, 0, ncat2-1, 0, npara2-1);
	free_dmatrix(Cr1Star, 0, ncat1-1, 0, npara1-1);
	free_dmatrix(Cr2Star, 0, ncat2-1, 0, npara2-1);
	free_dmatrix(B1, 0, npara1-1, 0, ncat1-1);
	free_dmatrix(B2, 0, npara2-1, 0, ncat2-1);
	free_dmatrix(UR, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(UP, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(UPStar, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(Utemp, 0, ncatx-1, 0, npara1-1);
	free_dmatrix(VR, 0, ncaty-1, 0, npara1-1);
	free_dmatrix(US, 0, ncatx-1, 0, npara2-1);
	free_dmatrix(VS, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(UQ, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(UQStar, 0, ncaty-1, 0, npara2-1);
	free_dmatrix(Vtemp, 0, ncaty-1, 0, npara2-1);
	free_dvector(pl, 0, ncatx-1);
	free_dvector(vPl, 0, npara1-1);
	free_dvector(ql, 0, ncaty-1);
	free_dvector(vQl, 0, npara2-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);
	free_dvector(FrUR, 0, npara1-1);
	free_dvector(GsVR, 0, npara1-1);
	free_dvector(FrUS, 0, npara2-1);
	free_dvector(GsVS, 0, npara2-1);
	free_dmatrix(fitbdist1, 0, ncatx-1, 0, ncatv-1);
	free_dmatrix(fitbdist2, 0, ncaty-1, 0, ncatv-1);
}

/*--------------------------------------------------------------------------
	KernelEquateNEATPS
	
	functionality:

	Computes kernel equating function for the NEAT design 
	with post stratification method based on procedures in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 3/29/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         wts for group taking form X.
 
	output:
 		Equatedx    a vector containing the equated score 
--------------------------------------------------------------------------*/
void KernelEquateNEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx)
{
	int i, j, ncat1, ncat2, ncatx, ncaty, ncatv;
	long np1, np2;
	double *vP1, *vP2, **M1, **N1, **M2, **N2, *r1, *s1, *r2, *s2, *r, *s,  *v1, *v2;
	/* r1 is for X taken first(first group), r2 is for X taken second (second group),
	   s1 is for Y taken first(second group), s2 is for Y taken second (first group)
	   vP1 is for first group, vP2 is for second group */
	double *vP121, *vP212;
	double hx, hy;
	double *scoresx, *scoresy;
	double **fitbdist1, **fitbdist2;

	ncatx = bivar1->nsx;
	ncaty = bivar2->nsx;
	ncatv = bivar1->nsv;
	np1 = bivar1->num_persons;
	np2 = bivar2->num_persons;
	scoresx = dvector(0, ncatx);
	for (i=0;i<ncatx;i++) scoresx[i] = bivar1->minx + i * bivar1->incx;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar2->minx + i * bivar2->incx;
	ncat1 = ncatx * ncatv;
	ncat2 = ncaty * ncatv;

	vP1 = dvector(0, ncat1-1);
	vP2 = dvector(0, ncat2-1);
	vP121 = dvector(0, ncat1-1);
	vP212 = dvector(0, ncat2-1);
	M1 = dmatrix(0, ncatx-1, 0, ncat1-1);
	N1 = dmatrix(0, ncatv-1, 0, ncat1-1);
	M2 = dmatrix(0, ncaty-1, 0, ncat2-1);
	N2 = dmatrix(0, ncatv-1, 0, ncat2-1);
	r1 = dvector(0, ncatx-1);
	s1 = dvector(0, ncaty-1);
	r2 = dvector(0, ncatx-1);
	s2 = dvector(0, ncaty-1);
	r = dvector(0, ncatx-1);
	s = dvector(0, ncaty-1);
	v1 = dvector(0, ncatv-1);
	v2 = dvector(0, ncatv-1);
    fitbdist1 = dmatrix(0, ncatx-1, 0, ncatv-1);
    fitbdist2 = dmatrix(0, ncaty-1, 0, ncatv-1);


	for(j = 0; j <ncatv; j++) {
		for(i = 0; i <ncatx; i++) 
			fitbdist1[i][j] = bivar1->bfd[i][j] / np1;
		for(i = 0; i <ncaty; i++) 
			fitbdist2[i][j] = bivar2->bfd[i][j] / np2;
	}

	vPMN(ncatx, ncatv, fitbdist1, vP1, M1, N1);
	vPMN(ncaty, ncatv, fitbdist2, vP2, M2, N2);
	MatrixMultiVector(ncatx, ncat1, M1, vP1, r1);
	MatrixMultiVector(ncatv, ncat1, N1, vP1, v1);
	MatrixMultiVector(ncaty, ncat2, M2, vP2, s2);
	MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);


	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncatx; i++) {
			vP121[i+j*ncatx] = vP1[i+j*ncatx] * v2[j] / v1[j];
		/*	for (k=0; k<npara1; k++) 
				Cr1Star[i+j*ncatx][k] = Cr1[i+j*ncatx][k] * v2[j] / v1[j]; */
		}
	}

	for (j=0; j<ncatv; j++) {
		for (i=0; i<ncaty; i++) {
			vP212[i+j*ncaty] = vP2[i+j*ncaty] * v1[j] / v2[j];
		/*	for (k=0; k<npara2; k++) 
				Cr2Star[i+j*ncatx][k] =  Cr2[i+j*ncatx][k] * v1[j] / v2[j]; */
		}
	}
	MatrixMultiVector(ncatx, ncat1, M1, vP121, r2);
	MatrixMultiVector(ncatx, ncat2, M2, vP212, s1);


	for (i=0; i<ncatx; i++) 
		r[i] = wts * r1[i] + (1 - wts) * r2[i];
	for (i=0; i<ncaty; i++) 
		s[i] = wts * s1[i] + (1 - wts) * s2[i];

	hx = Optimalh(ncatx, scoresx, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 
	
	KernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, Equatedx); 

    	
	free_dvector(vP1, 0, ncat1-1);
	free_dvector(vP2, 0, ncat2-1);
	free_dvector(vP121, 0, ncat1-1);
	free_dvector(vP212, 0, ncat2-1);
	free_dmatrix(M1, 0, ncatx-1, 0, ncat1-1);
	free_dmatrix(N1, 0, ncatv-1, 0, ncat1-1);
	free_dmatrix(M2, 0, ncaty-1, 0, ncat2-1);
	free_dmatrix(N2, 0, ncatv-1, 0, ncat2-1);
	free_dvector(r1, 0, ncatx-1);
	free_dvector(s1, 0, ncaty-1);
	free_dvector(r2, 0, ncatx-1);
	free_dvector(s2, 0, ncaty-1);
	free_dvector(r, 0, ncatx-1);
	free_dvector(s, 0, ncaty-1);
	free_dvector(v1, 0, ncatv-1);
	free_dvector(v2, 0, ncatv-1);
	free_dmatrix(fitbdist1, 0, ncatx-1, 0, ncatv-1);
	free_dmatrix(fitbdist2, 0, ncaty-1, 0, ncatv-1);
}



/*--------------------------------------------------------------------------
	KernelEquateSEENEATChn
	
	functionality:

	Computes kernel equating function and SEE for the NEAT design 
	with chained equating method based on procedures in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 4/1/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         wts for group taking form X.
 
	output:
 		Equatedx    a vector containing the equated score 
		SEE         a vector containing standard error of equating 
--------------------------------------------------------------------------*/
void KernelEquateSEENEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *Equatedx, double *SEE)
{
	int i, j, k, ncat, ncatx, ncaty, ncatv, npara, cu, cv, cuv;
	long np;
	double *SEEA, *EquatedxA, *SEEY;
	double **fitbdist;
	double *vP, **M, **N, *r, *s, **U, **V;
	double hv, hy;
	double **Cr, **B;
	double *Fr, *Gs, GpQ, HpQ, *FrU, *GsV;
	double cdfv;
	double *scoresv, *scoresy;
	int *intery, *interv;
	int **cpm2;

	cuv = bivar2->cuv;
	cu = bivar2->cu;
	cv = bivar2->cv;
	cpm2 = imatrix(0, 1, 0, cuv);
    for(i=0;i<cuv;i++) for(j=0;j<2;j++) cpm2[j][i] = bivar2->cpm[i][j];

	ncatx = bivar1->nsx;
	SEEA = dvector(0, ncatx-1);
	EquatedxA = dvector(0, ncatx-1);
	SEEY = dvector(0, ncatx-1);

    KernelEquateSEESG(bivar1, EquatedxA, SEEA);	

	ncaty = bivar2->nsx;
	ncatv = bivar2->nsv;
	npara = bivar2->cu + bivar2->cv + bivar2->cuv;
	intery = cpm2[0];
	interv = cpm2[1];
	np = bivar2->num_persons;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar2->minx + i * bivar2->incx;
	scoresv = dvector(0, ncatv);
	for (i=0;i<ncatv;i++) scoresv[i] = bivar2->minv + i * bivar2->incv;
	ncat = ncatv * ncaty;

	Cr = dmatrix(0, ncat-1, 0, npara-1);
	B = dmatrix(0, npara-1, 0, ncat-1);
	U = dmatrix(0, ncatv-1, 0, npara-1);
	V = dmatrix(0, ncaty-1, 0, npara-1);
	Fr = dvector(0, ncatv-1);
	Gs = dvector(0, ncaty-1);
	FrU = dvector(0, npara-1);
	GsV = dvector(0, npara-1);
 	fitbdist = dmatrix(0, ncatv-1, 0, ncaty-1);

	for(i = 0; i <ncatv; i++) 
		for(j = 0; j <ncaty; j++) fitbdist[i][j] = bivar2->bfd[j][i] / np;


    /*The following code set a natrual design matrix corresponding to a natrual basis of 
       polynomials. This B matrix is the transpose of the design matrix in log-linear model
       First, assign natrual design matrix first for rows corresponding to x */
	
	for (i=0; i<cu; i++)
	{
		for (j=0; j<ncatv; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatv)+j] = pow((double)j,(double)(i+1));
			}
		}
	}

    /* then assign values to rows corresponding to y */
	for (i=cu; i<cu+cv; i++)
	{
		for (j=0; j<ncatv; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i][k*(ncatv)+j] = pow((double)k,(double)(i-cu+1));
			}
		}
	}

    /* assign value to the last columns corresponding to the interaction terms */
	for (i=0; i<cuv; i++)
	{
		for (j=0; j<ncatv; j++)
		{
			for (k=0; k<ncaty; k++)
			{
				B[i+cu+cv][k*(ncatv)+j] = pow((double)j, (double)interv[i]) * 
                                          pow((double)k, (double)intery[i]);
			}
		}
	}

	vP = dvector(0, ncat-1);
	M = dmatrix(0, ncatv-1, 0, ncat-1);
	N = dmatrix(0, ncaty-1, 0, ncat-1);
	r = dvector(0, ncatv-1);
	s = dvector(0, ncaty-1);

	vPMN(ncatv, ncaty, fitbdist, vP, M, N);
	MatrixMultiVector(ncatv, ncat, M, vP, r);
	MatrixMultiVector(ncaty, ncat, N, vP, s);

	ComputeCmatrixGen(ncat, npara, np, B, vP, Cr);
	MatrixMultiMatrix(ncatv, ncat, npara, M, Cr, U);
	MatrixMultiMatrix(ncaty, ncat, npara, N, Cr, V);

	hv = Optimalh(ncatv, scoresv, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 

	for (i=0; i<ncatx; i++) {
		cdfv = KernelContinuCdf(ncatv, scoresv, r, hv, EquatedxA[i]);
		Equatedx[i] = KernelInverseCdf(ncaty, scoresy, s, hy, cdfv);
	}

	for (i=0; i<ncatx; i++)	{
        PartialFPartialr(ncatv, scoresv, r, hv, Fr, EquatedxA[i]);
        PartialFPartialr(ncaty, scoresy, s, hy, Gs, Equatedx[i]);
        HpQ = KernelContinuPdf(ncatv, scoresv, r, hv, EquatedxA[i]);
        GpQ = KernelContinuPdf(ncaty, scoresy, s, hy, Equatedx[i]);
		VectorMultiMatrix(ncatv, npara, Fr, U, FrU);
		VectorMultiMatrix(ncaty, npara, Gs, V, GsV);
		for (j=0; j<npara; j++) FrU[j] -= GsV[j];
		SEEY[i] = sqrt(VectorNormSq(npara, FrU)) / GpQ;
		SEE[i] = sqrt(DSQR(HpQ / GpQ * SEEA[i]) + DSQR(SEEY[i]));
	}
 
	free_dmatrix(Cr, 0, ncat-1, 0, npara-1);
	free_dmatrix(B, 0, ncat-1, 0, ncat-1);
	free_dmatrix(U, 0, ncatx-1, 0, npara-1);
	free_dmatrix(V, 0, ncaty-1, 0, npara-1);
	free_dvector(Fr, 0, ncatx-1);
	free_dvector(Gs, 0, ncaty-1);
	free_dvector(FrU, 0, npara-1);
	free_dvector(GsV, 0, npara-1);
 	free_dmatrix(fitbdist, 0, ncatx-1, 0, ncaty-1);
	free_dvector(SEEA, 0, ncatx-1);
	free_dvector(EquatedxA, 0, ncatx-1);
	free_dvector(SEEY, 0, ncatx-1);
}


/*--------------------------------------------------------------------------
	KernelEquateNEATChn
	
	functionality:

	Computes kernel equating function and for the NEAT design 
	with chained equating method based on procedures in von Davier, 
	Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
	
	author: Tianyou Wang 4/1/2005.
	
	input:
		bivar1      information about bivariate distribution for form X
		bivar2      information about bivariate distribution for form Y
		wts         wts for group taking form X.
 
	output:
 		Equatedx    a vector containing the equated score 
--------------------------------------------------------------------------*/
void KernelEquateNEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *Equatedx)
{
	int i, j, ncat, ncatx, ncaty, ncatv;
	long np;
	double *EquatedxA;
	double **fitbdist;
	double *vP, **M, **N, *r, *s;
	double hv, hy;
	double cdfv;
	double *scoresv, *scoresy;
	
	ncatx = bivar1->nsx;
	EquatedxA = dvector(0, ncatx-1);

    KernelEquateSG(bivar1, EquatedxA);	

	ncaty = bivar2->nsx;
	ncatv = bivar2->nsv;
	np = bivar2->num_persons;
	scoresy = dvector(0, ncaty);
	for (i=0;i<ncaty;i++) scoresy[i] = bivar2->minx + i * bivar2->incx;
	scoresv = dvector(0, ncatv);
	for (i=0;i<ncatv;i++) scoresv[i] = bivar2->minv + i * bivar2->incv;
	ncat = ncatv * ncaty;

 	fitbdist = dmatrix(0, ncatv-1, 0, ncaty-1);

	for(i = 0; i <ncatv; i++) 
		for(j = 0; j <ncaty; j++) fitbdist[i][j] = bivar2->bfd[j][i] / np;

	vP = dvector(0, ncat-1);
	M = dmatrix(0, ncatv-1, 0, ncat-1);
	N = dmatrix(0, ncaty-1, 0, ncat-1);
	r = dvector(0, ncatv-1);
	s = dvector(0, ncaty-1);

	vPMN(ncatv, ncaty, fitbdist, vP, M, N);
	MatrixMultiVector(ncatv, ncat, M, vP, r);
	MatrixMultiVector(ncaty, ncat, N, vP, s);


	hv = Optimalh(ncatv, scoresv, r, 1); 
	hy = Optimalh(ncaty, scoresy, s, 1); 

	for (i=0; i<ncatx; i++) {
		cdfv = KernelContinuCdf(ncatv, scoresv, r, hv, EquatedxA[i]);
		Equatedx[i] = KernelInverseCdf(ncaty, scoresy, s, hy, cdfv);
	}

 	free_dmatrix(fitbdist, 0, ncatx-1, 0, ncaty-1);
	free_dvector(EquatedxA, 0, ncatx-1);
	free_dvector(vP, 0, ncat-1);
	free_dmatrix(M, 0, ncatv-1, 0, ncat-1);
	free_dmatrix(N, 0, ncaty-1, 0, ncat-1);
	free_dvector(r, 0, ncatv-1);
	free_dvector(s, 0, ncaty-1);

}

/************************************************************************/

void Wrapper_RK(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for doing kernel equating with RG design
    and log-linear smoothing smoothing
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'R'(random groups)
    method = 'E'(equipercentile equating)
    smoothing = 'K' (kernel 'smoothing')  
    *x = pointer to struct USTATS (new form)
    *y = pointer to struct USTATS (old form)
    *ullx = pointer to struct ULL_SMOOTH (new form)
    *ully = pointer to struct ULL_SMOOTH (old form)
    rep = replication number for bootstrap; should be set to 0
          for actual equating;  
    
  Output
    
    struct PDATA *inall:   populates selected values of inall 
    
    struct ERAW_RESULTS *r: populates            

      **eraw: equated raw scores;          
              method (rows) by raw score (columns) matrix
              of equated scores. Here there is only one method.
              So, memory allocated for eraw[][] is: 
              eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                             (loc(x->max,x->min,x>-inc)]
              because we are getting equated raw scores for x 
      **mts:  moments for equated raw scores           
      
  NOTE: If Wrapper_RK() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 

  Function calls other than C or NR utilities:                   
    KernelEquateSEERG()
    MomentsFromFD()  
                                                
  Tianyou Wang

  Date of last revision: 6/30/08       
*/
{ 
                          /* method name --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};         
  double maxx, maxy, *scoresx, *scoresy;
  int i;
  double SEE[200];

  scoresx = dvector(0, ullx->ns);
  scoresy = dvector(0, ully->ns);
  maxx = ullx->min + (ullx->ns - 1) * ullx->inc;
  maxy = ully->min + (ully->ns - 1) * ully->inc;
  for(i=0;i<ullx->ns;i++) scoresx[i] = i;
  for(i=0;i<ully->ns;i++) scoresy[i] = i;

  inall->rep = rep;               /* should be set to 0 for actual equating */
                    /* counting of replications done in Wrapper_Bootstrap() */ 
                    
  /* allocation and assignments for struct PDATA inall
     Note that for every assignment of the form inall->(var) = x->(var)
     or inall->(var) = y->(var), values vary depending on whether x or y 
	 is for actual equating or a bootstrap sample; all other values are 
	 the same for the actual equating and a bootstrap sample */
  
  if(inall->rep == 0){     /* no assignment or stor alloc for bootstrap reps */
    strcpy(inall->xfname,x->fname);
    strcpy(inall->yfname,y->fname);
    inall->x = x;
    inall->y = y;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    strcpy(inall->names[0],names[0]);
 
    inall->min = x->min;  
    inall->max = x->max;
    inall->inc = x->inc;
    inall->fdx = x->fd;
    inall->n = x->n;

    inall->ullx = ullx;
    inall->ully = ully;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,0,0,loc(inall->max,inall->min,inall->inc)); 
    r->mts = dmatrix(0,0,0,3);                    /* 0,3 is for the 4 moments */
  }
   
/* Compute equating results */
  KernelEquateSEERG(ullx->ns, ullx->c, ullx->num_persons, scoresx, ullx->density,  
	                ully->ns, ully->c, ully->num_persons, scoresy, ully->density,  
					r->eraw[0], SEE); 

  /* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

} 


/********************************************************************************/

void Print_RK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print RK results (RG; equipercentile; kernel with log-linear smoothing)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: None
                                                
  Tianyou Wang 

  Date of last revision: 6/30/08   
*/
  int i,j;
  
  fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
    
  fprintf(fp,"Kernel Equating with Random Groups Design\n"
	  	     "and Polynomial Log-Linear Smoothing");

  fprintf(fp,"\n\nLog-Linear Smoothing for new form %c: ",inall->x->id);
  fprintf(fp,"\n   number of degrees of smoothing = %2d",inall->ullx->c);
  fprintf(fp,"\nLog-Linear Smoothing for old form %c: ",inall->y->id);
  fprintf(fp,"\n   number of degrees of smoothing = %2d",inall->ully->c);
 
  fprintf(fp,"\n\nInput file for new form %c: %s\n",inall->x->id,inall->xfname);
  fprintf(fp,"Input file for old form %c: %s\n\n",inall->y->id,inall->yfname);

  for(i=1;i<=30;i++) fprintf(fp,"-");

 /* following code set up for any number of methods */

  fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) fprintf(fp," ");
  fprintf(fp,"Equated Raw Scores");
  fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"   Method %d:",j);
  fprintf(fp,"\nRaw Score (%c)  ",inall->x->id);
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"  %s",inall->names[j]);
  fprintf(fp,"\n\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][0]);
  fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][1]);
  fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][2]);
  fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][3]);
  
  fprintf(fp,"\n\n");
  for(j=1;j<=63;j++) fprintf(fp,"*");

}


/*******************************************************************************/

void Wrapper_SK(char design, char method, char smoothing,  struct BSTATS *xy, 
               struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall, 
			   struct ERAW_RESULTS *r)  
/*
  Wrapper for doing kernel equating with SG design
  and log-linear smoothing.  
  
  NOTE: This is for the SG design in which x and y do not share any items in 
  common, which means that functionally this is the external anchor case.  
  The bivariate log-linear smoothing procedure needs to know this. So, when
  Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
  anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
  (with anchor set to 0) and Wrapper_SL() can still be used, but convergence 
  of the smoothing algorithm may be compromised because of dependencies between
  x and y. 
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'S' (single group)
    method = 'E' (equipercentile equating)
    smoothing = 'K' (kernel 'smoothing')  
    xy = struct BSTATS 
	bllxy = struct BLL_SMOOTH
    rep = replication number for bootstrap; should be set to 0
          for actual equating;  
    
    NOTE: it is assumed that the first variable
          in xy is indeed x and the second variable is y.
		  For bllxy, data memebers with '_x' are for x,
		  and data members with '_v' are for y.  This somewhat
		  inconsistent notation arises because BLL_SMOOTH is
		  usually used with the CG design in which the variables
		  are more naturally designated x (or u) and v.

  Output
    
    struct PDATA inall   populates selected values of inall  

    struct ERAW_RESULTS *r: populates            

      **eraw: equated raw scores;          
              method (rows) by raw score (columns) matrix
              of equated scores. Here there is only one method.
              So, memory allocated for eraw[][] is: 
              eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                             (loc(x->max,x->min,x>-inc)]
              because we are getting equated raw scores for x 
      **mts:  moments for equated raw scores            
      
  NOTE: If Wrapper_SK() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 
                                            
  Function calls other than C or NR utilities:
    KernelEquateSG()
    MomentsFromFD()  
                                                
  Tianyou Wang

  Date of last revision: 6/30/08   
*/
{ 
                       /* method names --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};                    
  
  inall->rep = rep;                /* should be set to 0 for actual equating */
                   /* counting of replications done in Wrapper_Bootstrap() */ 
                    
 /* Allocation and assignments for struct PDATA inall>
    Note that for every assignment of the form inall->(var) = x->(var)
    or inall->(var) = y->(var), values vary depending on whether x or y 
	is for actual equating or a bootstrap sample; all other values are 
	the same for the actual equating and a bootstrap sample */
  
  if(inall->rep == 0){     /* no assignment or stor alloc for bootstrap reps */
    strcpy(inall->xyfname,xy->fname);
    inall->xy = xy;
	inall->bllxy = bllxy;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
	inall->anchor = 0;  /* implicitly, anchor is external for biv log-linear
						                        smoothing with the SG design */
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    strcpy(inall->names[0],names[0]);
 
    inall->min = xy->min1;  
    inall->max = xy->max1;
    inall->inc = xy->inc1;
    inall->fdx = xy->fd1;
    inall->n = xy->n;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,0,0,loc(inall->max,inall->min,inall->inc)); 
    r->mts = dmatrix(0,0,0,3);                    /* 0,3 is for the 4 moments */
  }
   
/* Compute equating results. Put x on scale of y.
   Note that in struct xy, '1' designates x and '2' designates y; 
   in struct bllxy, '_x' designates x and '_v' designates y. So: 
     xy->ns2 = number of score categories for y
	 xy->min2 = minimum score for y
	 xy->inc2 = increment for y
	 bllxy->crfd_v = log-linear smoothed cum rel fd for y
	 xy->ns1 = number of score categories for x
     bllxy->prd_x = log-linear smoothed PR dist for x
	 r->eraw[0] = y equivalents for x (output) */
  
  KernelEquateSG(bllxy, r->eraw[0]);
 /* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

  return;

}

/*******************************************************************************/

void Print_SK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print SK results (SG; equipercentile; kernel log-linear smoothing)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: None
                                                
  Tianyou Wang 

  Date of last revision: 6/30/08   
*/
{
  int i,j;
  
  fprintf(fp,"\n\n%s\n\n",tt);
  fprintf(fp,"Input filename:  %s\n\n",inall->xy->fname);
  
  if(inall->rep>0)  fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
    
  fprintf(fp,"kernel Equating with Single Group Design\n"
	  	     "and Polynomial Log-Linear Smoothing\n\n");

  fprintf(fp,"Bivariate Log-Linear Smoothing using the Multinomial Model\n\n");

  fprintf(fp,"Number of score categories for %c  = %d\n",
	         inall->xy->id1,inall->bllxy->nsu);
  fprintf(fp,"Number of score categories for %c  = %d\n",
	         inall->xy->id2,inall->bllxy->nsv);
  fprintf(fp," Total number of score categories = %d\n\n",inall->bllxy->ns);

  fprintf(fp,"Number of persons (i.e., total of frequencies) = %7d\n\n",
	         inall->bllxy->num_persons);

  fprintf(fp,"Polynomial degree for %c = %d\n",inall->xy->id1,inall->bllxy->cu);
  fprintf(fp,"Polynomial degree for %c = %d\n",inall->xy->id2,inall->bllxy->cv);
  fprintf(fp,"Number of cross-products = %d\n",inall->bllxy->cuv);
  if(inall->bllxy->cuv>0){
    fprintf(fp,"Cross-Product moments: ");
    for(i=0;i<inall->bllxy->cuv;i++) fprintf(fp,"  (%c%d,%c%d)%c",
                            inall->xy->id1,inall->bllxy->cpm[i][0],
							inall->xy->id2,inall->bllxy->cpm[i][1],
                            (inall->bllxy->cuv==i+1) ? ' ' : ',');
  }

  fprintf(fp,"\n\n");
  for(i=1;i<=30;i++) fprintf(fp,"-");

 /* following code set up for any number of methods */

  fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) fprintf(fp," ");
  fprintf(fp,"Equated Raw Scores");
  fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"   Method %d:",j);
  fprintf(fp,"\nRaw Score (%c)  ",inall->xy->id1);
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"  %s",inall->names[j]);
  fprintf(fp,"\n\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][0]);
  fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][1]);
  fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][2]);
  fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][3]);
  
  fprintf(fp,"\n\n");
  for(j=1;j<=63;j++) fprintf(fp,"*");

  return;

}



/*******************************************************************************/

void Wrapper_BK(char design, char method, char smoothing,  double wtsx, 
				double wtsy, struct BSTATS *xy1, struct BSTATS *xy2, 
				struct BLL_SMOOTH *bllxy1, struct BLL_SMOOTH *bllxy2, 
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r)  
/*
  Wrapper for doing kernel equating with SG with counter-balance design 
  and log-linear smoothing.  
  
  NOTE: This is for the SG with counter balance design in which x and y 
  in both group 1 and group 2 do not share any items in 
  common, which means that functionally this is the external anchor case.  
  The bivariate log-linear smoothing procedure needs to know this. So, when
  Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
  anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
  (with anchor set to 0) and Wrapper_SL() can still be used, but convergence 
  of the smoothing algorithm may be compromised because of dependencies between
  x and y. (To date my experience does not suggest this is a problem.)
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'B' (single group w/ counter balance)
    method = 'E' (equipercentile equating)
    smoothing = 'K' (kernel 'smoothing')  
	wtsx = weight of x for population 1
	wtsy = weight of y for population 2
    xy1 = struct BSTATS for group 1
    xy2 = struct BSTATS for group 2
	bllxy1 = struct BLL_SMOOTH for group 1
 	bllxy2 = struct BLL_SMOOTH for group 2
    rep = replication number for bootstrap; should be set to 0
          for actual equating;  
    
    NOTE: it is assumed that the first variable
          in xy is indeed x and the second variable is y.
		  For bllxy, data memebers with '_x' are for x,
		  and data members with '_v' are for y.  This somewhat
		  inconsistent notation arises because BLL_SMOOTH is
		  usually used with the CG design in which the variables
		  are more naturally designated x (or u) and v.

  Output
    
    struct PDATA inall   populates selected values of inall  

    struct ERAW_RESULTS *r: populates            

      **eraw: equated raw scores;          
              method (rows) by raw score (columns) matrix
              of equated scores. Here there is only one method.
              So, memory allocated for eraw[][] is: 
              eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                             (loc(x->max,x->min,x>-inc)]
              because we are getting equated raw scores for x 
      **mts:  moments for equated raw scores            
      
  NOTE: If Wrapper_SK() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 
                                            
  Function calls other than C or NR utilities:
    KernelEquateSEECB()
    MomentsFromFD()  
                                                
  Tianyou Wang

  Date of last revision: 6/30/08   
*/
{ 
                       /* method names --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};     
  double wts[2];
  double SEE[200];

  wts[0]=wtsx; wts[1]=wtsy;
  
  inall->rep = rep;                /* should be set to 0 for actual equating */
                   /* counting of replications done in Wrapper_Bootstrap() */ 
                    
 /* Allocation and assignments for struct PDATA inall>
    Note that for every assignment of the form inall->(var) = x->(var)
    or inall->(var) = y->(var), values vary depending on whether x or y 
	is for actual equating or a bootstrap sample; all other values are 
	the same for the actual equating and a bootstrap sample */
  
  if(inall->rep == 0){     /* no assignment or stor alloc for bootstrap reps */
/*    strcpy(inall->xyfname,xy->fname);
    inall->xy = xy;
	inall->bllxy = bllxy; */
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
	inall->anchor = 0;  /* implicitly, anchor is external for biv log-linear
						                        smoothing with the SG design */
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    strcpy(inall->names[0],names[0]);
 
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,0,0,loc(inall->max,inall->min,inall->inc)); 
    r->mts = dmatrix(0,0,0,3);                    /* 0,3 is for the 4 moments */
  }
   
/* Compute equating results. Put x on scale of y.
   Note that in struct xy, '1' designates x and '2' designates y; 
   in struct bllxy, '_x' designates x and '_v' designates y. So: 
     xy->ns2 = number of score categories for y
	 xy->min2 = minimum score for y
	 xy->inc2 = increment for y
	 bllxy->crfd_v = log-linear smoothed cum rel fd for y
	 xy->ns1 = number of score categories for x
     bllxy->prd_x = log-linear smoothed PR dist for x
	 r->eraw[0] = y equivalents for x (output) */
  
  KernelEquateSEECB(bllxy1, bllxy2, wts, r->eraw[0], SEE);
  
 /* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

  return;

}

/*******************************************************************************/


void Wrapper_CK(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for conitnuized log-linear equating for CG design with log-linear smoothing. 
  Equipercentile equating includes frequency estimation with 
  Braun-Holland (linear) results, modified frequency estimation with 
  Braun-Holland (linear) results, and chained equipercentile equating
    
  Assumes that in xv, score 1 is for x and score 2 is for v
  Assumes that in yv, score 1 is for y and score 2 is for v
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep
  
  Input:
  
    design = 'C' (CINEG)

    method:  'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
             'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
             'G' = FE + BH-FE + MFE + BH-MFE
             'C' = Chained
			 'H' = FE + BH-FE + Chained
             'A' = FE + BH-FE + MFE + BH-MFE + Chained
              
    smoothing = 'K' (kernel 'smoothing') 

    w1 = weight for pop. 1 (associated with xv)
         [0,1] except that for any number outside this 
         range, proportional weights are used -- i.e.,
         w1 = xv->n/(xv->n + yv->n)
    anchor = 0 --> external; 1 --> internal
    rv1 = reliability of common items for population 1 
          (set to 0 for all methods except 'F', 'G, and 'A')
    rv2 = reliability of common items for population 2
          (set to 0 for all methods except 'F', 'G, and 'A')
    xv = struct BSTATS
    yv = struct BSTATS 
    bllxv = struct BLL_SMOOTH; uses brfd, prd_x,  and crfd_v 
	bllyv = struct BLL_SMOOTH; uses brfd, crfd_x, and crfd_v
	        (Note that bllyv->crfd_x is really crfd for y in pop 2)
    rep = replication number for bootstrap; should be set to 0
          for actual equating; 
   
    NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be obtained 
    
  Output:
    
    struct PDATA inall:   populates selected values of inall 
    
    struct ERAW_RESULTS r          
 
      a[] = slopes for Braun-Holland
      b[] = intercepts for Braun-Holland
      eraw[][]:  equated raw scores
      mts[][]:  moments for equated raw scores   
          
      NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
            of equated scores; memory allocated here;
            eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
            because we are getting equated raw scores for x.
            eraw[][] stored in this "row"  manner so that 
            Equated_ss() can be used directly on 
            eraw[k] where k is the method number  
          
  NOTE: Whenever method differs, there must be different structures
        passed as struct PDATA and struct ERAW_RESULTS 
    
  NOTE: If Wrapper_CK() is called in a bootstrap loop, then in
        the calling function struct ERAW_RESULTS must be different
        from struct ERAW_RESULTS for the actual equating. 
                                            
  Function calls other than C or NR utilities:
    KernelEquateNEATPS()
    KernelEquateNEATChn()
    runerror()
                                               
  Tianyou Wang

  Date of last revision: 6/30/08   
*/
{ 
  int i;
  double *ptr;                                      /* pointer for eraw[] */
                       /* method names --- 10 characters; right justified */
  char *names[] ={"        FE", "       MFE", "  ChainedE"};                   
  
  inall->rep = rep;               /* should be set to 0 for actual equating. */
                    /* Counting of replications done in Wrapper_Bootstrap(), 
             which is why this statement cannot be in the if statement below */ 
                    
  /* allocation and assignments for inall
     Note that for every assignment of the form inall->(var) = xv->(var)
	 or inall->(var) = yv->(var) values vary depending on whether xv or yv
	 is for actual equating or a bootstrap sample; all other values are 
	 the same for the actual equating and a bootstrap sample */
  
  if(inall->rep == 0){   /* no assignment or stor alloc for bootstrap reps */
    strcpy(inall->xfname,xv->fname);
    strcpy(inall->yfname,yv->fname);
    inall->xv = xv;
    inall->yv = yv;
    inall->bllxv = bllxv;
	inall->bllyv = bllyv;
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
    inall->w1 = (w1<0 || w1>1) ? (double) (xv->n)/(xv->n + yv->n) : w1; 
                                   /* proportional wts if w1 outside [0,1] */          
    inall->anchor = anchor;
    inall->rv1 = rv1;
    inall->rv2 = rv2;

    if((method == 'F' || method == 'G' || method == 'A') && 
       (rv1 == 0 || rv2 == 0))
       runerror("\nMFE cannot be conducted since rv1 == 0 or rv2 == 0");
    
    inall->names = cmatrix(0,4,0,11);             /* maximum of five names */
 
    if(method == 'E'){                                    /* method == 'E' */
      inall->nm = 1;
      strcpy(inall->names[0],names[0]);
    }
    else if(method == 'F'){                               /* method == 'F' */
      inall->nm = 1;
      strcpy(inall->names[0],names[1]);
    }
    else if(method == 'G'){                               /* method == 'G' */
      inall->nm = 2;
      for(i=0;i<=1;i++) strcpy(inall->names[i],names[i]);
    }
    else if(method == 'C'){                               /* method == 'C' */
      inall->nm = 1;
      strcpy(inall->names[0],names[2]);
    }
	else if(method == 'H'){                               /* method == 'H' */
      inall->nm = 2;
	  strcpy(inall->names[0],names[0]);
      strcpy(inall->names[1],names[2]);
    }
    else{                                                 /* method == 'A' */
      inall->nm = 3;
      for(i=0;i<=2;i++) strcpy(inall->names[i],names[i]);
    }  

    inall->min = xv->min1;  
    inall->max = xv->max1;
    inall->inc = xv->inc1;
    inall->fdx = xv->fd1;
    inall->n = xv->n;
  }
                                                         
/* allocation and assignments for r */ 

  if(inall->rep <= 1){         /* no storage allocation for bootstrap reps >1 */
    r->eraw = dmatrix(0,inall->nm-1,0,loc(xv->max1,xv->min1,xv->inc1)); 
    r->mts = dmatrix(0,inall->nm-1,0,3);          /* 0,3 is for the 4 moments */
    r->fxs = dmatrix(0,1,0,loc(xv->max1,xv->min1,xv->inc1));
    r->gys = dmatrix(0,1,0,loc(yv->max1,yv->min1,yv->inc1));
  }
   
/* Equipercentile results, including Braun-Holland (BH) linear. 
   Note: For FE syn densities are in fxs[0] and gys[0]
         For MFE syn densities are in fxs[1] and gys[1] 
         For BH under FE, slope in a[0] and intercept in b[0]
         For BH under MFE, slope in a[1] and intercept in b[1] */

                                           /* FE + BH-FE in positions 0 and 1*/
  if(method == 'E' || method == 'G' || method == 'A' || method == 'H')     
	KernelEquateNEATPS(bllxv, bllyv, inall->w1, r->eraw[0]);

  if(method == 'C') ptr = r->eraw[0];
  else if(method == 'A') ptr = r->eraw[2];
  else if(method == 'H') ptr = r->eraw[1];
  else ptr = NULL;

  if(ptr)
    KernelEquateNEATChn(bllxv, bllyv, ptr);
                        
/* get moments */

  for(i=0;i<=inall->nm-1;i++) 
    MomentsFromFD(xv->min1,xv->max1,xv->inc1,r->eraw[i],xv->fd1,r->mts[i]);
  
  return;
} 

/****************************************************************************/

void Print_CK(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print CK results (CG; log-linear smoothed kernel equating)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: 
    score() 
                                            
  Tianyou Wang

  Date of last revision: 6/30/08     
*/
{
  int i,j,
      jlast = (inall->method=='A') ? 75 : 63;

  char x = inall->xv->id1,                            /* x = u[] + v for new form */
       v = inall->xv->id2,                           /* common items for new form */
       y = inall->yv->id1,                            /* y = w[] + z for old form */
       z = inall->yv->id2,                           /* common items for old form */
       u[3],                  /* string designating non-common items for new form */
       w[3];                  /* string designating non-common items for old form */

  if(inall->anchor){u[0] = x; u[1] = '\''; u[2] = '\0';}  /* internal anchor "x'" */
  else {u[0] = ' '; u[1] = x; u[2] = '\0';}               /* external anchor " x" */

  if(inall->anchor){w[0] = y; w[1] = '\''; w[2] = '\0';}  /* internal anchor "y'" */
  else {w[0] = ' '; w[1] = y; w[2] = '\0';}               /* external anchor " y" */
  
  fprintf(fp,"\n\n%s\n\n",tt);
  
  if(inall->rep>0)  fprintf(fp,"Bootstrap relication %d\n\n",inall->rep);
  
  fprintf(fp,"Kernel Equating with CINEG Design"
             " (%s Anchor; w1 = %7.5f)\n"
			 "and Polynomial Log-Linear Smoothing\n\n",
             (inall->anchor != 0) ? "Internal" : "External", inall->w1); 
  fprintf(fp,"Input file for %c%c and pop 1: %s\n",x,v,inall->xfname);
  fprintf(fp,"Input file for %c%c and pop 2: %s\n\n",y,z,inall->yfname);

  if(inall->anchor)
	  fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score; and\n"
				 "(2) variables for old form and pop 2 are identified as\n"
	             "    %c, %s, and %c, where %c = %s + %c; i.e., %c is the \n"
	             "    total score, %s is the non-common-item score, and\n"
			     "    %c is the common-item score.\n\n",
			     x,u,v,x,u,v,x,u,v,  y,w,z,y,w,z,y,w,z);
  else
	  fprintf(fp,"NOTE: In this output:\n"
	             "(1) variables for new form and pop 1 are identified as\n"
				 "    %c and %c, where %c does not include %c;\n"
				 "(2) variables for old form and pop 2 are identified as\n"
				 "    %c and %c, where %c does not include %c.\n\n",
			     x,v,x,v, y,z,y,z);

  /* information about x and v in pop 1 */

  fprintf(fp,"Number of persons for %c = %7d\n\n",x,inall->bllxv->num_persons);
  fprintf(fp,"Number of score categories for %c = %d\n",x,inall->bllxv->nsx);
  fprintf(fp,"Number of score categories for %c = %d\n",v,inall->bllxv->nsv);

  fprintf(fp,"Polynomial degree for %s = %d\n",u,inall->bllxv->cu);
  fprintf(fp,"Polynomial degree for %c = %d\n",v,inall->bllxv->cv);
  fprintf(fp,"Number of cross-products %s%c = %d\n",u,v,inall->bllxv->cuv);
  if(inall->bllxv->cuv>0) fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllxv->cuv;i++) fprintf(fp,"  (%s%d,%c%d)%c",
                               u,inall->bllxv->cpm[i][0],v,inall->bllxv->cpm[i][1],
                               (inall->bllxv->cuv==i+1) ? ' ' : ',');

 /* information about y and v in pop 2 */

  fprintf(fp,"\n\nNumber of persons for %c = %7d\n\n",y,inall->bllyv->num_persons);
  fprintf(fp,"Number of score categories for %c = %d\n",y,inall->bllyv->nsx);
  fprintf(fp,"Number of score categories for %c = %d\n",z,inall->bllyv->nsv);

  fprintf(fp,"Polynomial degree for %s = %d\n",w,inall->bllyv->cu);
  fprintf(fp,"Polynomial degree for %c = %d\n",z,inall->bllyv->cv);
  fprintf(fp,"Number of cross-products %s%c = %d\n",w,z,inall->bllyv->cuv);
  if(inall->bllyv->cuv>0) fprintf(fp,"  Cross-Product moments: ");
  for(i=0;i<inall->bllyv->cuv;i++) fprintf(fp,"  (%s%d,%c%d)%c",
                               w,inall->bllyv->cpm[i][0],z,inall->bllyv->cpm[i][1],
                               (inall->bllyv->cuv==i+1) ? ' ' : ',');
  
  /* print equated raw scores for methods */

  if(inall->method == 'E' || inall->method == 'F' || 
	 inall->method == 'G' || inall->method == 'A')
     for(i=1;i<=jlast;i++) fprintf(fp,"-");

  /* following code set up for any number of methods */

  fprintf(fp,"\n\n");
  for(j=1;j<=15+(inall->nm*12-18)/2;j++) fprintf(fp," ");
  fprintf(fp,"Equated Raw Scores");
  fprintf(fp,"\n\n               ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"   Method %d:",j);
  fprintf(fp,"\nRaw Score (%c)  ",x);
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"  %s",inall->names[j]);
  fprintf(fp,"\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    fprintf(fp,"\n %12.5f  ",score(i,inall->min,inall->inc));
    for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->eraw[j][i]);
  }
  
  fprintf(fp,"\n\n         Mean  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][0]);
  fprintf(fp,"\n         S.D.  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][1]);
  fprintf(fp,"\n         Skew  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][2]);
  fprintf(fp,"\n         Kurt  ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"%12.5f",r->mts[j][3]);
  
  fprintf(fp,"\n\n");
  for(j=1;j<=jlast;j++) fprintf(fp,"*");

}

/*****************************************************************************
  var

  Functionality:

  Compute the variance of an array of data.

  author: Tianyou Wang 1/5/2005.

  Input:
		data        pointer to a double array
		n           number of elements in the array

  Output:
        return the variance
*****************************************************************************/

double var(double *data, unsigned long n)
{
	unsigned long j;
	double vv, ave, s;

	for (ave=0.0,j=0;j<n;j++) ave += data[j];
	ave /= n;
	vv = 0.0;
	for (j=0;j<n;j++) {
		s = data[j] - ave;
		vv += s*s;
	}
	vv /= (n-1);
	return vv;
}

/*****************************************************************************
  mean

  Functionality:

  Compute the mean of an array of data.

  author: Tianyou Wang 1/5/2005.

  Input:
		data        pointer to a double array
		n           number of elements in the array

  Output:
        return the mean
*****************************************************************************/

double mean(double *data, unsigned long n)
{
	unsigned long j;
	double ave;

	for (ave=0.0,j=0;j<n;j++) ave += data[j];
	ave /= n;
	return ave;
}

