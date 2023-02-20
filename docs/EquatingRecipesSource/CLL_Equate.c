/*

CLL_Equate.c  File for continuized log-linear equating

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
#include <math.h>
#include <stdlib.h>
#include "ERutilities.h"
#include "NRutilities.h" 
#include "LogLinear.h"
#include "CG_EquiEquate.h"
#include "CLL_Equate.h"
#include "matrix.h"
/*--------------------------------------------------------------------------
	CLLEGPdf
	
	functionality

	Computes the continuous pdf based on the fitted loglinear model
	parameters for a discrete distribution
	
	Author: Tianyou Wang 10/29/2004.
	
	input
		min         lower limit of the distribution
		max         upper limit of the distribution
        npara       number of parameters
        para		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis.
		x           a particular score for which the pdf and cdf is 
		            generated   
		nc          normalizing constant 

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLEGPdf(double min, double max, int npara, double *para,  
					double x, double nc)
{
	double pdf;

	pdf = ExpPolynomial(npara, para, x);

	pdf /= nc;

	return pdf;

}


/*--------------------------------------------------------------------------
	CLLEGCdf
	
	functionality

	Computes the continuous pdf and cdf based on the fitted loglinear model
	parameters for a discrete distribution
	
	Author: Tianyou Wang 10/29/2004.
	
	input
		min         lower limit of the distribution
		max         upper limit of the distribution
        npara       number of parameters
        para		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis.
		x           a particular score for which the pdf and cdf is 
		            generated            
		nc          normalizing constant 

 
	output
 		the function returns the smoothed cdf
--------------------------------------------------------------------------*/

double CLLEGCdf(double min, double max, int npara, double *para,  
					double x, double nc)
{
	double cdf;

	cdf = GaussianQuadrature64(&ExpPolynomial, min, x, npara, para);

	cdf /= nc;

	return cdf;

}

/*--------------------------------------------------------------------------
	GaussianQuadrature16
	
	functionality

	Computes the numerical integration using Gaussian quadrature with 
	16 points
	
	Author: Tianyou Wang 10/29/2004.
	
	input
        (*func)()   pointer to a function which is the integrand
        a   		lower limit of integral
        b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature16(PTR_FTN2 func, double a, double b, int npara,
						  double *para)
{
	int j;
	double xr, xm, dx, s;
	static double x[]={0.09501250983764, 0.28160355077926, 0.45801677765723,
        0.61787624440264, 0.755404408355, 0.86563120238783, 
		0.94457502307323, 0.98940093499165};
	static double w[]={0.18945061045507, 0.18260341504492, 0.169156519395,
        0.14959598881658, 0.12462897125553, 0.09515851168249, 
		0.06225352393865, 0.02715245941175};

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j=0; j< 8; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

/*--------------------------------------------------------------------------
	GaussianQuadrature32
	
	functionality

	Computes the numerical integration using Gaussian quadrature with 
	32 points
	
	Author: Tianyou Wang 10/29/2004.
	
	input
        (*func)()   pointer to a function which is the integrand
        a   		lower limit of integral
        b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature32(PTR_FTN2 func, double a, double b, int npara,
						  double *para)
{
	int j;
	double xr, xm, dx, s;
	static double x[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	xm = 0.5 * (b + a);
	xr = 0.5 * (b - a);
	s = 0;
	for (j=0; j< 16; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

/*--------------------------------------------------------------------------
	GaussianQuadrature64
	
	functionality

	Computes the numerical integration using Gaussian quadrature with 
	64 points
	
	Author: Tianyou Wang 10/29/2004.
	
	input
        (*func)()   pointer to a function which is the integrand
        a   		lower limit of integral
        b   		upper limit of integral
		npara       number of parameters for the integrand
		para        vector of parameters for the integrand
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double GaussianQuadrature64(PTR_FTN2 func, double a, double b, int npara,
						  double *para)
{
	int j;
	double xr, xm, dx, s;
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
	s = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx)) + func(npara, para, (xm - dx)));
	}

	return s *= xr;
}

double GaussianQuadrature64i(PTR_FTN3 func, double a, double b, int npara,
						  double *para,  int i)
{
	int j;
	double xr, xm, dx, s;
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
	s = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s += w[j] * (func(npara, para, (xm + dx), i) + func(npara, para, (xm - dx), i));
	}

	return s *= xr;
}


/*--------------------------------------------------------------------------
	ExpPolynomial
	
	functionality

	Computes the exponential function of a polynomial (the fitted mean of 
	the loglinear model)
	
	Author: Tianyou Wang 10/29/2004.
	
	input
        npara       the number of polynomial coefficients
        *para  		the vector that contains the parameters
        x   		the plugged in value
 
	output
        The function returns exponential function of a polynomial.
	--------------------------------------------------------------------------*/

double ExpPolynomial(int npara, double *para, double x)
{
	int j;
	double s=0;

	for (j=0; j<npara; j++) 
		s += para[j] * pow(x, j);

	return exp(s);
}

/*--------------------------------------------------------------------------
	ExpPolynomialxi
	
	functionality

	Computes the exponential function of a polynomial (the fitted mean of 
	the loglinear model)
	
	Author: Tianyou Wang 10/29/2004.
	
	input
        npara       the number of polynomial coefficients
        *para  		the vector that contains the parameters
        x   		the plugged in value
 
	output
        The function returns exponential function of a polynomial.
	--------------------------------------------------------------------------*/

double ExpPolynomialxi(int npara, double *para, double x, int i)
{
	int j;
	double s=0;

	for (j=0; j<npara; j++) 
		s += para[j] * pow(x, j);

	s = exp(s);
	s *= pow(x, i);

	return s;
}


/*------------------------------------------------------------------------------
	CalcLLContinuMoments	
	
	functionality

	calculates mean, sd, skewness and kurtosis for a continuous distribution.

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

void CalcLLContinuMoments(PTR_FTN4 pdf, double a, double b, int npara,
						  double *para, double *moments)
{
	int j;
	double xr, xm, dx, s1, s2, s3, s4;
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
	s1 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s1 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * (xm + dx) + 
			pdf(a, b, npara, para, (xm - dx)) * (xm-dx));
	}
	s1 *= xr;

	s2 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s2 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * pow(xm + dx - s1, 2) + 
			pdf(a, b, npara, para, (xm - dx))* pow(xm - dx - s1, 2));
	}
	s2 *= xr;

	s3 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s3 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * pow(xm + dx - s1, 3) + 
			pdf(a, b, npara, para, (xm - dx))* pow(xm - dx - s1, 3));
	}
	s3 *= xr / pow(s2, 1.5);

	s4 = 0;
	for (j=0; j< 32; j++) {
		dx = xr * x[j];
		s4 += w[j] * (pdf(a, b, npara, para, (xm + dx)) * pow(xm + dx - s1, 4) + 
			pdf(a, b, npara, para, (xm - dx))* pow(xm - dx - s1, 4));
	}
	s4 *= xr / pow(s2, 2);

	moments[0] = s1;
	moments[1] = sqrt(s2);
	moments[2] = s3;
	moments[3] = s4;
}	


/*--------------------------------------------------------------------------
	CLLEquateEG
	
	functionality:

	Computes equating function based on continuized Log-linear cdf in Wang 	(2005). 
    	
	author: Tianyou Wang 1/5/2005.
	
	input:
		minx        lower limit of the distribution for the new form
		maxx        upper limit of the distribution for the new form
        nparax      number of parameters for the new form
        paraxx		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the new form.
		miny        lower limit of the distribution for the old form
		maxy        upper limit of the distribution for the old form
        nparay      number of parameters for the old form
        paraxy		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the old form.
		nDistCatx   Number of discrete score categories for the new form
		scoresx     vector containing the discrete scores for the new form

	output:
 		Equatedx   a vector containing the equated score 
--------------------------------------------------------------------------*/
void CLLEquateEG(double minx, double maxx, int nparax, double *parax,  
			   double miny, double maxy, int nparay, double *paray,	
			   int nDistCatx, double *scoresx, double *Equatedx)
{
	int i;
	double cdfx, ncx, ncy;

	ncx = GaussianQuadrature64(&ExpPolynomial, minx, maxx, nparax, parax);
	ncy = GaussianQuadrature64(&ExpPolynomial, miny, maxy, nparay, paray);

	for (i=0; i<nDistCatx; i++) {
		cdfx = CLLEGCdf(minx, maxx, nparax, parax, scoresx[i], ncx);
		Equatedx[i] = CLLInverseCdf(miny, maxy, nparay, paray, cdfx, ncy);
	}

}

/*--------------------------------------------------------------------------
	CLLInverseCdf
	
	functionality:

	Computes the inverse of the cdf for the continuized log-linear cdf in Wang 
	(2005). 
	
	author: Tianyou Wang 1/5/2005.
	
	input:
		min         lower limit of the distribution
		max         upper limit of the distribution
        npara       number of parameters
        para		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural basis.
		cdf         a particular cdf for which the score is found
 
	output:
 		The function returns the inverse of cdf
--------------------------------------------------------------------------*/
double CLLInverseCdf(double min, double max, int npara, double *para, double cdf, double nc)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .000001;

	lb = min;
	ub = max;
	cdfl = CLLEGCdf(min, max, npara, para, lb, nc);
	cdfu = CLLEGCdf(min, max, npara, para, ub, nc);
	if (cdf<cdfl) 
		return min;
	else if (cdf>cdfu)
		return max;
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = CLLEGCdf(min, max, npara, para, half, nc);  
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
	CLLBivPdf
	
	functionality

	Computes the continuous bivariate pdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivPdf(struct BLL_SMOOTH *bivar, double x, double y)
{
	static double integ;
	double minx, maxx, miny, maxy, pdf;

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;

	integ = BivGaussianQuadrature64(&BivExpPolynomial, bivar, minx, maxx, miny, maxy);
	pdf = BivExpPolynomial(bivar, x, y)/ integ;


	return pdf;

}

/*--------------------------------------------------------------------------
	CLLBivCdf
	
	functionality

	Computes the continuous bivariate cdf based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLBivCdf(struct BLL_SMOOTH *bivar, double x, double y, double nc)
{
	double minx, maxx, miny, maxy, cdf;

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;

	cdf = BivGaussianQuadrature32(&BivExpPolynomial, bivar, minx, x, miny, y);

	cdf /= nc;  /* bivar->num_persons; */

	return cdf;

}

/*--------------------------------------------------------------------------
	CLLMargYPdf
	
	functionality

	Computes the continuous marginal pdf for Y based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		y           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargYPdf(struct BLL_SMOOTH *bivar, double y, double nc)
{
	int j;
	double minx, maxx, miny, maxy;
	double xr, xm, dx, s;
	static double x[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;
	xm = 0.5 * (maxx + minx);
	xr = 0.5 * (maxx - minx);
	s = 0;
	for (j=0; j< 16; j++) {
		dx = xr * x[j];
		s += w[j] * (BivExpPolynomial(bivar, xm + dx, y) + BivExpPolynomial(bivar, xm - dx, y));
	}

	s *= xr;



	s /= nc; 


	return s;

}

/*--------------------------------------------------------------------------
	CLLMargXPdf
	
	functionality

	Computes the continuous marginal pdf for X based on the fitted bivariate 
	loglinear model	parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLMargXPdf(struct BLL_SMOOTH *bivar, double x, double nc)
{
	int j;
	double minx, maxx, miny, maxy;
	double yr, ym, dy, s;
	static double y[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minx = bivar->minx - .5;
	maxx = bivar->minx + (bivar->nsx - 1) * bivar->incx + .5;
	miny = bivar->minv - .5;
	maxy = bivar->minv + (bivar->nsv - 1) * bivar->incv + .5;
	ym = 0.5 * (maxy + miny);
	yr = 0.5 * (maxy - miny);
	s = 0;
	for (j=0; j< 16; j++) {
		dy = yr * y[j];
		s += w[j] * (BivExpPolynomial(bivar, x, ym + dy) + BivExpPolynomial(bivar, x, ym - dy));
	}

	s *= yr;

	s /= nc;

	return s;
}

/*--------------------------------------------------------------------------
	BivExpPolynomial
	
	functionality

	Computes the bivariate exponential function of a polynomial (the fitted mean of 
	the loglinear model)
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        x   		the plugged in value for X
        y   		the plugged in value for Y
 
	output
        The function returns exponential function of a polynomial.
	--------------------------------------------------------------------------*/

double BivExpPolynomial(struct BLL_SMOOTH *bivar, double x, double y)
{
	int j, nx, ny, nxbyy;
	double s=0;

	/* when internal anchor, f(x,v) = f(u,v) for u = x - v and u is within valid range */
	if (bivar->anchor==1) {
		x -= y; 
		if (x < bivar->minu) 
			return 0;
		if (x > (bivar->minu + (bivar->nsu - 1) * bivar->incu)) 
			return 0;
	}
	nx = bivar->cu;
	ny = bivar->cv;
	nxbyy = bivar->cuv;
	s = bivar->ap;
	for (j=1; j<=nx; j++) 
		s += bivar->Beta[j-1] * pow(x, j);
	for (j=1; j<=ny; j++) 
		s += bivar->Beta[nx+j-1] * pow(y, j);
	for (j=1; j<=nxbyy; j++) 
		s += bivar->Beta[nx+ny+j-1] * pow(x, bivar->cpm[j-1][0]) * pow(y, bivar->cpm[j-1][1]);

	return exp(s);
}

/*--------------------------------------------------------------------------
	BivGaussianQuadrature64
	
	functionality

	Computes the tow-dimensional numerical integration using Gaussian quadrature with 
	64 points
	
	Author: Tianyou Wang 4/29/2007.
	
	input
        (*func)()   pointer to a function which is the integrand
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        ax   		lower limit of x variable
        bx   		upper limit of x variable
        ay   		lower limit of y variable
        by   		upper limit of y variable
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double BivGaussianQuadrature64(PTR_FTN5 func, struct BLL_SMOOTH *bivar, double ax, 
							   double bx, double ay, double by)
{
	int i, j;
	double xr, xm, dx, yr, ym, dy, si, sj1, sj2;
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

	xm = 0.5 * (bx + ax);
	xr = 0.5 * (bx - ax);
	ym = 0.5 * (by + ay);
	yr = 0.5 * (by - ay);
	si = 0;
	for (i=0; i< 32; i++) {
		dx = xr * x[i];
		sj1 = 0;
		sj2 = 0;
		for (j=0; j< 32; j++) {
			dy = yr * x[j];
			sj1 += yr * w[j] * (func(bivar, (xm + dx), (ym + dy)) + func(bivar, (xm + dx), (ym - dy)));
		}
		for (j=0; j< 32; j++) {
			dy = yr * x[j];
			sj2 += yr * w[j] * (func(bivar, (xm - dx), (ym + dy)) + func(bivar, (xm - dx), (ym - dy)));
		}
		si += w[i] * (sj1  + sj2);
	}

	si *=  xr ;

	return si;
}

/*--------------------------------------------------------------------------
	BivGaussianQuadrature32
	
	functionality

	Computes the tow-dimensional numerical integration using Gaussian quadrature with 
	64 points
	
	Author: Tianyou Wang 4/29/2007.
	
	input
        (*func)()   pointer to a function which is the integrand
		bivar       the structure that contain bivariate distribution parameters
		            defined in "MLBivLogLin.h"
        ax   		lower limit of x variable
        bx   		upper limit of x variable
        ay   		lower limit of y variable
        by   		upper limit of y variable
 
	output
        The function returns the integrated value.
	--------------------------------------------------------------------------*/

double BivGaussianQuadrature32(PTR_FTN5 func, struct BLL_SMOOTH *bivar, double ax, double bx,
							   double ay, double by)
{
	int i, j;
	double xr, xm, dx, yr, ym, dy, si, sj1, sj2;
	static double x[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	xm = 0.5 * (bx + ax);
	xr = 0.5 * (bx - ax);
	ym = 0.5 * (by + ay);
	yr = 0.5 * (by - ay);
	si = 0;
	for (i=0; i< 16; i++) {
		dx = xr * x[i];
		sj1 = 0;
		sj2 = 0;
		for (j=0; j< 16; j++) {
			dy = yr * x[j];
			sj1 += yr * w[j] * (func(bivar, (xm + dx), (ym + dy)) + func(bivar, (xm + dx), (ym - dy)));
		}
		for (j=0; j< 16; j++) {
			dy = yr * x[j];
			sj2 += yr * w[j] * (func(bivar, (xm - dx), (ym + dy)) + func(bivar, (xm - dx), (ym - dy)));
		}
		si += w[i] * (sj1  + sj2);
	}

	si *=  xr ;

	return si;
}

/*--------------------------------------------------------------------------
	CLLNEATPSMargPdf
	
	functionality

	Computes the continuous marginal pdf for X for the synthetic population
	using the basic assumption of the frequency estimation based on the fitted 
	bivariate loglinear model parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wts         weights for population 1 in the synthetic population
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLNEATPSMargPdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *fv, 
						double wts, double x)
{
	int j;
	static double nc;
	double minv, maxv, minx, maxx;
	double vr, vm, dv, s, fxv, fx1, fx2;
	static double v[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minv = bivar1->minv - .5;
	maxv = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	minx = bivar1->minx - .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	vm = 0.5 * (maxv + minv);
	vr = 0.5 * (maxv - minv);
	nc = BivGaussianQuadrature64(&BivExpPolynomial, bivar1, minx, maxx, minv, maxv);
	s = 0;
	for (j=0; j< 16; j++) {
		dv = vr * v[j];

		/* nc was the last parameter in CLLBivPdf(), but not in function;
		   nc eliminated by rlb on 8-18-09=8 */

		fxv = CLLBivPdf(bivar1, x,  vm + dv) / fv[j] * fv[j+32]       
		    + CLLBivPdf(bivar1, x,  vm - dv) / fv[j+16] * fv[j+48];
		s += w[j] * fxv;
	}

	fx2 = s * vr;

	fx1 = CLLMargXPdf(bivar1, x, nc);

	//integ1 = BivGaussianQuadrature64(&BivExpPolynomial, bivar1, minx1, maxx1, minv, maxv);
	//integ2 = BivGaussianQuadrature64(&BivExpPolynomial, bivar2, minx2, maxx2, minv, maxv);


	return fx1 * wts  + fx2 * (1-wts) ;
}

/*--------------------------------------------------------------------------
	CLLNEATPSMargCdf
	
	functionality

	Computes the continuous marginal cdf for X for the synthetic population
	using the basic assumption of the frequency estimation based on the fitted 
	bivariate loglinear model parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wts         weights for population 1 in the synthetic population
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLNEATPSMargCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
								 double x)
{
	int j;
	static double nc1, nc2;
	double minv, maxv, minx, maxx, miny, maxy;
	double vr, vm, dv, xr, xm, dx, fx, s;
	double fv[4*16];  /* an array containing 4 arrays of V distribution at quadrature points,
					     1: f1v plus, 2: f1v minxs, 3: f2v plus, 4: f2v minxs, */
//	double sum1=0, sum2=0;

	static double xx[]={0.04830766568774, 0.14447196158280, 0.23928736225214, 
        0.33186860228213, 0.42135127613064, 0.50689990893223, 
        0.58771575724076, 0.66304426693022, 0.73218211874029, 
        0.79448379596794, 0.84936761373257, 0.89632115576605, 
        0.93490607593774, 0.96476225558751, 0.98561151154527, 
        0.99726386184948};
	static double w[]={0.09654008851473, 0.09563872007927, 0.09384439908080, 
        0.09117387869576, 0.08765209300440, 0.08331192422695, 
        0.07819389578707, 0.07234579410885, 0.06582222277636, 
        0.05868409347854, 0.05099805926238, 0.04283589802223,
        0.03427386291302, 0.02539206530926, 0.01627439473091, 
        0.00701861000947};

	minx = bivar1->minx - .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	miny = bivar2->minx - .5;
	maxy = bivar2->minx + bivar2->incx * (bivar2->nsx - 1) + .5;
	xm = 0.5 * (x + minx);
	xr = 0.5 * (x - minx);
	minv = bivar1->minv - .5;
	maxv = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	vm = 0.5 * (maxv + minv);
	vr = 0.5 * (maxv - minv);
	s = 0;
	nc1 = BivGaussianQuadrature64(&BivExpPolynomial, bivar1, minx, maxx, minv, maxv);
	nc2 = BivGaussianQuadrature64(&BivExpPolynomial, bivar2, miny, maxy, minv, maxv);
	for (j=0; j< 16; j++) {
		dv = vr * xx[j];
		fv[j] = CLLMargYPdf(bivar1, vm + dv, nc1);
		fv[j+16] = CLLMargYPdf(bivar1, vm - dv, nc1);
		fv[j+32] = CLLMargYPdf(bivar2, vm + dv, nc2);
		fv[j+48] = CLLMargYPdf(bivar2, vm - dv, nc2);
	}
	for (j=0; j< 16; j++) {
		dx = xr * xx[j];
		fx = CLLNEATPSMargPdf(bivar1, bivar2, fv, wts, xm + dx) + 
			 CLLNEATPSMargPdf(bivar1, bivar2, fv, wts, xm - dx);
		s += w[j] * fx;
	}

	return s *= xr;
}

/*--------------------------------------------------------------------------
	CLLNEATPSInverseCdf
	
	functionality

	Computes the continuous marginal cdf for X for the synthetic population
	using the basic assumption of the frequency estimation based on the fitted 
	bivariate loglinear model parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wts         weights for population 1 in the synthetic population
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

double CLLNEATPSInverseCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, 
										double wts, double cdf)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .00001;

	lb = bivar1->minx-0.5;
	ub = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	cdfl = CLLNEATPSMargCdf(bivar1, bivar2, wts, lb);
	cdfu = CLLNEATPSMargCdf(bivar1, bivar2, wts, ub);
	if (cdf<cdfl) 
		return lb;
	else if (cdf>cdfu)
		return ub;
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = CLLNEATPSMargCdf(bivar1, bivar2, wts, half);  
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
	CLLEquateNEATPS
	
	functionality

	Computes the equating fucntion for the NEAT design using the poststratification
	(frequency estimation) + CLL method
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wts         weight for the X group in the synthetic population            

 
	output
 		Equatedx    The equating function 
--------------------------------------------------------------------------*/

int CLLEquateNEATPS(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wts, 
						   double *Equatedx)
{
	int i, nDistCatx;
	static double integ;

	double cdfx;

	nDistCatx = bivar1->nsx;
	for (i=0; i<nDistCatx; i++) {
		cdfx = CLLNEATPSMargCdf(bivar1, bivar2, wts, bivar1->minx + bivar1->incx * i);
		Equatedx[i] = CLLNEATPSInverseCdf(bivar2, bivar1, 1-wts, cdfx);
	}

	return i;

}

/*--------------------------------------------------------------------------
	CLLMargInverseYCdf
	
	functionality

	Computes the inverse of continuous marginal cdf for Y for a bivariate distribution.
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		ycdf        the cdf for y that is used to compute the inverse            

 
	output
 		the function returns y score correspdonding the the ycdf
--------------------------------------------------------------------------*/

double CLLMargInverseYCdf(struct BLL_SMOOTH *bivar1, double ycdf, double nc)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .00001;
	double maxx;

	lb = bivar1->minv - .5;
	ub = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	cdfl = CLLBivCdf(bivar1, maxx, lb, nc);
	cdfu = CLLBivCdf(bivar1, maxx, ub, nc);
	if (ycdf<cdfl) 
		return lb;
	else if (ycdf>cdfu)
		return ub;
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = CLLBivCdf(bivar1, maxx, half, nc);  
            absdif = fabs(ycdf - cdfhalf);
            if (absdif < eps)
                return half;
			else if (cdfhalf<ycdf) 
				lb = half;
			else
				ub = half;
        } while (niter <=200);
	}

	return half;
}

/*--------------------------------------------------------------------------
	CLLMargInverseXCdf
	
	functionality

	Computes the inverse of continuous marginal cdf for X for a bivariate distribution.
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		xcdf        the cdf for x that is used to compute the inverse            

 
	output
 		the function returns y score correspdonding the the ycdf
--------------------------------------------------------------------------*/

double CLLMargInverseXCdf(struct BLL_SMOOTH *bivar1, double xcdf, double nc)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .00001;
	double maxy;

	lb = bivar1->minx - .5;
	ub = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	maxy = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	cdfl = CLLBivCdf(bivar1, lb, maxy, nc);
	cdfu = CLLBivCdf(bivar1, ub, maxy, nc);
	if (xcdf<cdfl) 
		return lb;
	else if (xcdf>cdfu)
		return ub;
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = CLLBivCdf(bivar1, half, maxy, nc);  
            absdif = fabs(xcdf - cdfhalf);
            if (absdif < eps)
                return half;
			else if (cdfhalf<xcdf) 
				lb = half;
			else
				ub = half;
        } while (niter <=200);
	}

	return half;
}

/*--------------------------------------------------------------------------
	CLLMargInverseCBYCdf
	
	functionality

	Computes the inverse of continuous marginal cdf for Y for the weighted
	bivariate distribution of the CB design.
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wtsy        the weight for population 1 and Y
		ycdf        the cdf for y that is used to compute the inverse            

 
	output
 		the function returns y score correspdonding the the ycdf
--------------------------------------------------------------------------*/

double CLLMargInverseCBYCdf(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wtsy, double ycdf, double nc1, double nc2)
{
	int niter=0;
	double ub, lb, half, cdfu, cdfl, cdfhalf;
	double absdif, eps = .00001;
	double maxx;

	lb = bivar1->minv - .5;
	ub = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	cdfl = wtsy * CLLBivCdf(bivar1, maxx, lb, nc1) + (1 - wtsy) * CLLBivCdf(bivar2, maxx, lb, nc2);
	cdfu = wtsy * CLLBivCdf(bivar1, maxx, ub, nc1) + (1 - wtsy) * CLLBivCdf(bivar2, maxx, ub, nc2);
	if (ycdf<cdfl) 
		return lb;
	else if (ycdf>cdfu)
		return ub;
    else
    {
        do
        {
            niter++;
            half = .5 * (lb + ub);
            cdfhalf = wtsy * CLLBivCdf(bivar1, maxx, half, nc1) + (1 - wtsy) * CLLBivCdf(bivar2, maxx, half, nc2);  
            absdif = fabs(ycdf - cdfhalf);
            if (absdif < eps)
                return half;
			else if (cdfhalf<ycdf) 
				lb = half;
			else
				ub = half;
        } while (niter <=200);
	}

	return half;
}


/*--------------------------------------------------------------------------
	CLLEquateNEATChn
	
	functionality

	Computes the equating fucntion for the NEAT design using the chained 
	equipercentile + CLL method
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
 
	output
 		Equatedx    The equating function 
--------------------------------------------------------------------------*/

int CLLEquateNEATChn(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double *Equatedx)
{
	int i, nDistCatx;
	static double nc1, nc2;
	double minv, maxv;
	double minx, maxx, miny, maxy;
	double cdfx, cdfv;
	double eqtempv;


	minx = bivar1->minx - .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	miny = bivar2->minx - .5;
	maxy = bivar2->minx + bivar2->incx * (bivar2->nsx - 1) + .5;
	minv = bivar1->minv - .5;
	maxv = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	nc1 = BivGaussianQuadrature64(&BivExpPolynomial, bivar1, minx, maxx, minv, maxv);
	nc2 = BivGaussianQuadrature64(&BivExpPolynomial, bivar2, miny, maxy, minv, maxv);
	maxv = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	maxy = bivar2->minx + bivar2->incx * (bivar2->nsx - 1) + .5;
	nDistCatx = bivar1->nsx;
	for (i=0; i<nDistCatx; i++) {
		cdfx = CLLBivCdf(bivar1, bivar1->minx + bivar1->incx * i, maxv, nc1);
		eqtempv = CLLMargInverseYCdf(bivar1, cdfx, nc1);
		cdfv = CLLBivCdf(bivar2, maxx, eqtempv, nc2);
		Equatedx[i] = CLLMargInverseXCdf(bivar2, cdfv, nc2);
	}

	return i;

}

/*--------------------------------------------------------------------------
	CLLEquateSG
	
	functionality

	Computes the equating fucntion for the single group design using the  
	equipercentile + CLL method
	
	Author: Tianyou Wang 7/10/2007.
	
	input
		bivar      the structure that contain bivariate distribution of X and 
		           Y (defined in "MLBivLogLin.h") 
 
	output
 		Equatedx    The equating function 
--------------------------------------------------------------------------*/

int CLLEquateSG(struct BLL_SMOOTH *bivar, double *Equatedx)
{
	int i, nDistCatx;
	static double nc;
	double minx, miny, maxx, maxy;
	double cdfx;

	minx = bivar->minx;
	miny = bivar->minv;
	maxx = minx + (bivar->nsx - 1) * bivar->incx + .5;
	maxy = miny + (bivar->nsv - 1) * bivar->incv + .5;
	nc = BivGaussianQuadrature64(&BivExpPolynomial, bivar, minx, maxx, miny, maxy);
	nDistCatx = bivar->nsx;
	for (i=0; i<nDistCatx; i++) {
		cdfx = CLLBivCdf(bivar, bivar->minx + bivar->incx * i, maxy, nc);
		Equatedx[i] = CLLMargInverseYCdf(bivar, cdfx, nc);

	}

	return i;

}

/*--------------------------------------------------------------------------
	CLLEquateCB
	
	functionality

	Computes the equating fucntion for the counter balance design using the  
	equipercentile + CLL method
	
	Author: Tianyou Wang 7/10/2007.
	
	input
		bivar1      the structure that contain bivariate distribution (X and Y) 
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution (X and Y)
		            (defined in "MLBivLogLin.h") for population 2
	    wtsx        weight for population 1 and X
		wtsy        weight for population 1 and Y
 
	output
 		Equatedx    The equating function 
--------------------------------------------------------------------------*/

int CLLEquateCB(struct BLL_SMOOTH *bivar1, struct BLL_SMOOTH *bivar2, double wtsx, double wtsy,
				double *Equatedx)
{
	int i, nDistCatx;
	double maxx, maxy;
	double cdfx, nc1, nc2;


	maxx = bivar1->minx + bivar1->incx * (bivar1->nsx - 1) + .5;
	maxy = bivar1->minv + bivar1->incv * (bivar1->nsv - 1) + .5;
	nDistCatx = bivar1->nsx;
	nc1 = BivGaussianQuadrature32(&BivExpPolynomial, bivar1, bivar1->minx, maxx, bivar1->minv, maxy);
	nc2 = BivGaussianQuadrature32(&BivExpPolynomial, bivar2, bivar1->minx, maxx, bivar1->minv, maxy);
	for (i=0; i<nDistCatx; i++) {
		cdfx = wtsx * CLLBivCdf(bivar1, bivar1->minx + bivar1->incx * i, maxy, nc1)
			+ (1 - wtsx) * CLLBivCdf(bivar2, bivar1->minx + bivar1->incx * i, maxy, nc2);
		Equatedx[i] = CLLMargInverseCBYCdf(bivar1, bivar2, wtsy, cdfx, nc1, nc2);
	}

	return i;

}

/*--------------------------------------------------------------------------
	CLLSEEEG
	
	functionality:

	Computes equating function based on continuized Log-linear cdf in Wang 	(2005). 
    	
	author: Tianyou Wang 2/22/2006.
	
	input:
        nparax      number of parameters for the new form
        paraxx		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the new form.
        nparay      number of parameters for the old form
        paraxy		a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the old form.
		nCatx		Number of discrete score categories for the new form
		Bx          design matrix from the new form loglinear model

	output:
 		SEE        a vector containing the standard error of equating
--------------------------------------------------------------------------*/
void CLLSEEEG(long npx, int nparax, double *parax, int nCatx, double *scoresx,
			  long npy, int nparay, double *paray,	int nCaty, double *scoresy, double *SEE)
{
	int i, j;
	double cdfx, pdfy;
	double **Sigmamx, **Sigmamy, *px, *py, **Sigmabetax, **Sigmabetay, **Bx, **By;
	double *eYbetax, *eYbetay;
	double integx, integx2, integy, integy2, minx, maxx, miny, maxy, Equatedx, integ1, integ2, integ3;

	Sigmamx = dmatrix(0, nCatx-1, 0, nCatx-1);
	Sigmamy = dmatrix(0, nCaty-1, 0, nCaty-1);
	Sigmabetax = dmatrix(0, nparax-1, 0, nparax-1);
	Sigmabetay = dmatrix(0, nparay-1, 0, nparay-1);
	Bx = dmatrix(0, nCatx-1, 0, nparax-1); 
	By = dmatrix(0, nCaty-1, 0, nparay-1); 
	px = dvector(0, nCatx-1);
	py = dvector(0, nCaty-1);
	eYbetax = dvector(0, nparax-1);
	eYbetay = dvector(0, nparay-1);

	for (i=0; i<nparax; i++)
	{
		for (j=0; j<nCatx; j++)
			if (j==0) 
				Bx[j][i] = 0;
			else 
			    Bx[j][i] = pow((double) j, (double) i);
	}
	Bx[0][0] = 1;   
	
	for (i=0; i<nparay; i++)
	{
		for (j=0; j<nCaty; j++)
			if (j==0) 
				By[j][i] = 0;
			else 
			    By[j][i] = pow((double) j, (double) i);
	}
	By[0][0] = 1;   
	
	for (i=0; i<nCatx; i++) {
		for (j=0; j<nCatx; j++) {
			Sigmamx[i][j] = -px[i] * px[j] / npx;
			if (i == j) 
				Sigmamx[i][j] += px[i];
		}
	}

	for (i=0; i<nCaty; i++) {
		for (j=0; j<nCaty; j++) {
			Sigmamy[i][j] = -py[i] * py[j] / npy;
			if (i == j) 
				Sigmamy[i][j] += py[i];
		}
	}

	PtAP0(nCatx, nparax, Bx, Sigmamx, Sigmabetax);
	PtAP0(nCaty, nparay, By, Sigmamy, Sigmabetay);

	/* ===== Before version 1.0 update =====  
	MatrixInverse0(nparax, Sigmabetax);
	MatrixInverse0(nparay, Sigmabetay);
	*/
	er_matrix_inverse(nparax, Sigmabetax);
	er_matrix_inverse(nparay, Sigmabetay);
	/* ===== End of version 1.0 update ===== */

	minx = scoresx[0] - .5;
	maxx = scoresx[nCatx -1 ] + .5;
	miny = scoresy[0] - .5;
	maxy = scoresy[nCaty -1 ] + .5;

	integx = GaussianQuadrature64(&ExpPolynomial, minx, maxx, nparax, parax);
	integx2 = integx * integx;
	integy = GaussianQuadrature64(&ExpPolynomial, miny, maxy, nparay, paray);
	integy2 = integy * integy;
	for (i=0; i<nCatx; i++) {
		cdfx = CLLEGCdf(minx, maxx, nparax, parax, scoresx[i], integx);
		Equatedx = CLLInverseCdf(miny, maxy, nparay, paray, cdfx, integy);
		pdfy = CLLEGPdf(miny, maxy, nparay, paray, Equatedx, integy);
		integ2 = GaussianQuadrature64(&ExpPolynomial, minx, scoresx[i], nparax, parax);
		for (j=0; j<nparax; j++) {
			integ1 = GaussianQuadrature64i(&ExpPolynomialxi, minx, scoresx[i], nparax, parax, j);
			integ3 = GaussianQuadrature64i(&ExpPolynomialxi, minx, maxx, nparax, parax, j);
			eYbetax[j] = integ1 * integx - integ2 * integ3;
			eYbetax[j] /= integx2;
			eYbetax[j] /= pdfy;
		}

		integ2 = GaussianQuadrature64(&ExpPolynomial, miny, Equatedx, nparay, paray);
		for (j=0; j<nparay; j++) {
			integ1 = GaussianQuadrature64i(&ExpPolynomialxi, miny, Equatedx, nparay, paray, j);
			integ3 = GaussianQuadrature64i(&ExpPolynomialxi, miny, maxy, nparay, paray, j);
			eYbetay[j] = integ2 * integ3 - integ1 * integy;
			eYbetay[j] /= integy2;
			eYbetay[j] /=  pdfy;
		}

		SEE[i] = vtAv0(nparax, eYbetax, Sigmabetax);
		SEE[i] += vtAv0(nparay, eYbetay, Sigmabetay);
		SEE[i] = sqrt(SEE[i]);
	}


}



/************************************************************************/

void Wrapper_RC(char design, char method, char smoothing,  
               struct USTATS *x, struct USTATS *y, 
               struct ULL_SMOOTH *ullx, struct ULL_SMOOTH *ully, int rep,
               struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for doing CLL equating with RG design
    and log-linear smoothing smoothing
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'R'(random groups)
    method = 'E'(equipercentile)
    smoothing = 'Z' (CLL 'smoothing')  
    *x = pointer to struct USTATS (new form)
    *y = pointer to struct USTATS (old form)
    *ullx = pointer to struct BB_SMOOTH (new form)
    *ully = pointer to struct BB_SMOOTH (old form)
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
      
  NOTE: If Wrapper_RC() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 

  Function calls other than C or NR utilities:                   
    CLLEquateEG()
    MomentsFromFD()  
                                                
  Tianyou Wang and Robert L. Brennan

  Date of last revision: 6/30/08       
*/
{ 
                          /* method name --- 10 characters; right justified */
  char *names[] ={"   y-equiv"};         
  double maxx, maxy, *scoresx, *scoresy;
  int i;
  double *parax, *paray;

  parax = dvector(0, ullx->c);
  paray = dvector(0, ully->c);
  parax[0] = ullx->ap;
  for(i=0;i<ullx->c;i++) parax[i+1] = ullx->Beta[i];
  paray[0] = ully->ap;
  for(i=0;i<ully->c;i++) paray[i+1] = ully->Beta[i];
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
  CLLEquateEG(ullx->min, maxx, ullx->c+1, parax,  
			   ully->min, maxy, ully->c+1, paray, 	
			   ullx->ns, scoresx, r->eraw[0]);
/*  CLLSEEEG(ullx->num_persons, ullx->c, ullx->Beta, ullx->ns, scoresx, 
	       ully->num_persons, ully->c, ully->Beta, ully->ns, scoresy, SEE); */
/* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

} 


/********************************************************************************/

void Print_RC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
{
/*
  print RC results (RG; equipercentile; log-linear smoothing)
  
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
    
  fprintf(fp,"CLL Equating with Random Groups Design\n"
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

void Wrapper_SC(char design, char method, char smoothing,  struct BSTATS *xy, 
               struct BLL_SMOOTH *bllxy, int rep, struct PDATA *inall, 
			   struct ERAW_RESULTS *r)  
/*
  Wrapper for doing CLL equating with SG design
  and log-linear smoothing.  
  
  NOTE: This is for the SG design in which x and y do not share any items in 
  common, which means that functionally this is the external anchor case.  
  The bivariate log-linear smoothing procedure needs to know this. So, when
  Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
  anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
  (with anchor set to 0) and Wrapper_SC() can still be used, but convergence 
  of the smoothing algorithm may be compromised because of dependencies between
  x and y.
    
  Assumes that equating puts raw scores for x on scale of y
  
  NOTE: This function is used (unaltered) for both actual equating and 
        equating done in Wrapper_Bootstrap().  Distinguishing between the
        two is the purpose of the variable rep

  Input
  
    design = 'S' (single group)
    method = 'E' (equipercentile)
    smoothing = 'Z' (CLL 'smoothing')  
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
      
  NOTE: If Wrapper_SC() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 
                                            
  Function calls other than C or NR utilities:
    CLLEquateSG()
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
  
  CLLEquateSG(bllxy, r->eraw[0]);

 /* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

  return;

}

/*******************************************************************************/

void Print_SC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print SC results (SG; equipercentile; log-linear smoothing)
  
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
    
  fprintf(fp,"CLL Equating with Single Group Design\n"
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

void Wrapper_BC(char design, char method, char smoothing,  double wtsx, 
				double wtsy, struct BSTATS *xy1, struct BSTATS *xy2, 
				struct BLL_SMOOTH *bllxy1, struct BLL_SMOOTH *bllxy2, 
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r)  
/*
  Wrapper for doing CLL equating for SG with counter-balance design 
  and log-linear smoothing.  
  
  NOTE: This is for the SG with counter balance design in which x and y 
  in both group 1 and group 2 do not share any items in 
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
  
    design = 'B' (single group)
    method = 'E' (equipercentile)
    smoothing = 'Z' (CLL 'smoothing')  
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
      
  NOTE: If Wrapper_BC() is called in a bootstrap loop,
        then in the calling function struct ERAW_RESULTS must
        be different from struct ERAW_RESULTS for the actual
        equating. 
                                            
  Function calls other than C or NR utilities:
    CLLEquateCB()
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
/*    strcpy(inall->xyfname,xy->fname);
    inall->xy = xy;
	inall->bllxy = bllxy; */
    inall->design = design;
    inall->method = method;
    inall->smoothing = smoothing;
	inall->anchor = 0;  /* implicitly, anchor is external for biv log-linear
						                        smoothing with the CB design */
    
    inall->nm = 1;                                                                    
    inall->names = cmatrix(0,inall->nm-1,0,11);    /* only one row/method, 0 */
    strcpy(inall->names[0],names[0]);
 
 /*   inall->min = xy->min1;  
    inall->max = xy->max1;
    inall->inc = xy->inc1;
    inall->fdx = xy->fd1;
    inall->n = xy->n; */
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
  
  CLLEquateCB(bllxy1, bllxy2, wtsx, wtsy, r->eraw[0]);
 /* get moments */               
   
  MomentsFromFD(inall->min,inall->max,inall->inc,r->eraw[0],inall->fdx,r->mts[0]);

  return;

}

/*******************************************************************************/

void Wrapper_CC(char design, char method, char smoothing, 
                double w1, int anchor, double rv1, double rv2,
                struct BSTATS *xv, struct BSTATS *yv, 
				struct BLL_SMOOTH *bllxv, struct BLL_SMOOTH *bllyv,
			    int rep, struct PDATA *inall, struct ERAW_RESULTS *r)              
/*
  Wrapper for CLL equating for CG design with log-linear smoothing. 
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

    method:  'E' = Frequency estimation (FE)
             'F' = Modified freq est (MFE) 
             'G' = FE +  MFE
             'C' = Chained
			 'H' = FE +  Chained
             'A' = FE + MFE + Chained
              
    smoothing = 'Z' (CLL 'smoothing') 

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
   
    NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
    
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
    
  NOTE: If Wrapper_CC() is called in a bootstrap loop, then in
        the calling function struct ERAW_RESULTS must be different
        from struct ERAW_RESULTS for the actual equating. 
                                            
  Function calls other than C or NR utilities:
    CLLEquateNEATPS()
    CLLEquateNEATChn()
    runerror()
                                               
  T. D. Wang

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
*/

                                           /* FE + BH-FE in positions 0 and 1*/
  if(method == 'E' || method == 'G' || method == 'A' || method == 'H')     
	i = CLLEquateNEATPS(bllxv, bllyv, inall->w1, r->eraw[0]);

  if(method == 'C') ptr = r->eraw[0];
  else if(method == 'A') ptr = r->eraw[2];
  else if(method == 'H') ptr = r->eraw[1];
  else ptr = NULL;

  if(ptr)
	CLLEquateNEATChn(bllxv, bllyv, ptr);
                        
/* get moments */

  for(i=0;i<=inall->nm-1;i++) 
    MomentsFromFD(xv->min1,xv->max1,xv->inc1,r->eraw[i],xv->fd1,r->mts[i]);
  
  return;
} 

/****************************************************************************/

void Print_CC(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r)
/*
  print CC results (CG; log-linear smoothed equipercentile equating)
  
  Input
    fp = file pointer for output
    tt[] = user supplied text identifier
    inall =  struct PDATA
    r = struct ERAW_RESULTS 

  Function calls other than C or NR utilities: 
    score() 
                                            
  T. D. Wang

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
  
  fprintf(fp,"CLL Equating with CINEG Design"
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
