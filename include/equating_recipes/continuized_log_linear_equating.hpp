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

#ifndef CONTINUIZED_LOG_LINEAR_EQUATING_HPP
#define CONTINUIZED_LOG_LINEAR_EQUATING_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Core>

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/log_linear_equating.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>

namespace EquatingRecipes {
  class ContinuizedLogLinearEquating {
  public:
  private:
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
    double cllEGPdf(const double& min,
                    const double& max,
                    const std::vector<double>& para,
                    const double& x,
                    const double& nc) {
      double pdf;

      pdf = expPolynomial(para, x);

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
    double cllEGCdf(const double& min,
                    const double& max,
                    const std::vector<double>& para,
                    const double& x,
                    const double& nc) {
      double cdf;

      cdf = gaussianQuadrature64(min, x, para);

      cdf /= nc;

      return cdf;
    }

    double gaussianQuadrature(const double& a,
                              const double& b,
                              const std::vector<double>& para,
                              const std::vector<double>& x,
                              const std::vector<double>& w) {
      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s = 0.0;

      for (size_t j = 0; j < x.size(); j++) {
        double dx = xr * x[j];
        s += w[j] * (expPolynomial(para, (xm + dx)) + expPolynomial(para, (xm - dx)));
      }

      s *= xr;

      return s;
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
    double gaussianQuadrature16(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.09501250983764, 0.28160355077926, 0.45801677765723,
                             0.61787624440264, 0.755404408355, 0.86563120238783,
                             0.94457502307323, 0.98940093499165};
      std::vector<double> w {0.18945061045507, 0.18260341504492, 0.169156519395,
                             0.14959598881658, 0.12462897125553, 0.09515851168249,
                             0.06225352393865, 0.02715245941175};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
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
    double gaussianQuadrature32(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
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
    double gaussianQuadrature64(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
    }

    double GaussianQuadrature64i(const double& a,
                                 const double& b,
                                 const std::vector<double>& para,
                                 const size_t& i) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s = 0;

      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s += w[j] * (expPolynomialxi(para, (xm + dx), i) +
                     expPolynomialxi(para, (xm - dx), i));
      }

      s *= xr;

      return s;
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
    double expPolynomial(const std::vector<double>& para,
                         const double& x) {
      double s = 0;

      for (size_t j = 0; j < para.size(); j++) {
        s += para[j] * std::pow(x, static_cast<double>(j));
      }

      s = std::exp(s);

      return s;
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
    double expPolynomialxi(const std::vector<double>& para,
                           const double& x,
                           const size_t& i) {
      double s = 0;

      for (size_t j = 0; j < para.size(); j++)
        s += para[j] * std::pow(x, j);

      s = exp(s);
      s *= pow(x, static_cast<double>(i));

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
    void CalcLLContinuMoments(std::function<double(const double&, const double&, const std::vector<double>&, const double&)> pdf,
                              const double& a,
                              const double& b,
                              const std::vector<double>& para,
                              std::vector<double>& moments) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s1 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s1 += w[j] * (pdf(a, b, para, (xm + dx)) * (xm + dx) +
                      pdf(a, b, para, (xm - dx)) * (xm - dx));
      }
      s1 *= xr;

      double s2 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s2 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 2) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 2));
      }
      s2 *= xr;

      double s3 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s3 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 3) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 3));
      }
      s3 *= xr / std::pow(s2, 1.5);

      double s4 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s4 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 4) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 4));
      }
      s4 *= xr / std::pow(s2, 2);

      moments = {s1,
                 std::sqrt(s2),
                 s3,
                 s4};
    }

    /*--------------------------------------------------------------------------
      CLLEquateEG
      
      functionality:

      Computes equating function based on continuized Log-linear cdf in Wang (2005). 
          
      author: Tianyou Wang 1/5/2005.
      
      input:
        minx        lower limit of the distribution for the new form
        maxx        upper limit of the distribution for the new form
        nparax      number of parameters for the new form
        paraxx		  a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the new form.
        miny        lower limit of the distribution for the old form
        maxy        upper limit of the distribution for the old form
        nparay      number of parameters for the old form
        paraxy		  a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the old form.
        nDistCatx   Number of discrete score categories for the new form
        scoresx     vector containing the discrete scores for the new form

      output:
        Equatedx   a vector containing the equated score 
    --------------------------------------------------------------------------*/
    std::vector<double> cllEquateEG(const double& minx,
                                    const double& maxx,
                                    const std::vector<double>& parax,
                                    const double& miny,
                                    const double& maxy,
                                    const std::vector<double>& paray,
                                    const std::vector<double>& scoresx) {
      double ncx = gaussianQuadrature64(minx, maxx, parax);
      double ncy = gaussianQuadrature64(miny, maxy, paray);
      std::vector<double> equatedx;

      std::for_each(scoresx.begin(),
                    scoresx.end(),
                    [&](const double& x) {
                      double cdfx = cllEGCdf(minx, maxx, parax, x, ncx);
                      double eqautedScore = cllInverseCdf(miny, maxy, paray, cdfx, ncy);

                      equatedx.push_back(eqautedScore);
                    });
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
    double cllInverseCdf(const double& min,
                         const double& max,
                         const std::vector<double>& para,
                         const double& cdf,
                         const double& nc) {
      double eps = .000001;

      double lb = min;
      double ub = max;
      double cdfl = cllEGCdf(min, max, para, lb, nc);
      double cdfu = cllEGCdf(min, max, para, ub, nc);
      double half;

      if (cdf < cdfl) {
        half = min;
      } else if (cdf > cdfu) {
        half = max;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = 0.5 * (lb + ub);
          double cdfhalf = cllEGCdf(min, max, para, half, nc);
          double absdif = std::abs(cdf - cdfhalf);

          if (absdif < eps) {
            break;
          } else if (cdfhalf < cdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
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
    double cllBivPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                     const double& x,
                     const double& y) {
      // static double integ;
      // double minx, maxx, miny, maxy, pdf;

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;

      double integ = bivGaussianQuadrature64(bivar, minx, maxx, miny, maxy);
      double pdf = bivExpPolynomial(bivar, x, y) / integ;

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
    double CLLBivCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                     const double& x,
                     const double& y,
                     const double& nc) {
      // double minx, maxx, miny, maxy, cdf;

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;

      double cdf = bivGaussianQuadrature32(bivar, minx, x, miny, y);

      cdf /= nc; /* bivar->num_persons; */

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
    double cllMargYPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing bivar,
                       const double& y,
                       const double& nc) {
      // int j;
      // double minx, maxx, miny, maxy;
      // double xr, xm, dx, s;
      std::vector<double> x {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;
      double xm = 0.5 * (maxx + minx);
      double xr = 0.5 * (maxx - minx);
      double s = 0;
      for (size_t j = 0; j < 16; j++) {
        double dx = xr * x[j];
        s += w[j] * (bivExpPolynomial(bivar, xm + dx, y) + bivExpPolynomial(bivar, xm - dx, y));
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
    double cllMargXPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                       const double& x,
                       const double& nc) {
      // int j;
      // double minx, maxx, miny, maxy;
      // double yr, ym, dy, s;
      std::vector<double> y {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;
      double ym = 0.5 * (maxy + miny);
      double yr = 0.5 * (maxy - miny);
      double s = 0;
      for (size_t j = 0; j < 16; j++) {
        double dy = yr * y[j];
        s += w[j] * (bivExpPolynomial(bivar, x, ym + dy) + bivExpPolynomial(bivar, x, ym - dy));
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
    double bivExpPolynomial(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                            const double& x,
                            const double& y) {
      // int j, nx, ny, nxbyy;
      // double s = 0;

      double xValue = x;
      double yValue = y;

      /* when internal anchor, f(x,v) = f(u,v) for u = x - v and u is within valid range */
      if (bivar.isInternalAnchor) {
        xValue -= yValue;

        if (xValue < bivar.minimumRawScoreU) {
          return 0;
        } else if (xValue > (bivar.minimumRawScoreU + static_cast<double>(bivar.numberOfScoresU - 1) * bivar.scoreIncrementU)) {
          return 0;
        }
      }

      size_t nx = bivar.numberOfDegreesOfSmoothingU;
      size_t ny = bivar.numberOfDegreesOfSmoothingV;
      size_t nxbyy = bivar.numberOfCrossProductMoments;

      double s = bivar.cllNormalizingConstant;

      for (size_t j = 1; j <= nx; j++) {
        s += bivar.betaCoefficients(j - 1) * std::pow(xValue, static_cast<double>(j));
      }

      for (j = 1; j <= ny; j++) {
        s += bivar.betaCoefficients(nx + j - 1) * std::pow(yValue, static_cast<double>(j));
      }

      for (j = 1; j <= nxbyy; j++) {
        s += bivar.betaCoefficients(nx + ny + j - 1) *
             std::pow(x, bivar.crossProductMoments(j - 1, 0)) *
             std::pow(y, bivar.crossProductMoments(j - 1, 1));
      }

      s = std::exp(s);

      return s;
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
  };
} // namespace EquatingRecipes

#endif