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

#ifndef KERNEL_EQUATING_HPP
#define KERNEL_EQUATING_HPP

#include <cmath>
#include <Eigen/Core>

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>

#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/log_linear_equating.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class KernelEquating {
  public:
    void Wrapper_RK(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const EquatingRecipes::Structures::UnivariateStatistics& x,
                    const EquatingRecipes::Structures::UnivariateStatistics& y,
                    const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& ullx,
                    const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& ully,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {}

    void Wrapper_SK(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const EquatingRecipes::Structures::BivariateStatistics& xy,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxv,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {}

    void Wrapper_CK(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const double& w1,
                    const bool& isInternalAnchor,
                    const double& reliabilityCommonItemsPopulation1,
                    const double& reliabilityCommonItemsPopulation2,
                    const EquatingRecipes::Structures::BivariateStatistics& xv,
                    const EquatingRecipes::Structures::BivariateStatistics& yv,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxv,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllyv,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {}

  private:
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
    double KernelContinuPdf(const size_t& ncat,
                            const Eigen::VectorXd& scores,
                            const Eigen::VectorXd& fd,
                            const double& hx,
                            const double& x) {
      double zscore; /*the standardized normal random variable */

      // double ax, sumsq = 0.0; /*ax, and second moment */
      // double temp1;           /*temporary variable */
      double pdf = 0;

      /*the mean and standard deviation of the unsmoothed distribution */
      double mu = scores.cwiseProduct(fd).sum();
      double sigma = (scores.cwiseProduct(scores)).cwiseProduct(fd).sum() - std::pow(mu, 2);

      double ax = std::sqrt(sigma / (sigma + hx * hx));

      for (size_t i = 0; i < ncat; i++) {
        double zscore = (x - ax * scores(i) - (1.0 - ax) * mu) / ax / hx;
        double temp1 = StdNormalPdf(zscore);
        pdf += fd(i) * temp1 / ax / hx;
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
    double KernelContinuCdf(const size_t& ncat,
                            const Eigen::VectorXd& scores,
                            const Eigen::VectorXd& fd,
                            const double& hx,
                            const double& x) {
      double mu = scores.cwiseProduct(fd).sum();
      double sigma = (scores.cwiseProduct(scores).cwiseProduct(fd)).sum() - std::pow(mu, 2);

      double ax = sqrt(sigma / (sigma + hx * hx));

      double cdf = 0;

      for (size_t i = 0; i < ncat; i++) {
        double zscore = (x - ax * scores(i) - (1.0 - ax) * mu) / ax / hx;
        double temp1 = StdNormalCdf(zscore);
        cdf += fd(i) * temp1;
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
    double StdNormalCdf(const double& x) {
      double sum = 0.0;
      double t, z;

      t = 1 / (1 + .2316419 * std::abs(x));
      z = 0.398942280401433 * std::exp(-0.5 * x * x);
      sum = 1 - z * t * (.319381530 + t * (-.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))));

      if (x >= 0) {
        return sum;
      } else {
        return 1 - sum;
      }
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
    double StdNormalPdf(const double& x) {
      return 0.398942280401433 * std::exp(-0.5 * x * x);
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
    double NormalPdf(const Eigen::VectorXd& para,
                     const double& x) {
      double mu;
      double sigma;
      if (para.size() == 2) {
        mu = para(0);
        sigma = para(1);
      } else {
        mu = 0.0;
        sigma = 1.0;
      }

      double pdf = 0.398942280401433 * std::exp(-0.5 * std::pow((x - mu) / sigma, 2)) / sigma;

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
    // void CalcKernelMoments(const double& a,
    //                        const double& b,
    //                        const size_t& ncat,
    //                        const Eigen::VectorXd& score,
    //                        const Eigen::VectorXd& fd,
    //                        const double& hx,
    //                        Eigen::VectorXd& moments) {
    //   size_t j;
    //   double xr, xm, dx, ss, s1, s2, s3, s4;
    //   std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
    //                          0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
    //                          0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
    //                          0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
    //                          0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
    //                          0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
    //                          0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
    //                          0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
    //                          0.99930504173577};
    //   std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
    //                          0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
    //                          0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
    //                          0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
    //                          0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
    //                          0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
    //                          0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
    //                          0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
    //                          0.00178328072170};

    //   xm = 0.5 * (b + a);
    //   xr = 0.5 * (b - a);

    //   ss = 0;
    //   for (j = 0; j < 32; j++) {
    //     dx = xr * x[j];
    //     ss += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) +
    //                   pdf(ncat, score, fd, hx, (xm - dx)));
    //   }
    //   ss *= xr;

    //   s1 = 0;
    //   for (j = 0; j < 32; j++) {
    //     dx = xr * x[j];
    //     s1 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * (xm + dx) +
    //                   pdf(ncat, score, fd, hx, (xm - dx)) / ss * (xm - dx));
    //   }
    //   s1 *= xr;

    //   s2 = 0;
    //   for (j = 0; j < 32; j++) {
    //     dx = xr * x[j];
    //     s2 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 2) +
    //                   pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 2));
    //   }
    //   s2 *= xr;

    //   s3 = 0;
    //   for (j = 0; j < 32; j++) {
    //     dx = xr * x[j];
    //     s3 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 3) +
    //                   pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 3));
    //   }
    //   s3 *= xr / pow(s2, 1.5);

    //   s4 = 0;
    //   for (j = 0; j < 32; j++) {
    //     dx = xr * x[j];
    //     s4 += w[j] * (pdf(ncat, score, fd, hx, (xm + dx)) / ss * pow(xm + dx - s1, 4) +
    //                   pdf(ncat, score, fd, hx, (xm - dx)) / ss * pow(xm - dx - s1, 4));
    //   }
    //   s4 *= xr / pow(s2, 2);

    //   moments[0] = s1;
    //   moments[1] = sqrt(s2);
    //   moments[2] = s3;
    //   moments[3] = s4;
    // }

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
    double Pen1(const size_t& ncat,
                const Eigen::VectorXd& scores,
                const Eigen::VectorXd& fd,
                const double& hx) {
      double sum = 0;

      for (size_t i = 0; i < ncat; i++) {
        double pdf = KernelContinuPdf(ncat, scores, fd, hx, scores(i));
        sum += (pdf - fd(i)) * (pdf - fd(i));
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
    double Pen2(const size_t& ncat,
                const Eigen::VectorXd& scores,
                const Eigen::VectorXd& fd,
                const double& hx) {}

    double KernelPdfDerivative(const size_t& ncat,
                               const Eigen::VectorXd& scores,
                               const Eigen::VectorXd& fd,
                               const double& hx,
                               const double& x) {
 
                               }

    double Pen(const size_t& ncat,
               const Eigen::VectorXd& scores,
               const Eigen::VectorXd& fd,
               const double& hx,
               const double& K) {}

    double Optimalh(const size_t& ncat,
                    const Eigen::VectorXd& scores,
                    const Eigen::VectorXd& fd,
                    const double& K) {}

    double KernelInverseCdf(const size_t& ncat,
                            const Eigen::VectorXd& scores,
                            const Eigen::VectorXd& fd,
                            const double& h,
                            const double& cdf) {}

    void KernelEquate(const size_t& ncatx,
                      const Eigen::VectorXd& scoresx,
                      const Eigen::VectorXd& fdx,
                      const double& hx,
                      const size_t& ncaty,
                      const Eigen::VectorXd& scoresy,
                      const Eigen::VectorXd& fdy,
                      const double& hy,
                      const Eigen::VectorXd& Equatedx) {}

    void ComputeCmatrix(const size_t& ncat,
                        const size_t& degree,
                        const long& np,
                        const Eigen::VectorXd& B,
                        const Eigen::VectorXd& fd,
                        const Eigen::VectorXd& Cr) {}

    void ComputeCmatrixGen(const size_t& ncat,
                           const size_t& degree,
                           const long& np,
                           const Eigen::VectorXd& B,
                           const Eigen::VectorXd& fd,
                           const Eigen::VectorXd& Cr) {}

    void PartialFPartialr(const size_t& ncat,
                          const Eigen::VectorXd& scores,
                          const Eigen::VectorXd& fd,
                          const double& hx,
                          const Eigen::VectorXd& Fr,
                          const double& x) {}

    double FrCrSqNorm(const size_t& ncat,
                      const size_t& degree,
                      const Eigen::VectorXd& Fr,
                      const Eigen::VectorXd& Cr) {}

    void vPMN(const size_t& ncatx,
              const size_t& ncaty,
              const Eigen::VectorXd& bdist,
              const Eigen::VectorXd& vP,
              const Eigen::VectorXd& M,
              const Eigen::VectorXd& N) {}

    void vPT(const size_t& ncatx,
             const size_t& ncaty,
             const Eigen::VectorXd& bdist,
             const Eigen::VectorXd& vPP) {}

    void MatrixMultiVector(const size_t& nrow,
                           const size_t& ncol,
                           const Eigen::VectorXd& m,
                           const Eigen::VectorXd& v,
                           const Eigen::VectorXd& r) {}

    void MatrixMultiMatrix(const size_t& nrow1,
                           const size_t& ncol1,
                           const size_t& row2,
                           const size_t& ncol2,
                           const Eigen::VectorXd& m,
                           const Eigen::VectorXd& n,
                           const Eigen::VectorXd& r) {}

    void VectorMultiMatrix(const size_t& nrow,
                           const size_t& ncol,
                           const Eigen::VectorXd& v,
                           const Eigen::VectorXd& m,
                           const Eigen::VectorXd& r) {}

    double VectorNormSq(const size_t& ncat,
                        const Eigen::VectorXd& v) {}

    void KernelEquateSEERG(const size_t& ncatx,
                           const size_t& degreex,
                           long npx,
                           const Eigen::VectorXd& scoresx,
                           const Eigen::VectorXd& fdx,
                           const size_t& ncaty,
                           const size_t& degreey,
                           long npy,
                           const Eigen::VectorXd& scoresy,
                           const Eigen::VectorXd& fdy,
                           const Eigen::VectorXd& Equatedx,
                           const Eigen::VectorXd& SEE) {}

    void KernelEquateSEESG(struct BLL_SMOOTH* bivar,
                           const Eigen::VectorXd& Equatedx,
                           const Eigen::VectorXd& SEE) {}

    void KernelEquateSG(struct BLL_SMOOTH* bivar,
                        const Eigen::VectorXd& Equatedx) {}

    void KernelEquateSEECB(struct BLL_SMOOTH* bivar1,
                           struct BLL_SMOOTH* bivar2,
                           const Eigen::VectorXd& wts,
                           const Eigen::VectorXd& Equatedx,
                           const Eigen::VectorXd& SEE) {}

    void KernelEquateSEECB2(struct BLL_SMOOTH* bivar1,
                            struct BLL_SMOOTH* bivar2,
                            const Eigen::VectorXd& wts,
                            const Eigen::VectorXd& Equatedx,
                            const Eigen::VectorXd& SEE) {}

    void KernelEquateNEATPS(struct BLL_SMOOTH* bivar1,
                            struct BLL_SMOOTH* bivar2,
                            const double& wts,
                            const Eigen::VectorXd& Equatedx) {}

    void KernelBootStrapNEATPS(struct BLL_SMOOTH* bivar1,
                               struct BLL_SMOOTH* bivar2,
                               const double& wts,
                               const Eigen::VectorXd& BootStrapSEE,
                               const Eigen::VectorXd& MeanEq) {}

    void KernelBootStrapNEATPS2(struct BLL_SMOOTH* bivar1,
                                struct BLL_SMOOTH* bivar2,
                                const double& wts,
                                const Eigen::VectorXd& BootStrapSEE,
                                const Eigen::VectorXd& MeanEq) {}

    void KernelEquateSEENEATPS(struct BLL_SMOOTH* bivar1,
                               struct BLL_SMOOTH* bivar2,
                               const double& wts,
                               const Eigen::VectorXd& Equatedx,
                               const Eigen::VectorXd& SEE) {}

    void KernelEquateSEENEATChn(struct BLL_SMOOTH* bivar1,
                                struct BLL_SMOOTH* bivar2,
                                const Eigen::VectorXd& Equatedx,
                                const Eigen::VectorXd& SEE) {}

    void KernelEquateNEATChn(struct BLL_SMOOTH* bivar1,
                             struct BLL_SMOOTH* bivar2,
                             const Eigen::VectorXd& Equatedx) {}

    void KernelBootStrapNEATChn(struct BLL_SMOOTH* bivar1,
                                struct BLL_SMOOTH* bivar2,
                                const double& wts,
                                const Eigen::VectorXd& BootStrapSEE,
                                const Eigen::VectorXd& MeanEq) {}

    void KernelBootStrapNEATChn2(struct BLL_SMOOTH* bivar1,
                                 struct BLL_SMOOTH* bivar2,
                                 const double& wts,
                                 const Eigen::VectorXd& BootStrapSEE,
                                 const Eigen::VectorXd& MeanEq) {}
  } {}
} // namespace EquatingRecipes

#endif