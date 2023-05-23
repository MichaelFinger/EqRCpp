#ifndef BFGS_OPTIMIZER
#define BFGS_OPTIMIZER

#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <vector>
#include <string>

#include <Eigen/Core>

#include <equating_recipes/optimization_function.hpp>
#include <equating_recipes/structures/optimization_results.hpp>

#include <iostream>

namespace EquatingRecipes {
  class BFGSOptimizer {
  public:
    EquatingRecipes::Structures::OptimizationResults optimize(std::vector<double>& x,
                                                              std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction,
                                                              const std::optional<size_t>& maximumNumberOfIterations,
                                                              const std::optional<double>& maximumAbsoluteChangeInFunctionValue,
                                                              const std::optional<double>& maximumRelativeChangeInFunctionValue,
                                                              const std::optional<double>& maximumAbsoluteChangeInParameterValues,
                                                              const std::optional<double>& maximumRelativeChangeInParameterValues) {
      EquatingRecipes::Structures::OptimizationResults results;

      std::for_each(x.begin(),
                    x.end(),
                    [&](const double& value) {
                      results.parameterStartingValues.push_back(value);
                    });

      double minf;
      size_t numberOfIterations;

      Eigen::VectorXd xVect(x.size());
      for (size_t index = 0; index < x.size(); index++) {
        xVect(index) = x[index];
      }

      er_dfpmin(xVect,
                maximumRelativeChangeInParameterValues.value_or(0.0000000001),
                numberOfIterations,
                maximumNumberOfIterations.value_or(250),
                minf,
                optimizationFunction);

      results.functionValue = minf;

      for (size_t index = 0; index < x.size(); index++) {
        x[index] = xVect(index);
      }

      results.parameterEstimates = x;

      std::vector<double> grad;
      grad.resize(x.size());
      optimizationFunction->operator()(results.parameterEstimates,
                                       grad);

      results.gradientVector = grad;

      results.numberOfIterations = numberOfIterations;
      // results.maximumAbsoluteChangeInFunctionValue = opt.get_ftol_abs();
      // results.maximumRelativeChangeInFunctionValue = opt.get_ftol_rel();
      // results.maximumAbsoluteChangeInParameterValues = *iter;
      // results.maximumRelativeChangeInParameterValues = opt.get_xtol_rel();

      // results.resultCode = getOptimizationResultCode(result);
      // } catch (const std::exception& e) {
      //   std::cerr << e.what() << '\n';
      // }

      return results;
    }

  private:
    /*
      Purpose:
          Find lamda so that f(xold + lamda*sk) has decreased sufficiently. To check if
        the decrease is sufficient, we use the following condition

      Input:    
        xold[1..n]      n-dimensional point
        n               dimension of the array xold, gk, sk, xnew, f
      fvalue_old      f(xold)
      gk[1..n]
      sk[1..n]    the Newton direction
      xnew[1..n]      new point along the direction sk from xold
                      where the funtion ftn_ptr has decreased sufficiently.
      maxstep         maximum step length used to keep the function outside regions 
                      where it is undefined or subject to overflow.
      ftn_ptr         pointer to an input function f( )  
        

      Output:     
      fvalue_new      new function value
      return value = 0 if normal exit
                    = 1 if xnew is too close to xold
      Precondition:     

      Author: Jaehoon Seol
      Date: August 20, 2009
      Version: 1.0
      References :  
      Comments: 
          For detailed explanation of the algorithm, refer to p384-p385, NR.
    */
    int er_lnsrch(const Eigen::VectorXd& xold,
                  const Eigen::VectorXd& gk,
                  Eigen::VectorXd& sk,
                  const double& maxstep,
                  std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction,
                  Eigen::VectorXd& xnew) {
      int rcode;  /* return code value                        */
                  /* rcode = 0   on normal exit               */
                  /*       = 1   xnew is too close to xold    */
      double a,   /* coefficient  of gk(lamda)                */
          b,      /* coefficient  of gk(lamda)                */
          lamda0, /* temp. lamda                              */
          lamda1,
          lamda2,
          lamda1p2,
          lamda2p2,
          dgx0,
          scale,
          fvalue,
          fvalue2,
          fvalue_old, /* old function value                       */
          temp1,      /* temp. used to compute a                  */
          temp2,      /* temp. used to compute b                  */
          temp,       /* temp. variable                           */
          nrm2 = 0.0, /* L2 norm of sk                            */
          test = 0.0;
      double error = 2.0 * std::pow(10.0, -6.0); /* tolerance level           */
      const double alpha = std::pow(10.0, -7.0); /* constant alpha in 9.7.7   */
      size_t numberOfParameters = xold.size();

      rcode = 0;                    /* normal exit               */
      nrm2 = std::sqrt(sk.dot(sk)); /* compute L2-norm of sk     */
      if (nrm2 > maxstep) {         /* normalize sk vector       */
        sk *= maxstep / nrm2;
      }

      dgx0 = gk.dot(sk);

      for (size_t i = 0; i < numberOfParameters; i++) {
        temp = (std::abs(xold[i]) > 1.0 ? fabs(xold[i]) : 1.0);
        temp = std::abs(sk[i]) / temp;
        if (temp > test) {
          test = temp;
        }
      }

      error /= test;

      lamda1 = 1.0;
      fvalue_old = optimizationFunction->evaluateFunction(xold);
      fvalue = fvalue_old + alpha * lamda1 * dgx0 + 1.0;

      while (fvalue > fvalue_old + alpha * lamda1 * dgx0) {
        xnew = lamda1 * sk + xold;

        fvalue = optimizationFunction->evaluateFunction(xnew);

        /*  xnew is too close to xold   */
        if (lamda1 < error) {
          xnew = xold;
          rcode = 1;
          break;
        }

        /* backtracking: compute lamda0 */
        if (lamda1 != 1.0) {
          /* compute a & b in Eq. 9.7.13 */
          temp1 = fvalue - fvalue_old - lamda1 * dgx0;
          temp2 = fvalue2 - fvalue_old - lamda2 * dgx0;
          scale = 1.0 / (lamda1 - lamda2);
          lamda1p2 = std::pow(lamda1, -2.0);
          lamda2p2 = std::pow(lamda2, -2.0);
          a = scale * (temp1 * lamda1p2 - temp2 * lamda2p2);
          b = scale * (temp2 * lamda1 * lamda2p2 - temp1 * lamda2 * lamda1p2);
          if (a == 0.0) { /* use Eq. 9.7.11       */
            lamda0 = -dgx0 / (2.0 * b);
          } else { /* use Eq. 9.7.14       */
            /* check b*b-3.0*a*gk'(0) >= 0.0 */
            temp = std::pow(b, 2.0) - 3.0 * a * dgx0;

            if (temp < 0.0) {
              std::string msg = "Equating Recipes error in Source: er_lnsrch, Error: Roundoff problem";
              throw std::runtime_error(msg);
            }

            lamda0 = -b + std::sqrt(temp);
            lamda0 /= (3.0 * a);
          }

          if (lamda0 > 0.5 * lamda1) {
            lamda0 = 0.5 * lamda1;
          }
        } else {
          lamda0 = -dgx0 / (2.0 * (fvalue - fvalue_old - dgx0));
        }

        fvalue2 = fvalue;
        lamda2 = lamda1;
        lamda1 = (lamda0 > 0.1 * lamda1 ? lamda0 : 0.1 * lamda1);
      }

      return rcode;
    }

    /*
      Purpose:
          This is an implementation of BFGS that solves the minimization problem for a real 
        function f : R^n -> R of n variable.

      Input:    
        xold[1..n]      n-dimensional starting point 
        n               dimension of the array startx
      error           tolerance level
      ftp_ptr         function pointer to f( )
      dftn_ptr        function pointer to D[f( )]
        

      Output:     
      xold[1..n]        xold is over-written by the final minimizer. 

      Precondition:     

      Author: Jaehoon Seol  
      Date: August 20, 2009
      Version: 1.0
      References :  
      Comments: 
    */
    void er_dfpmin(Eigen::VectorXd& xold,
                   const double& error,
                   size_t& numiter,
                   const size_t& maximumNumberOfIterations,
                   double& fvalue,
                   std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction) {
      int i;
      // num = 0; /* counts number of iterations      */
      double scale1,
          scale2,
          maxstep,
          nrm2;
      Eigen::VectorXd gk;    /* gk = D[f(x_(k+1))-f(x_k)]       */
      Eigen::VectorXd qk;    /* qk = g_(k+1)-g_k                */
      Eigen::VectorXd hddf;  /* hddf = hessian * qk             */
      Eigen::MatrixXd matxh; /* hessian matrix                  */
      Eigen::VectorXd xnew;  /* x_(k+1)                         */
      Eigen::VectorXd sk;    /* sk = Hk * gk, search direction  */
      Eigen::VectorXd pk;    /* pk = xnew-xold                  */

      Eigen::Index numberOfParameters = xold.size();

      pk.resize(numberOfParameters);
      qk.resize(numberOfParameters);
      gk.resize(numberOfParameters);
      sk.resize(numberOfParameters);
      hddf.resize(numberOfParameters);
      xnew.resize(numberOfParameters);
      matxh.resize(numberOfParameters, numberOfParameters);

      gk = optimizationFunction->evaluateGradient(xold);

      /* set hessian matrix by identity matrix */
      matxh.setIdentity(numberOfParameters, numberOfParameters);

      nrm2 = std::sqrt(xold.dot(xold));

      if (nrm2 > static_cast<double>(numberOfParameters)) {
        maxstep = 100.0 * nrm2;
      } else {
        maxstep = 100.0 * static_cast<double>(numberOfParameters);
      }

      // size_t maximumNumberOfIterations = 250;
      size_t num;
      for (num = 0; num <= maximumNumberOfIterations; num++) {
        sk = -1.0 * matxh * gk;

        er_lnsrch(xold, gk, sk, maxstep, optimizationFunction, xnew); /* determine x(k+1)*/

        sk = -1.0 * xold + xnew; // y = alpha*x +y

        fvalue = optimizationFunction->evaluateFunction(xold);
        xold = xnew;
        qk = -1.0 * gk;

        gk = optimizationFunction->evaluateGradient(xold);

        qk += gk;

        /* stopping condition */
        std::cout << "Iteration = " << num << ", "
                  << "sk.dot(sk) = " << std::sqrt(sk.dot(sk)) << ", "
                  << "gk.dot(gk) = " << std::sqrt(gk.dot(gk)) << "\n";

        if (std::sqrt(sk.dot(sk)) < std::pow(10.0, -8.0) ||
            std::sqrt(gk.dot(gk)) < error) {
          break;
        }

        hddf = matxh * qk; /* update hessian matrix   */
        scale1 = qk.dot(sk);
        scale2 = qk.dot(hddf);

        if (scale1 * scale1 > 3.0e-8 * qk.dot(qk) * sk.dot(sk)) {
          qk = (1.0 / scale1) * sk;
          qk = (-1.0 / scale2) * hddf + qk;

          // matx + scale * vect* vect^t
          matxh += (1.0 / scale1) * sk * (sk.transpose());
          matxh += (-1.0 / scale2) * hddf * (hddf.transpose());
          matxh += scale2 * qk * qk.transpose();
        }
      }

      numiter = num;
    }
  };
} // namespace EquatingRecipes

#endif