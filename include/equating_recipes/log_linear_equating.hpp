/* 	
  LogLinear.c  code for log-linear equating

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

  File contains functions used for:
      (a) univariate log-linear fitting under the multinomial model,
          which populates struct ULL_SMOOTH 
      (b) bivariate log-linear fitting under the multinomial model,
          which populates struct BLL_SMOOTH 

  Code follows procedures (and usually the notation) in Holland and
  Thayer (1987), which is sometimes abbreviated H&T in comments.

  A somewhat unique feature of these functions is that the user 
  can select from among various convergence criteria, including
  criteria based on moments for the actual raw scores (e.g., 
  mean, sd, skew, kurt, etc for the scores obtained using
  score() in ERutilities.c).  
*/

#ifndef LOG_LINEAR_EQUATING_HPP
#define LOG_LINEAR_EQUATING_HPP

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/score_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class LogLinearEquating {
  public:
  private:
    /*
      Create design matrix.  Code and variable names are for a
      bivariate u*v distribution, where u is rows and v is columns.
      Note that x = u + v = total score, with v = common-item score 
      and u = non-common-item score.  We never create a design matrix
      with x and v where v is internal to x, because there would be
      a great deal of collinearity. 
      
      Input
        nsu = number of score categories for u
        minu = minimum score for u
        incu = increment for u
        nsv = number of score categories for v
        minv = minimum score for v
        incv = increment for v
        cu = number of degrees of smoothing for u
        cv = number of degrees of smoothing for v
        cuv = number of cross-product moments
        cpm[cuv-1][2] = zero-offset matrix designating cross-
                        product moments.  Example: let cuv = 3,
                        and the desired cross-product moments be 
                        (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                        Then cpm[0] = 1,1; cpm[1] = 1,2; and
                        cpm[2] = 2,1.  
        scale: 0 --> no scaling; 
              1 --> scale such that each column of B has
                    sum (elements) = 0 and sum (elements^2) = 1
                    
      Output
        B_raw[][] = zero-offset design matrix for raw scores; 
                    #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                    space already allocated for B_raw 
        B[][]     = zero-offset design matrix used for solution;
                    involves scaling if scale==1 
                    #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                    space already allocated for B 
        
      NOTE: For univariate smoothing, set 
            nsv=0, cv=0, cuv=0, cpm = NULL.  In this case, 
            obviously, u plays a generic role 
            (i.e., any single variable such as x or y)
            
      NOTE: Usually for bivariate smoothing in equating, set cuv=1
            and cpm[0] = 1,1. Otherwise, the conditional 
            distributions are not necessarily stochastically ordered
            (see Rosenbaum & Thayer, 1987, p. 46).

      NOTE: Using minu!=0 and incu!=1 changes both the B and 
            B_raw matrices; similarly for minv and incv.  
            If convergence not achieved, it may be wise to set
            minu = 0 and incu = 1 (same for minv and incv)

      Function calls other than C or NR utilities:  
        score()
      runerror()

      R. L. Brennan

      Date of last revision: 6/30/08

    */
    void designMatrix(const size_t& numberOfScoresU,
                      const double& minimumScoreU,
                      const double& scoreIncrementU,
                      const size_t& numberOfScoresV,
                      const double& minimumScoreV,
                      const double& scoreIncrementV,
                      const size_t& numberOfDegreesOfSmoothingU,
                      const size_t& numberOfDegreesOfSmoothingV,
                      const size_t& numberOfCrossProductMoments,
                      const Eigen::MatrixXi& crossProductMomentDesignations,
                      const bool& useStandardizedScale,
                      Eigen::MatrixXd& rawScoreDesignMatrix,
                      Eigen::MatrixXd& solutionDesignMatrix) {
      size_t numberOfDesignMatrixColumns = numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV + numberOfCrossProductMoments; /* # columns in design matrix */
      size_t numberOfDesigmMatrixRows = (numberOfScoresV > 0) ? numberOfScoresU * numberOfScoresV : numberOfScoresU;                /* # rows in design mat */

      /* univariate LL smoothing */
      if (numberOfScoresV == 0) {
        for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          double scoreU = EquatingRecipes::Utilities::getScore(scoreLocationU, minimumScoreU, scoreIncrementU);

          for (size_t degreesOfSmoothingUIndex = 1; degreesOfSmoothingUIndex <= numberOfDegreesOfSmoothingU; degreesOfSmoothingUIndex++) {
            rawScoreDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1) = std::pow(scoreU, static_cast<double>(degreesOfSmoothingUIndex));
            solutionDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1) = rawScoreDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1);
          }
        }
      } else if (numberOfDegreesOfSmoothingV == 0) {
        throw std::runtime_error("Number of degrees of smoothing for V or number of score categories for V is misspecified.");
      } else {
        /* bivariate LL smoothing */

        /* u polynomials */
        size_t cell = 0;
        for (size_t scoreLocationU; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          double scoreU = EquatingRecipes::Utilities::getScore(scoreLocationU, minimumScoreU, scoreIncrementU);

          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            for (size_t degreesOfSmoothingUIndex = 1; degreesOfSmoothingUIndex <= numberOfDegreesOfSmoothingU; degreesOfSmoothingUIndex++) {
              rawScoreDesignMatrix(cell, degreesOfSmoothingUIndex - 1) = std::pow(scoreU, static_cast<double>(degreesOfSmoothingUIndex));
              cell++;
            }
          }

          /* v polynomials */
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            double scoreV = EquatingRecipes::Utilities::getScore(scoreLocationV,
                                                                 minimumScoreV,
                                                                 scoreIncrementV);

            for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
              for (size_t degreesOfSmoothingVIndex = 1; degreesOfSmoothingVIndex <= numberOfDegreesOfSmoothingV; degreesOfSmoothingVIndex++) {
                rawScoreDesignMatrix(scoreLocationV + scoreLocationU * numberOfScoresV,
                                     numberOfDegreesOfSmoothingU + scoreLocationV - 1) = std::pow(scoreV, static_cast<double>(degreesOfSmoothingVIndex));
              }
            }
          }

          if (numberOfDegreesOfSmoothingV > 0 && (crossProductMomentDesignations.rows() == 0 || crossProductMomentDesignations.cols() == 0)) {
            throw std::runtime_error("Number of degrees of smoothing for V or cross-product moment designations misspecified.");
          }
        }

        /* cross-product polynomials */
        cell = 0;

        for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            for (size_t degreesOfSmoothingVIndex = 0; degreesOfSmoothingVIndex < numberOfDegreesOfSmoothingV; degreesOfSmoothingVIndex++) {
              rawScoreDesignMatrix(cell, numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV + degreesOfSmoothingVIndex) = rawScoreDesignMatrix(cell, crossProductMomentDesignations(degreesOfSmoothingVIndex, 0) - 1) *
                                                                                                                                 rawScoreDesignMatrix(cell, numberOfDegreesOfSmoothingU + crossProductMomentDesignations(degreesOfSmoothingVIndex, 1) - 1);

              cell++;
            }
          }
        }

        /* copy B_raw to B */
        for (size_t rowIndex = 0; rowIndex < numberOfDesigmMatrixRows; rowIndex++) {
          for (size_t columnIndex = 0; columnIndex < numberOfDesignMatrixColumns; columnIndex++) {
            solutionDesignMatrix(rowIndex, columnIndex) = rawScoreDesignMatrix(rowIndex, columnIndex);
          }
        }
      }

      /* scaling such that for each column of B,
        sum (elements) = 0 and sum (elements^2) = 1 */

      if (useStandardizedScale) {
        Eigen::VectorXd columnMeans = rawScoreDesignMatrix.colwise().mean();
        Eigen::VectorXd columnSumSquaredDeviations = (rawScoreDesignMatrix.cwiseProduct(rawScoreDesignMatrix)).colwise().sum() -
                                                     static_cast<double>(numberOfDesignMatrixColumns) * columnMeans.cwiseProduct(columnMeans);

        for (size_t columnIndex = 0; columnIndex < numberOfDesignMatrixColumns; columnIndex++) {
          Eigen::VectorXd vectorOfMean = Eigen::VectorXd::Constant(numberOfDesigmMatrixRows, columnMeans(columnIndex));

          rawScoreDesignMatrix.col(columnIndex) = (rawScoreDesignMatrix.col(columnIndex) - vectorOfMean) / std::sqrt(columnSumSquaredDeviations(columnIndex));
        }
      }
    }

    /*
      Get nct[] from bfd[][]

      Convert bivariate fd (xv->bfd[][]) to a vector nct[] with row j
      elements followed by row j+1 elements.  Note that for an 
      internal anchor the rows are x = u + v and the cols are v,
      which means that there are structural zeros, and we want nct[]
      to contain only the u*v elements.  See example below in which
      - indicates a structural 0 and the within-matrix numbers are 
      the cell locations in nct[].

                    v                                   v
                0  1  2  3                          0  1  2  3
          ---------------                      --------------
          0 |  0  -  -  -                      0|  0  1  2  3
          1 |  4  1  -  -    nsx = 9           1|  4  5  6  7
          2 |  8  5  2  -    nsv = 4 -->     u 2|  8  9 10 11  
          3 | 12  9  6  3    nsu = 6           3| 12 13 14 15
        x 4 | 16 13 10  7                      4| 16 17 18 19
          5 | 20 17 14 11                      5| 20 21 22 23
          6 |  - 21 18 15
          7 |  -  - 22 19
          8 |  -  -  - 23
    

      Input
        anchor : 0 --> external; 1 --> internal
        nsx = number of score categories for total scores (x)
        nsv = number of score categories for common-item scores (v)
        bfd[][] = bivariate freq dist for x and v

      Output
        nct[] = vector version of bfd[][], where nct[] is "collaped",
                as discussed above, if anchor is internal.
                Assumes space already allocated for nct[]
    
      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08      
    */
    void getNctBfd(const bool& isInternalAnchor,
                   const size_t& numberOfScoresX,
                   const size_t& numberOfScoresV,
                   const Eigen::MatrixXd& bivariateFreqDist,
                   Eigen::VectorXd& nct) {
      if (!isInternalAnchor) {
        /* external anchor */

        size_t cells = 0;
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            nct(cells++) = bivariateFreqDist(scoreLocationX, scoreLocationV);
          }
        }
      } else {
        /* internal anchor */

        size_t numberOfScoresU = numberOfScoresX - numberOfScoresV + 1; /* number of score categories for non-common items */

        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              continue;
            } else {
              nct((scoreLocationX - scoreLocationV) * numberOfScoresV + scoreLocationV) = bivariateFreqDist(scoreLocationX, scoreLocationV);
            }
          }
        }
      }
    }

    /*
      Get bfd[][] from mct[]

      For an external anchor, directly convert the row major 
      vector mct[nsu*nsv] to bfd[nsu][nsv]. For an 
      internal anchor, convert mct[nsu*nsv] to
      bfd[nsu+nsv][nsv] be adding structural zeros.
      See comments for get_nct_bfd().

      Input
        anchor : 0 --> external; 1 --> internal
        nsx = number of score categories for total scores (x)
        nsv = number of score categories for common-item scores (v)
        mct[] = row-major vector for fitted bivariate frequencies
                for non-common items by common items

      Output
        bfd[][] = fitted bivariate frequencies for total scores by
                  common-item scores. That is mct[] is mapped into 
                  bfd[][] with structural zeros added if anchor is 
                  internal (see comments for get_nct_bfd()). 
                  Assumes space already allocated for bfd[][]
    
      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void getBfdMct(const bool& isInternalAnchor,
                   const size_t& numberOfScoresX,
                   const size_t& numberOfScoresV,
                   const Eigen::VectorXd& mct,
                   Eigen::MatrixXd& bivariateFreqDist) {
      if (!isInternalAnchor) {
        /* external anchor */

        size_t cells = 0;
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            bivariateFreqDist(scoreLocationX, scoreLocationV) = mct(cells++);
          }
        }
      } else {
        /* internal anchor */

        size_t numberOfScoresU = numberOfScoresX - numberOfScoresV + 1; /* number of score categories for non-common items */
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              bivariateFreqDist(scoreLocationX, scoreLocationV) = 0.0;
            } else {
              bivariateFreqDist(scoreLocationX, scoreLocationV) = mct((scoreLocationX - scoreLocationV) * numberOfScoresV + scoreLocationV);
            }
          }
        }
      }
    }

    /*
      Input:
        B[][] = design matrix (ns x nc)
        m[]   = fitted frequencies (ns x 1)

      Output:
        BtSmB[][] = Bt x Sm x B (nc x nc)
                  = minus the 2nd derivative of log-likelihood
                  (Eq. 22 and 32 in Holland & Thayer, 1987)

      Removed From Input:
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix
        N     = total of all frequencies

      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::MatrixXd get_BtSmB(const Eigen::MatrixXd& designMatrix,
                              const Eigen::VectorXd& fittedFrequencies) {
      Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(designMatrix.cols(),
                                                      designMatrix.cols());

      Eigen::VectorXd weightedSumDesignMatrixColumns =
          designMatrix.cwiseProduct(fittedFrequencies.replicate(designMatrix.cols(), 1)).colwise().sum();

      double fittedFrequencySum = fittedFrequencies.sum();

      for (size_t rowIndex = 0; rowIndex < designMatrix.cols(); rowIndex++) {
        for (size_t columnIndex = rowIndex; columnIndex < designMatrix.cols(); columnIndex++) {
          hessian(rowIndex, columnIndex) =
              (designMatrix.col(rowIndex) * designMatrix.col(columnIndex) * fittedFrequencies).sum();

          hessian(rowIndex, columnIndex) /= fittedFrequencySum;
          hessian(columnIndex, rowIndex) = hessian(rowIndex, columnIndex);
        }
      }

      return hessian;
    }

    /*
      Input:
        B[][] = design matrix (ns x nc)
        n[]   = actual frequencies (ns x 1)
        m[]   = fitted frequencies (ns x 1)

      Output:
        Btnm[] = 1st derivative of log-likelihood (nc x 1)
                (Eq. 19 and 33 in Holland & Thayer, 1987)

      Removed From Input:
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix

      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::VectorXd get_Btnm(const Eigen::MatrixXd& designMatrix,
                             const Eigen::VectorXd& observedFrequencies,
                             const Eigen::VectorXd& fittedFrequencies) {
      Eigen::VectorXd gradient(designMatrix.cols());

      for (size_t columnIndex = 0; columnIndex < designMatrix.cols(); columnIndex++) {
        gradient(columnIndex) = (designMatrix.col(columnIndex) * (observedFrequencies - fittedFrequencies)).sum();
      }

      return gradient;
    }

    /*
      Input:
        B[][] = zero-offset design matrix (ns x nc) 
        n[]   = zero-offset actual frequencies (ns x 1)
        
        fp    = output file pointer for debugging
                (NULL --> no output)

      Output:
        Beta0[] = zero-offset initial values of Beta (nc x 1);
                  space already allocated;
                  based on Eq 49 (and next line) of Holland and
                    Thayer (1987) -- abbreviated H&T below

      Removed From Input:
        N     = total of frequencies
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix

      Function calls other than C or NR utilities: 
        Print_vector()
        Print_matrix() 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::VectorXd getBeta0(const Eigen::MatrixXd& designMatrix,
                  const Eigen::VectorXd& observedFrequencies,
                  const bool& debug) {
      size_t numberOfRows = designMatrix.rows();
      size_t numberOfColumns = designMatrix.cols();
      double observedFrequenciesSum = observedFrequencies.sum();
     
      /* get a; 0.8 can be changed to any value in (0,1) */
      double rho = 0.8;
      
      Eigen::VectorXd a = (rho * observedFrequencies) + Eigen::VectorXd::Constant(numberOfRows, observedFrequenciesSum / static_cast<double>(numberOfRows));
      
      if (debug) {
        std::cout << fmt::format("a debug:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(a, false));
      }

      /* get BtSaB -- first term on left side of Eq 49 in H&T */
      Eigen::MatrixXd BtSaB = get_BtSmB(designMatrix, a);

      if (debug) {
        std::cout << fmt::format("BtSaB debug:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(BtSaB));
      }

      /*  get right side of Eq 49 in H&T, which is computed using
      Equation 38 in Holland and Thayer (2000) with all mu = 0 */
      double aLogA = a.array().log().sum();

      Eigen::VectorXd BtSaloga(numberOfColumns);

      for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
        double bALogA = (designMatrix.col(columnIndex).array() * a.array() * a.array().log()).sum();
        double bA = (designMatrix.col(columnIndex) * a).sum();

        BtSaloga(columnIndex) = bALogA - (bA * aLogA / observedFrequenciesSum);
      }

      /* get Beta0 using NR ludcmp() and lubksb() */

      /* right[] is one-offset B in Ax=B */ 
      /* left[][] is one-offset A in Ax=B */
      Eigen::MatrixXd left = Eigen::MatrixXd::Zero(numberOfColumns + 1, numberOfColumns + 1);
      Eigen::VectorXd right = Eigen::VectorXd::Zero(numberOfColumns + 1);
      

      right(Eigen::seq(1, numberOfColumns)) = BtSaloga;

      left.block(1, 1, numberOfColumns, numberOfColumns) = BtSaB;

      Eigen::PartialPivLU<Eigen::MatrixXd> luDecomp = left.partialPivLu();
      Eigen::MatrixXd luMatrix = luDecomp.matrixLU();
      Eigen::VectorXd solution = luDecomp.solve(right);
      
      /* Beta0 is zero-offset solution */
      Eigen::VectorXd beta0(numberOfColumns + 1);
      beta0(Eigen::seq(1, numberOfColumns + 1)) = solution;

      return beta0;
    }
  };
} // namespace EquatingRecipes

#endif