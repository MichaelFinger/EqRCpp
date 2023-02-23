/*
  IRT True Score Equating for Mixed-Format Tests

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

  For details about partial derivatives, refer to the following report:

  Kim, S., & Kolen, M.J. (2005). Methods for obtaining a common scale under
      unidimensional IRT models: A technical review and further extensions
      (Iowa Testing Programs Occasional Paper, No. 52). The University of
      Iowa.

  Note: All pointers are of 0-offset unless otherwise indicated.
*/

#ifndef IRT_EQUATING_HPP
#define IRT_EQUATING_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_control.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>

namespace EquatingRecipes {
  class IRTEquating {
  public:
    /*---------------------------------------------------------------------------
      Functionality:
        Performs IRT true score equating with mixed-format tests.
      
      Input:
        Handle      A pointer to control the computing environment for equating
                    The same type of structure as used for IRT scale transformation
                    is used.
        NewItems    A pointer to designate items on the new form (0-offset)
        OldItems    A pointer to designate items on the old form (0-offset)
        nScores     Nnumber of new form scores for which old form equivalents
                    are to be computed
        newScores   New form scores for which old form equivalents are to be
                    computed. These scores are assumed to be equally spaced;
                    the difference between consecutive elements is assumed to be
                    constant (0-offset).

      Output:
        eqvOld      vector of old form equivalents of new form scores; 0-offset
        theta       vector of theta values which produce the vector of integer new
                    form true scores.
                    
        *NewMin     the lowest possible score on the new form.
        *OldMin     the lowest possible score on the old form.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    short trueScoreEquating(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                            const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                            const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                            const size_t& numberOfScores,
                            const Eigen::VectorXd& newScores,
                            Eigen::VectorXd& oldFormEquivalents,
                            Eigen::VectorXd& theta,
                            double& lowestObservableScoreNewForm,
                            double& lowestObservableScoreOldForm) {
      double trueMinNew = 0.0;
      double trueMinOld = 0.0;
      double oldScoreMax = 0.0; /* Maximum old form score */
      double newTestMin = 0;    /* Minimum new form test score */
      double oldTestMin = 0;    /* Minimum old form test score */

      this->controlHandle = handle;
      this->controlNewItems = newItems;

      std::for_each(newItems.begin(),
                    newItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& newItem) {
                      if (newItem.irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                        trueMinNew += newItem.scoringFunctionValues(0) * (1.0 - newItem.c(1) +
                                                                          newItem.scoringFunctionValues(1) * newItem.c(1));
                      } else {
                        trueMinNew += newItem.scoringFunctionValues(0);
                      }

                      newTestMin += newItem.scoringFunctionValues(0);
                    });

      std::for_each(oldItems.begin(),
                    oldItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& oldItem) {
                      if (oldItem.irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                        trueMinOld += (oldItem.scoringFunctionValues(0) + (oldItem.scoringFunctionValues(1)) * oldItem.c(1));
                      } else {
                        trueMinOld += oldItem.scoringFunctionValues(0);
                      }

                      size_t maximumCategory = oldItem.numberOfCategories;
                      oldScoreMax += oldItem.scoringFunctionValues(maximumCategory - 1);
                      oldTestMin += oldItem.scoringFunctionValues(0);
                    });

      /* Find index of largest score below sum of c's */
      size_t chance;
      for (chance = 0; chance < numberOfScores - 1; chance++) {
        if (newScores(chance) > trueMinNew) {
          break;
        }
      }

      /* Compute slope and intercept of interpolating function for 
        new form scores below chance */
      double slope = trueMinNew - newTestMin;

      if (slope > 0) {
        slope = (trueMinOld - oldTestMin) / slope;
      }

      double intercept = oldTestMin - slope * newTestMin;

      double xh = 99;
      double xl = -99;

      this->trueS = newScores(chance); /* Check if true score at theta=-99 is greater      */
      if (f_mix(xl) < 0.0) {           /* than computed largest score less than the        */
        ++chance;                      /* the sum of the c's, if so increase chance index. */
      }

      /* Assign score equivalents */
      for (size_t scoreLocation = 0; scoreLocation <= chance; scoreLocation++) {
        oldFormEquivalents(scoreLocation) = intercept + newScores(scoreLocation) * slope;
        theta(scoreLocation) = xl;
      }

      double x1 = -99.0;
      double x2 = 99.0;

      for (size_t scoreLocation = chance + 1; scoreLocation < numberOfScores - 1; scoreLocation++) {
        this->trueS = newScores(scoreLocation);

        for (size_t r = 0; r <= 10; r++) {
          theta(scoreLocation) = er_rtsafe(x1, x2, 0.00001);

          if (theta(scoreLocation) == -9999.0) {
            theta(scoreLocation) = -99.0; /* tianyou added this line 7/18/08 */
          }

          if (theta(scoreLocation) != -9999.0 &&
              theta(scoreLocation) != -99.0 &&
              theta(scoreLocation) != 99.0) {
            break;
          } else {
            x1 = -1.0 * std::pow(2.0, static_cast<double>(r));
            x2 = std::pow(2.0, static_cast<double>(r));
          }
        }

        oldFormEquivalents(scoreLocation) = trueScore(oldItems,
                                                      theta(scoreLocation));
      }

      /* Convert maximum score on new form to maximum score on old form */
      oldFormEquivalents(numberOfScores - 1) = oldScoreMax;
      theta(numberOfScores - 1) = xh;

      lowestObservableScoreNewForm = trueMinNew;
      lowestObservableScoreOldForm = trueMinOld;

      return 0;
    }

  private:
    EquatingRecipes::Structures::IRTScaleTransformationControl controlHandle;
    std::vector<EquatingRecipes::Structures::ItemSpecification> controlNewItems;
    double trueS;

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double trueScore(const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                     const double& theta) {
      double expectedRawScore = 0.0;

      std::for_each(items.begin(),
                    items.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      for (size_t categoryIndex = 0; categoryIndex < item.scoringFunctionValues.size(); categoryIndex++) {
                        double probResponse = itemResponseFunction(item,
                                                                   categoryIndex,
                                                                   theta);

                        expectedRawScore += probResponse * item.scoringFunctionValues(categoryIndex);
                      }
                    });

      return expectedRawScore;
    }

    short IRTmixObsEq(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle, 
    const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
		const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems, 
    const double& wNew, 
    const double& wOld,
    EquatingRecipes::Structures::IRTFittedDistribution& newForm,
    EquatingRecipes::Structures::IRTFittedDistribution& oldForm,
		EquatingRecipes::Structures::IRTEquatingResults& irtEquatingResults) {
      return 0;
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double f_mix(const double& theta) {
      double expectedRawScore = trueScore(this->controlNewItems,
                                          theta);

      return (this->trueS - expectedRawScore);
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double f_mixDer(const double& theta) {
      // int j, k;
      // double pd, Wjk, v;

      double v = 0.0;
      std::for_each(this->controlNewItems.begin(),
                    this->controlNewItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      for (size_t categoryIndex = 0; categoryIndex < item.scoringFunctionValues.size(); categoryIndex++) {
                        double pd = PdCCCoverTheta(item, categoryIndex, theta);

                        v += pd * item.scoringFunctionValues(categoryIndex);
                      }
                    });

      v *= -1.0;

      return v;
    }

    double PdCCCoverTheta(const EquatingRecipes::Structures::ItemSpecification& item,
                          const size_t& categoryIndex,
                          const double& theta) {
      return 0;
    }

    double itemResponseFunction(const EquatingRecipes::Structures::ItemSpecification& item,
                                const size_t& categoryIndex,
                                const double& theta) {
      return -1;
    }

    /*
      Purpose:  
          This function implements bisection method to find x such that

              f(x)==0.0.

      Input:             
          funcd   function pointer whose root we are trying to find.
        x0      left end of the function domain 
        x1      right end of the function domain

      Output:                
          x that satisfies f(x)==0.0
      
      Precondition:                                                     
        f(x) is defined on [x0,x1] with f(x0)*f(x1) < = 0.
    
      Author: Jaehoon Seol
      Date: August 20, 2009
      Version: 1.0
      References : 
          J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
        David S. Watkins, Fundamentals of Matrix Computations,2002
      Comments:
    */
    double er_rtsafe(double x0, double x1, double error) {
      double left = testResponseFunction(x0);
      double right = testResponseFunction(x1);

      /* assert that there is a root between x1 and x2 */
      if (left * right > 0) { /* function values have the same sign */
        std::string msg = "Equating Recipes error occured\n";
        msg.append("Source: er_find_root, Error: no root found exists on the interval\n");
        throw std::runtime_error(msg);
      }

      double diff = std::abs(right - left);

      while (diff > error) {
        double mid = (left + right) / 2.0;
        double side1 = testResponseFunction(mid);
        double side2 = testResponseFunction(right);

        if (side1 * side2 <= 0) {
          left = mid;
        } else {
          right = mid;
        }

        diff = std::abs(right - left);
      }

      return mid;
    }
  };
} // namespace EquatingRecipes

#endif