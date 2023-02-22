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
    short trueScoreEquating(EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                            std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                            std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                            const size_t& numberOfScores,
                            const Eigen::VectorXd& newScores,
                            Eigen::VectorXd& oldFormEquivalents,
                            Eigen::VectorXd& thetaValues,
                            double& lowestObservableScoreNewForm,
                            double& lowestObservableScoreOldForm) {
      double trueMinNew = 0.0;
      double trueMinOld = 0.0;
      double oldScoreMax = 0.0; /* Maximum old form score */
      double newTestMin = 0;    /* Minimum new form test score */
      double oldTestMin = 0;    /* Minimum old form test score */

      // int j, r, CatMax;
      // long chance;
      // double slope;
      // double intercept;
      // double xh;
      // double xl;
      // double x1;
      // double x2;

      // double newScoreMin;       /* minimum new form score */
      // double newScoreInc;       /* Increment between consecutive new form scores */

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
      for (size_t chance = 0; chance < numberOfScores - 1; chance++) {
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
        eqvOld(scoreLocation) = intercept + newScores(scoreLocation) * slope;
        theta(scoreLocation) = xl;
      }

      double x1 = -99.0;
      double x2 = 99.0;

      for (size_t scoreLocation = chance + 1; scoreLocation < numberOfScores - 1; scoreLocation++) {
        this->trueS = newScores(scoreLocation);

        for (size_t r = 0; r <= 10; r++) {
          theta[j] = er_rtsafe(funcd, x1, x2, 0.00001);
          /* ===== End of version 1.0 update ===== */
          if (theta[j] == -9999.0)
            theta[j] = -99.0; /* tianyou added this line 7/18/08 */
          if (theta[j] != -9999.0 && theta[j] != -99.0 && theta[j] != 99.0)
            break;
          else {
            x1 = -pow(2.0, (double)r);
            x2 = pow(2.0, (double)r);
          }
        }
        eqvOld[j] = trueScore(OldItems, Handle->OldItemNum, theta[j]);
      }

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
    double testResponseFunction(const double& theta,
                                const std::vector<EquatingRecipes::Structures::ItemSpecification>& items) {
      double testResponseFunctionValue = 0.0;

      std::for_each(items.begin(),
                    items.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      for (size_t categoryIndex = 0; categoryIndex < item.numberOfCategories; categoryIndex++) {
                        double probResponse = probabilityItemResponse(item, categoryIndex, theta);
                        testResponseFunctionValue += item.scoringFunctionValues(categoryIndex) * probResponse;
                      }
                    });

      return (this->trueS - testResponseFunctionValue);
    }

    double probabilityItemResponse(const EquatingRecipes::Structures::ItemSpecification& item,
                                   const size_t& categoryIndex,
                                   const double& theta) {
      return -1;
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double testResponseFunctionDerivative(const double& theta) {
      // int j, k;
      // double pd, Wjk, v;

      // v = 0.0;
      // for (j = 1; j <= ContHandle->NewItemNum; j++) {
      //   for (k = 1; k <= ContNewItems[j].CatNum; k++) {
      //     pd = PdCCCoverTheta(&ContNewItems[j], k, theta);
      //     Wjk = ContNewItems[j].ScoreFunc[k];
      //     v += Wjk * pd;
      //   }
      // }

      // return -v;

      return 0;
    }

//     void funcd(double x, double* f, double* fd)
//     /*
//   Author: Seonghoon Kim
//   Date of last revision 9/25/08
// */
//     {
//       *f = f_mix(x);
//       *fd = f_mixDer(x);

//       return;
//     }
  };
} // namespace EquatingRecipes

#endif