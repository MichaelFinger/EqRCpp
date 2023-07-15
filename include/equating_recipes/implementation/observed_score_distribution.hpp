#ifndef OBSERVED_SCORE_DISTRIBUTION_HPP
#define OBSERVED_SCORE_DISTRIBUTION_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/distributions/normal.hpp>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

namespace EquatingRecipes {
  namespace Implementation {
    class ObservedScoreDistribution {
    public:
      /*------------------------------------------------------------------------------
      Functionality:
        Calculates the conditional distribution of total score given theta
        for IRT models. The IRT models include the 3PL, LGR, GPC, and NR models.
        Uses Hanson's (1994) generalization of the Lord-Wingersky (1982) recursive
        algorithm.

      Input:
        theta:   examinee ability
        Items:   pointer to designate items on a test form (0-offset)
        n :      number of items on a test form
        MaxCat:  the maximum number of categories across items
        MaxScrP: maximum number of score points

      Output:
        nscr:    number of score categories for observed score distribution
        xscr:    vector of scores associated with each score category
                (are consecutive integers); 0-offset
        xnew:    vector of probabilities associated with each score category;
                0-offset

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      enum class IRTModel {
        THREE_PARAMETER,
        GRADED_RESPONSE,
        GENERALIZED_PARTIAL_CREDIT,
        GRADED_RESPONSE_RATING_SCALE,
        GENERALIZED_PARTIAL_CREDIT_RATING_SCALE
      };

      enum class Density {
        LOGISTIC,
        NORMAL
      };

      struct ItemParameters {
        size_t itemIndex;
        double slope;
        double location;
        double lowerAsymptote = 0;
        Eigen::VectorXd thresholds;
        Eigen::VectorXd categoryParameters;
        Density density;
        IRTModel irtModel;
        double scalingConstant = 1.0;
        bool locationIsIntercept = false;
      };

      Eigen::VectorXd itemResponseProbabilities(const ItemParameters& itemParameters,
                                                const double& theta) {
        Eigen::VectorXd respProbs;

        if (itemParameters.irtModel == IRTModel::THREE_PARAMETER) {
          respProbs.setZero(2);
          double deviate;
          if (itemParameters.locationIsIntercept) {
            deviate = itemParameters.scalingConstant * (itemParameters.slope * theta + itemParameters.location);
          } else {
            deviate = itemParameters.scalingConstant * (itemParameters.slope * (theta - itemParameters.location));
          }

          respProbs(1) = itemParameters.lowerAsymptote + (1.0 - itemParameters.lowerAsymptote) * getCDF(deviate, itemParameters.density);
          respProbs(0) = 1.0 - respProbs(1);

        } else if (itemParameters.irtModel == IRTModel::GRADED_RESPONSE) {
          size_t numberOfScoreCategories = itemParameters.thresholds.size() + 1;
          respProbs.setZero(numberOfScoreCategories);

          Eigen::VectorXd cumProbs(numberOfScoreCategories + 1);
          cumProbs(0) = 1.0;
          cumProbs(numberOfScoreCategories) = 0.0;

          for (size_t itemResponseIndex = 1; itemResponseIndex < numberOfScoreCategories; itemResponseIndex++) {
            double deviate;

            if (itemParameters.locationIsIntercept) {
              deviate = itemParameters.scalingConstant * (itemParameters.slope * theta + itemParameters.thresholds(itemResponseIndex - 1));
            } else {
              deviate = itemParameters.scalingConstant * (itemParameters.slope * (theta - itemParameters.thresholds(itemResponseIndex - 1)));
            }

            cumProbs(itemResponseIndex) = getCDF(deviate, itemParameters.density);

            respProbs(itemResponseIndex) = cumProbs(itemResponseIndex - 1) - cumProbs(itemResponseIndex);
          }

        } else if (itemParameters.irtModel == IRTModel::GENERALIZED_PARTIAL_CREDIT) {
        } else if (itemParameters.irtModel == IRTModel::GRADED_RESPONSE_RATING_SCALE) {
        } else if (itemParameters.irtModel == IRTModel::GENERALIZED_PARTIAL_CREDIT_RATING_SCALE) {
        }

        size_t numberOfScoreCategories;

        if (itemParameters.thresholds.size() >= 1) {
          numberOfScoreCategories = itemParameters.thresholds.size() + 1;
        } else if (itemParameters.categoryParameters.size() >= 1) {
          numberOfScoreCategories = itemParameters.categoryParameters.size();
        } else {
          numberOfScoreCategories = 2;
        }

        return respProbs;
      }

      void ObsDistGivenTheta(const double& theta,
                             const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                             const size_t& numberOfItemsOnForm,
                             const size_t& maximumCategoryIndex,
                             const size_t& maximumScorePoint,
                             size_t& numberOfScores,
                             Eigen::VectorXd& scores,
                             Eigen::VectorXd& xnew) {
        // int i, j, k, index;
        // int mino, maxo, minn, maxn;
        Eigen::VectorXd xitem(maximumCategoryIndex + 1); /* zero-offset, but not use xitem[0] */
        Eigen::VectorXd xold(maximumScorePoint);         /* zero-offset */

        EquatingRecipes::Implementation::IRTModelFunctions irtModelFunctions;

        /* calculates probabilities for Item 1 */
        for (size_t categoryIndex = 0; categoryIndex < items[0].numberOfCategories; categoryIndex++) {
          xitem(categoryIndex) = irtModelFunctions.itemResponseFunction(items[0], categoryIndex, theta);
        }

        double mino = items[0].scoringFunctionValues(0);
        double maxo = items[0].scoringFunctionValues(items[0].numberOfCategories - 1);
        double minn = mino;
        double maxn = maxo;

        xold.setZero();

        for (size_t categoryIndex = 0; categoryIndex < items[0].numberOfCategories; categoryIndex++) {
          size_t index = static_cast<size_t>(items[0].scoringFunctionValues(categoryIndex) - minn);
          xold(index) = xitem(categoryIndex); /* mino associated with index of 0 */
        }                                     /* mino does vary; see below      */

        xnew = xold;

        if (numberOfItemsOnForm == 1) {
          size_t maxMinusMin = static_cast<size_t>(maxn - minn);

          for (size_t score = 0; score <= maxMinusMin; score++) {
            scores(score) = static_cast<double>(score) + minn;
          }

          numberOfScores = static_cast<size_t>(maxn - minn + 1);

          return;
        }

        /* updates distribution for items 2 through nitems */
        for (size_t itemIndex = 1; itemIndex < numberOfItemsOnForm; itemIndex++) {
          for (size_t categoryIndex = 0; categoryIndex < items[itemIndex].numberOfCategories; categoryIndex++) {
            xitem(categoryIndex) = irtModelFunctions.itemResponseFunction(items[itemIndex],
                                                                          categoryIndex,
                                                                          theta);
          }

          recurs(mino,
                 maxo,
                 xold,
                 items[itemIndex].numberOfCategories,
                 items[itemIndex].scoringFunctionValues,
                 xitem,
                 minn,
                 maxn,
                 xnew);

          mino = minn;
          maxo = maxn;

          size_t maxIndex = static_cast<size_t>(maxn - minn);

          xold(Eigen::seq(0, maxIndex)) = xnew(Eigen::seq(0, maxIndex));
        }

        size_t maxIndex = static_cast<size_t>(maxn - minn);
        for (size_t index = 0; index <= maxIndex; index++) {
          scores(index) = static_cast<double>(index) + minn;
        }

        numberOfScores = static_cast<size_t>(maxn - minn + 1);
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Updates a distribution of scores using Hanson's (1994) generalization of
        the Lord-Wingersky (1982) formula.

        Assumes that test scores are consecutive integers
        Assumes that item scores are consecutive integers
        Assumes that test and item scores are sorted from low to high.

      Input (r-1 added items):
        mino : minimum integer score for old distribution
        maxo : maximum integer score for old distribution
        xold : probability array for each score point from mino to maxo
              (f_r-1 (x|theta)); 0-offset
        mitem: number of distinct score points for the added (new) item
        iitem: values of the scoring function for the added item; 0-offset
        xitem: values of category response functions for the added item; 0-offset

      Output (r added items):
        minn : minimum integer score for new (updated) distribution
        maxn : maximum integer score for new (updated) distribution
        xnew : probability array for each score point from minn to maxn
              (f_r(x|theta)); 0-offset

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      void recurs(const int& mino,
                  const int& maxo,
                  const Eigen::VectorXd& xold,
                  const size_t& numberOfCategories,
                  const Eigen::VectorXd& iitem,
                  const Eigen::VectorXd& xitem,
                  double& minn,
                  double& maxn,
                  Eigen::VectorXd& xnew) {
        minn = mino + iitem(0);
        maxn = maxo + iitem(numberOfCategories - 1);

        for (double i = minn; i <= maxn; i += 1.0) {
          double in = i - minn;

          xnew(static_cast<size_t>(in)) = 0.0;

          for (size_t j = 0; j < numberOfCategories; j++) {
            double io = i - iitem(j) - mino;

            if (io >= 0 && io <= maxo - mino) {
              xnew(static_cast<size_t>(in)) += xold(static_cast<size_t>(io)) * xitem(j);
            }
          }
        }
      }

    // private:
      boost::math::normal_distribution<> normalDist;

      double getCDF(const double& deviate, const Density& density) {
        double cdfValue;

        if (density == Density::LOGISTIC) {
          cdfValue = getLogisticCDF(deviate);
        } else {
          cdfValue = getNormalCDF(deviate);
        }

        return cdfValue;
      }

      double getLogisticCDF(const double& deviate) {
        double cdfValue = 1.0 / (1.0 + std::exp(-1.0 * deviate));
        return cdfValue;
      }

      double getNormalCDF(const double& deviate) {
        double cdfValue = boost::math::cdf(normalDist, deviate);
        return cdfValue;
      }

      Eigen::VectorXd probResp3P(const double& slope,
                                 const double& location,
                                 const double& lowerAsymptote,
                                 const double& theta,
                                 const double& scalingConstant,
                                 const Density& density,
                                 const bool& locationIsIntercept) {
        Eigen::VectorXd respProbs(2);

        double deviate;

        if (locationIsIntercept) {
            deviate = (slope * theta + location);
          } else {
            deviate = slope * (theta - location);
          }

        if (density == Density::LOGISTIC) {
          deviate *= scalingConstant;
        }
        
        respProbs(1) = lowerAsymptote + (1.0 - lowerAsymptote) * getCDF(deviate, density);
        respProbs(0) = 1.0 - respProbs(1);

        return respProbs;
      }

      Eigen::VectorXd probRespGR(const double& slope,
                                 const Eigen::VectorXd& thresholds,
                                 const double& theta,
                                 const double& scalingConstant,
                                 const Density& density,
                                 const bool& thresholdIsIntercept) {
        size_t numberOfCategories = thresholds.size() + 1;

        Eigen::VectorXd cumProbs(numberOfCategories + 1);
        Eigen::VectorXd respProbs(numberOfCategories);

        cumProbs(0) = 1.0;
        cumProbs(numberOfCategories) = 0.0;
        
        for (size_t itemResp = 1; itemResp < numberOfCategories; itemResp++) {
          std::cout << itemResp << "\n";

          double deviate;

          if (thresholdIsIntercept) {
            deviate = slope * theta + thresholds(itemResp - 1);
          } else {
            deviate = slope * (theta - thresholds(itemResp - 1));
          }

          if (density == Density::LOGISTIC) {
            deviate *= scalingConstant;
          }
        
          cumProbs(itemResp) = getCDF(deviate, density);

          respProbs(itemResp - 1) = cumProbs(itemResp - 1) - cumProbs(itemResp);
        }

        respProbs(numberOfCategories - 1) = cumProbs(numberOfCategories - 1) - cumProbs(numberOfCategories);

        return respProbs;
      }
    };
  } // namespace Implementation
} // namespace EquatingRecipes

#endif