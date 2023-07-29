#ifndef TESTS_EXAMPLES_OBS_SCORE_DIST_HPP
#define TESTS_EXAMPLES_OBS_SCORE_DIST_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>

#include <boost/math/distributions/normal.hpp>

#include <datasets/dummyV3Items.hpp>
#include <datasets/item_specs_file.hpp>
#include <datasets/quadrature_file.hpp>
#include <equating_recipes/analyses/irt_scale_transformation.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <datasets/item_pair_file.hpp>

#include <equating_recipes/implementation/irt_equating.hpp>

#include "datasets/lsat6_items.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      class ObservedScoreDistribution {
      public:
        void operator()() {
          EquatingRecipes::Implementation::IRTEquating irtEquating;

          EquatingRecipes::Tests::Examples::Datasets::LSAT6Items lsat6Items;

          double minTheta = -6.0;
          double maxTheta = 6.0;
          size_t numberOfQuadPoints = 61;

          std::vector<EquatingRecipes::Structures::ItemSpecification> items = lsat6Items();
          size_t numberOfItemsOnTestForm = 5;
          size_t maximumScorePoint = 5;
          Eigen::VectorXd quadraturePoints(numberOfQuadPoints);
          Eigen::VectorXd quadratureWeights(numberOfQuadPoints);
          size_t numberOfScores = 6;
          Eigen::VectorXd scores(numberOfScores);
          Eigen::VectorXd marginalResponseProbabilities(numberOfScores);
          
          double inc = (maxTheta - minTheta) / (static_cast<double>(numberOfQuadPoints - 1));
          boost::math::normal_distribution<double> normalDist(0.0, 1.0);

          for (size_t quadIndex = 0; quadIndex < numberOfQuadPoints; quadIndex++) {
            quadraturePoints(quadIndex) = minTheta + static_cast<double>(quadIndex) * inc;
            quadratureWeights(quadIndex) = pdf(normalDist, quadraturePoints(quadIndex));
          }

          quadratureWeights /= quadratureWeights.sum();

          size_t maximumNumberOfCategories = 2;

          Eigen::MatrixXd conditionalObservedScoreDistribution(numberOfQuadPoints, numberOfScores);
          Eigen::VectorXd xnew;

          for (size_t quadraturePointIndex = 0; quadraturePointIndex < quadraturePoints.size(); quadraturePointIndex++) {
            double theta = quadraturePoints(quadraturePointIndex);

            irtEquating.ObsDistGivenTheta(theta,
                              items,
                              numberOfItemsOnTestForm,
                              maximumNumberOfCategories,
                              maximumScorePoint,
                              numberOfScores,
                              scores,
                              xnew);

            conditionalObservedScoreDistribution.row(quadraturePointIndex) = xnew;
          }

          std::cout.precision(20);
          std::cout << "Conditional Observed Score Distribution:\n";
          for (size_t rowIndex = 0; rowIndex < conditionalObservedScoreDistribution.rows(); rowIndex++) {
            for (size_t columnIndex = 0; columnIndex < conditionalObservedScoreDistribution.cols(); columnIndex++) {
              std::cout << conditionalObservedScoreDistribution(rowIndex, columnIndex) << ", ";
            }

            std::cout << "\n";
          }


          irtEquating.irtMixObsDist(items,
                         numberOfItemsOnTestForm,
                         maximumScorePoint,
                         quadraturePoints,
                         quadratureWeights,
                         numberOfScores,
                         scores,
                         marginalResponseProbabilities);

          std::cout << "Marg Prob: " << marginalResponseProbabilities << "\n";
        }
      };
    }
  }
}

#endif