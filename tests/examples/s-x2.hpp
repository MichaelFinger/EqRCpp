#ifndef TESTS_EXAMPLES_S_X2_HPP
#define TESTS_EXAMPLES_S_X2_HPP

#include <iostream>
#include <Eigen/Core>
#include <equating_recipes/implementation/observed_score_distribution.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      class S_X2 {
      public:
        void operator()() {
          double theta = 1.0;
          size_t numberOfItems = 5;

          std::optional<Eigen::VectorXd> scoringFunction;

          EquatingRecipes::Structures::ItemSpecification item1 = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(1, 0.8250552, -3.3607460, 1.0, scoringFunction, true);
          EquatingRecipes::Structures::ItemSpecification item2 = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(2, 0.7230608, -1.3695513, 1.0, scoringFunction, true);
          EquatingRecipes::Structures::ItemSpecification item3 = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(3, 0.8899989, -0.2798928, 1.0, scoringFunction, true);
          EquatingRecipes::Structures::ItemSpecification item4 = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(4, 0.6886588, -1.8657302, 1.0, scoringFunction, true);
          EquatingRecipes::Structures::ItemSpecification item5 = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(5, 0.6575904, -3.1229745, 1.0, scoringFunction, true);

          std::vector<EquatingRecipes::Structures::ItemSpecification> items {item1,
                                                                             item2,
                                                                             item3,
                                                                             item4,
                                                                             item5};

          size_t numberOfScores;
          Eigen::VectorXd scores;
          Eigen::VectorXd xnew;

          EquatingRecipes::Implementation::ObservedScoreDistribution observedScoreDistribution;

          observedScoreDistribution.ObsDistGivenTheta(theta,
                                                      items,
                                                      numberOfScores,
                                                      scores,
                                                      xnew);

          std::cout << numberOfScores << "\n";
          std::cout << scores << "\n";
          std::cout << xnew << "\n";
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif