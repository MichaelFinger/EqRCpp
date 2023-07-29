#ifndef TESTS_EXAMPLES_DATASETS_LSAT6_ITEMS_HPP
#define TESTS_EXAMPLES_DATASETS_LSAT6_ITEMS_HPP

#include <vector>
#include <Eigen/Core>
#include <equating_recipes/structures/item_specification.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class LSAT6Items {
        public:
          std::vector<EquatingRecipes::Structures::ItemSpecification> operator()() {
            std::vector<double> a {0.485128357458009,
                                   0.4245383922115,
                                   0.523767089092307,
                                   0.404365650412287,
                                   0.385782845075016};

            std::vector<double> b {-3.35871448703245,
                                   -1.37033461537717,
                                   -0.279531663951713,
                                   -1.8666908185003,
                                   -3.12694187447321};

            std::vector<EquatingRecipes::Structures::ItemSpecification> items;

            for (size_t itemIndex = 0; itemIndex < a.size(); itemIndex++) {
              std::optional<Eigen::VectorXd> scoringFunctionValues;
              EquatingRecipes::Structures::ItemSpecification item = EquatingRecipes::Structures::ItemSpecification::buildTwoParameterLogistic(itemIndex + 1,
                                                                                                                                              a[itemIndex],
                                                                                                                                              b[itemIndex],
                                                                                                                                              1.702,
                                                                                                                                              scoringFunctionValues,
                                                                                                                                              true);

              items.push_back(item);
            }

            return items;
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes

#endif