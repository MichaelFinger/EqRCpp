#ifndef TESTS_EXAMPLES_DATASETS_YCT_MATH_HPP
#define TESTS_EXAMPLES_DATASETS_YCT_MATH_HPP

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class YctMath {
        public:
          EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;

          YctMath() {
            this->configure();
          }

        private:
          void configure() {            
            EquatingRecipes::JSON::JsonDocument jsonDoc;
            jsonDoc.fromTextFile("./resources/json/datasets/yctmath.txt.json");
            nlohmann::json j = jsonDoc.json;

            Eigen::MatrixXd matrix;
            j.get_to(matrix);

            for (size_t rowIndex = 0; rowIndex < matrix.rows(); rowIndex++) {
              EquatingRecipes::Structures::RawToScaledScoreTable::Entry entry;
              entry.rawScore = matrix(rowIndex, 0);
              entry.scaledScore = matrix(rowIndex, 1);
              entry.scoreLocation = rowIndex;

              rawToScaledScoreTable.lookup[rowIndex] = entry;
            }
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes
#endif