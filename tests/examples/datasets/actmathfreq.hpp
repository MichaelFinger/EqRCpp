#ifndef TESTS_EXAMPLES_DATASETS_ACT_MATH_FREQ_HPP
#define TESTS_EXAMPLES_DATASETS_ACT_MATH_FREQ_HPP

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>


namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class ACTMathFreq {
        public:
          Eigen::VectorXd rawScores;
          Eigen::VectorXd freqX;
          Eigen::VectorXd freqY;

          ACTMathFreq() {
            this->configure();
          }

          Eigen::VectorXd getParameterEstimates() {
            Eigen::VectorXd parameterEstimates(4);

            parameterEstimates << 0.990464749260068,
                1.800493018403728,
                0.219225988720759,
                1.000000000000000;

            return parameterEstimates;
          }

        private:
          void configure() {            
            EquatingRecipes::JSON::JsonDocument jsonDoc;
            jsonDoc.fromTextFile("./resources/json/datasets/actmathfreq.dat.json");
            nlohmann::json j = jsonDoc.json;

            Eigen::MatrixXd matrix;
            j.get_to(matrix);

            this->rawScores = matrix.col(0);
            this->freqX = matrix.col(1);
            this->freqY = matrix.col(2);
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes
#endif