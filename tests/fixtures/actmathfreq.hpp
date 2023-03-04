#ifndef FIXTURES_ACT_MATH_FREQ_HPP
#define FIXTURES_ACT_MATH_FREQ_HPP

#include <Eigen/Core>
#include <nlohmann/json.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Fixtures {
      class ACTMathFreq {
      public:
        Eigen::VectorXd rawScores;
        Eigen::VectorXd freqX;
        Eigen::VectorXd freqY;

        Eigen::VectorXd getParameterEstimates() {
          Eigen::VectorXd parameterEstimates(4);

          parameterEstimates << 0.990464749260068,
              1.800493018403728,
              0.219225988720759,
              1.000000000000000;

          return parameterEstimates;
        }

        void configure(const nlohmann::json& j) {
          size_t numberOfRows = j.size();
          size_t numberOfColumns = j[0].size();

          Eigen::MatrixXd matrix(numberOfRows, numberOfColumns);

          for (size_t rowIndex = 0; rowIndex < numberOfRows; rowIndex++) {
            nlohmann::json jsonRow = j[rowIndex];

            for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
              matrix(rowIndex, columnIndex) = jsonRow[columnIndex];
            }
          }

          this->rawScores = matrix.col(0);
          this->freqX = matrix.col(1);
          this->freqY = matrix.col(2);
        }
      };
    } // namespace Fixtures
  }   // namespace Tests
} // namespace EquatingRecipes
#endif