#ifndef FIXTURES_ACT_MATH_FREQ_HPP
#define FIXTURES_ACT_MATH_FREQ_HPP

#include <string>
#include <Eigen/Core>

#include "equating_recipe_input_dataset.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Fixtures {
      class ACTMathFreq : public EquatingRecipeInputDataset {
      public:
        struct ACTMathFreqData {
          Eigen::VectorXd rawScores;
          Eigen::VectorXd freqX;
          Eigen::VectorXd freqY;
        };

        ACTMathFreqData data;
        Eigen::VectorXd getParameterEstimates() {
          Eigen::VectorXd parameterEstimates(4);

          parameterEstimates << 0.990464749260068,
              1.800493018403728,
              0.219225988720759,
              1.000000000000000;

          return parameterEstimates;
        }

      protected:
        size_t numberOfRows = 0;
        size_t currentRowIndex;

        void configure() override {
          this->filename = "actmathfreq.dat";
          this->currentRowIndex = 0;
          numberOfRows = this->getNumberOfLines();
        }

        void processLine(const std::string& lineRead) override {
          std::cout << lineRead << "\n";

          std::vector<std::string> values = splitLine(lineRead);

          if (values.size() == 3) {
            double rawScore = std::stod(values[0]);
            double freqX = std::stod(values[1]);
            double freqY = std::stod(values[2]);

            this->data.rawScores(currentRowIndex) = rawScore;
            this->data.freqX(currentRowIndex) = freqX;
            this->data.freqY(currentRowIndex) = freqY;
          }

          currentRowIndex++;
        }

        // static ACTMathFreqData getFreqDist() {
        //   std::string filename = "actmathfreq.dat";

        //   std::ifstream ifs(filename);

        //   std::string lineRead;
        //   size_t numberOfRows = 0;

        //   while (ifs) {
        //     std::getline(ifs, lineRead);

        //     if (lineRead.empty() || lineRead.size() < 20) {
        //       break;
        //     }

        //     numberOfRows++;
        //   }

        //   ifs.close();

        //   ifs.open(filename);

        //   ACTMathFreqData actMathFreqData;
        //   actMathFreqData.rawScores.setZero(numberOfRows);
        //   actMathFreqData.freqX.setZero(numberOfRows);
        //   actMathFreqData.freqY.setZero(numberOfRows);

        //   size_t rowIndex = 0;
        //   while (ifs && rowIndex < numberOfRows) {
        //     std::getline(ifs, lineRead);

        //     if (lineRead.empty()) {
        //       break;
        //     }

        //     actMathFreqData.rawScores(rowIndex) = std::stod(lineRead.substr(0, 2));
        //     actMathFreqData.freqX(rowIndex) = std::stod(lineRead.substr(11, 3));
        //     actMathFreqData.freqY(rowIndex) = std::stod(lineRead.substr(22, 3));

        //     rowIndex++;
        //   }

        //   return actMathFreqData;
        // }

        // Eigen::VectorXd get4ParameterBetaEstimatesX() {
        //   Eigen::VectorXd parameterEstimates;

        //   parameterEstimates(0) = 0.990464749260068;
        //   parameterEstimates(1) = 1.800493018403728;
        //   parameterEstimates(2) = 0.219225988720759;
        //   parameterEstimates(3) = 1.000000000000000;

        //   return parameterEstimates;
        // }
      };
    } // namespace Fixtures
  }   // namespace Tests
} // namespace EquatingRecipes
#endif