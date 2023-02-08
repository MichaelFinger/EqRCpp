#ifndef FIXTURES_ACT_MATH_FREQ_HPP
#define FIXTURES_ACT_MATH_FREQ_HPP

#include <fstream>
#include <string>
#include <Eigen/Core>

namespace Tests {
  namespace Fixtures {
    struct ACTMathFreqData {
      Eigen::VectorXd rawScores;
      Eigen::VectorXd freqX;
      Eigen::VectorXd freqY;
    };

    struct ACTMathFreq {
      static ACTMathFreqData getFreqDist() {
        std::string filename = "actmathfreq.dat";

        std::ifstream ifs(filename);

        std::string lineRead;
        size_t numberOfRows = 0;

        while (ifs) {
          std::getline(ifs, lineRead);

          if (lineRead.empty() || lineRead.size() < 20) {
            break;
          }

          numberOfRows++;
        }

        ifs.close();

        ifs.open(filename);

        ACTMathFreqData actMathFreqData;
        actMathFreqData.rawScores.setZero(numberOfRows);
        actMathFreqData.freqX.setZero(numberOfRows);
        actMathFreqData.freqY.setZero(numberOfRows);

        size_t rowIndex = 0;
        while (ifs && rowIndex < numberOfRows) {
          std::getline(ifs, lineRead);

          if (lineRead.empty()) {
            break;
          }

          actMathFreqData.rawScores(rowIndex) = std::stod(lineRead.substr(0, 2));
          actMathFreqData.freqX(rowIndex) = std::stod(lineRead.substr(11, 3));
          actMathFreqData.freqY(rowIndex) = std::stod(lineRead.substr(22, 3));

          rowIndex++;
        }

        return actMathFreqData;
      }

      static Eigen::VectorXd get4ParameterBetaEstimatesX() {
        Eigen::VectorXd parameterEstimates(4);

        parameterEstimates(0) = 0.990464749260068;
        parameterEstimates(1) = 1.800493018403728;
        parameterEstimates(2) = 0.219225988720759;
        parameterEstimates(3) = 1.000000000000000;

        return parameterEstimates;
      }
    };
  } // namespace Fixtures
} // namespace Tests

#endif