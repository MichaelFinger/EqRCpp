#ifndef TESTS_EXAMPLES_DATASETS_QUADRATURE_FILE_HPP
#define TESTS_EXAMPLES_DATASETS_QUADRATURE_FILE_HPP

#include <iostream>
#include <istream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

#include <boost/algorithm/string.hpp>

#include <equating_recipes/structures/quadrature.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class QuadratureFile {
        public:
          EquatingRecipes::Structures::Quadrature quadrature;

          void import(const std::string& filename) {
            importQuadratureFile(filename);
          }

        private:
          void importQuadratureFile(const std::string& filename) {
            std::vector<std::string> linesRead;

            std::ifstream ifs(filename,
                              std::ios::in);

            for (std::string lineRead; std::getline(ifs, lineRead);) {
              linesRead.push_back(lineRead);
            }

            ifs.close();

            size_t lineIndex = 0;

            size_t numberOfQuadraturePoints = static_cast<size_t>(std::stoi(linesRead[0]));
            quadrature.thetaValues.resize(numberOfQuadraturePoints);
            quadrature.thetaWeights.resize(numberOfQuadraturePoints);

            for (size_t quadIndex = 0; quadIndex < numberOfQuadraturePoints; quadIndex++) {
              lineIndex++;

              std::vector<std::string> values = splitString(linesRead[lineIndex]);

              quadrature.thetaValues(quadIndex) = std::stod(values[0]);
              quadrature.thetaWeights(quadIndex) = std::stod(values[1]);
            }
          }

          std::vector<std::string> splitString(const std::string& lineRead) {
            std::vector<std::string> values;

            boost::split(values, lineRead, boost::is_any_of(" "), boost::algorithm::token_compress_on);

            return values;
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes

#endif