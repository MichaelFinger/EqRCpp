#ifndef FIXTURES_MONDATX_HPP
#define FIXTURES_MONDATX_HPP

#include <fstream>
#include <string>
#include <Eigen/Core>

namespace Tests {
  namespace Fixtures {
    struct MondatX {
      static Eigen::MatrixXd jointRawScores() {
        std::string filename = "../tests/fixtures/mondatx.dat";

        std::ifstream ifs(filename);

        std::string lineRead;
        size_t numberOfRows = 0;

        while (ifs) {
          std::getline(ifs, lineRead);

          numberOfRows++;
        }

        ifs.close();

        ifs.open(filename);

        Eigen::MatrixXd rawScores(numberOfRows, 2);

        size_t rowIndex = 0;
        while (ifs) {
          std::getline(ifs, lineRead);

          if (lineRead.empty()) {
            break;
          }

          rawScores(rowIndex, 0) = std::stod(lineRead.substr(36, 2));
          rawScores(rowIndex, 1) = std::stod(lineRead.substr(38, 2));

          rowIndex++;
        }

        return rawScores;
      }
    };
  } // namespace Fixtures
} // namespace Tests

#endif