#ifndef TESTS_SOLVING_HPP
#define TESTS_SOLVING_HPP

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <fmt/core.h>
#include <equating_recipes/utilities.hpp>

namespace Tests {
  class Solving {
  public:
    void run() {
      // Eigen::MatrixXd A(3, 3);
      // Eigen::VectorXd b(3);
      // A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
      // b << 3, 3, 4;
      // std::cout << "Here is the matrix A:\n"
      //           << A << std::endl;
      // std::cout << "Here is the vector b:\n"
      //           << b << std::endl;
      // Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
      // std::cout << "The solution is:\n"
      //           << x << std::endl;

      // Eigen::MatrixXd lu = A.partialPivLu().matrixLU();

      // A x = B, solve for x
      // Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

      // std::cout << x << "\n";

      // Eigen::MatrixXd X = Eigen::MatrixXd::Random(10, 10);
      // Eigen::MatrixXd A = X + X.transpose();
      // std::cout << "Here is a random symmetric 4x4 matrix:\n"
      //           << A << "\n";

      // Eigen::Tridiagonalization<Eigen::MatrixXd> triOfA(A);
      // Eigen::MatrixXd pm = triOfA.packedMatrix();
      // std::cout << "The packed matrix M is:\n"
      //           << pm << "\n";
      // std::cout << "The diagonal and subdiagonal corresponds to the matrix T, which is:\n"
      //           << triOfA.matrixT() << "\n";

      Eigen::MatrixXd x = Eigen::MatrixXd::Random(11, 1);

      Eigen::MatrixXd h = x(Eigen::seq(1, 10), 0) - x(Eigen::seq(0, 9), 0);
      // std::cout << h << "\n";

      // Eigen::Tridiagonalization<Eigen::MatrixXd> triOfH(h * h.transpose());
      // std::cout << "The diagonal and subdiagonal corresponds to the matrix T, which is:\n"
      //           << triOfH.matrixT() << "\n";

      setupQt(true, x, 11, h);

      std::cout << x << "\n";
    }

  private:
    void setupQt(const bool& edist,
                 Eigen::MatrixXd& mqt,
                 const size_t n,
                 Eigen::MatrixXd& h) {
      mqt.setZero(n + 1, n - 1);

      double h0 = h(0);

      size_t maxIndex = (mqt.rows() <= mqt.cols()) ? mqt.rows() : mqt.cols();

      for (size_t index = 0; index < maxIndex; index++) {
        if (edist) {
          h0 = h(0, 0);

          if (index >= 1) {
            mqt(index - 1, index) = 1.0 / h0;
          }

          mqt(index, index) = -2.0 / h0;

          mqt(index + 1, index) = 1.0 / h0;
        } else {
          if (index >= 1) {
            mqt(index - 1, index) = 1.0 / h(index, 0);
          }

          mqt(index, index) = -1.0 * (1.0 / h(index - 1, 0) + 1.0 / h(index, 0));

          mqt(index + 1, index) = 1.0 / h(index, 0);
        }
      }
    }
  };
} // namespace Tests

#endif