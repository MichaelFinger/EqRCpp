#ifndef TESTS_MATRIX_DECOMP_HPP
#define TESTS_MATRIX_DECOMP_HPP

#include <iostream>
#include <Eigen/Dense>

namespace Tests {
  struct MatrixDecomp {
    void run() {
      Eigen::MatrixXd A(3, 3);
      Eigen::VectorXd b(3);
      A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
      b << 3, 3, 4;
      std::cout << "Here is the matrix A:\n"
                << A << std::endl;
      // std::cout << "Here is the vector b:\n"
      //           << b << std::endl;
      // Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
      // std::cout << "The solution is:\n"
      //           << x << std::endl;

      // Eigen::MatrixXd lu = A.partialPivLu().matrixLU();
      
      // A x = B, solve for x
      Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

      std::cout << x << "\n";
    }
  };
} // namespace Tests

#endif