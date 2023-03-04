#ifndef FIXTURES_MONDATX_HPP
#define FIXTURES_MONDATX_HPP

#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Fixtures {
      struct MondatX {
        Eigen::MatrixXd itemResponseMatrix;
        Eigen::MatrixXd rawScores;

        void configure(const nlohmann::json& j) {
          Eigen::MatrixXd matrix;
          j.get_to(matrix);

          this->itemResponseMatrix = matrix(Eigen::seq(0, matrix.rows() - 1), Eigen::seq(0, matrix.cols() - 3));
          this->rawScores = matrix(Eigen::seq(0, matrix.rows() - 1), Eigen::seq(matrix.cols() - 2, matrix.cols() - 1));
        }
      };
    } // namespace Fixtures
  }   // namespace Tests
} // namespace EquatingRecipes

#endif