#ifndef TESTS_EXAPLES_DATASETS_MONDATX_HPP
#define TESTS_EXAPLES_DATASETS_MONDATX_HPP

#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class MondatX {
        public:
          std::string datasetName = "Mondat X";
          Eigen::MatrixXd itemResponseMatrix;
          Eigen::MatrixXd rawScores;

          MondatX() {
            this->configure();
          }

        private:
          void configure() {
            EquatingRecipes::JSON::JsonDocument jsonDoc;
            jsonDoc.fromTextFile("./resources/json/datasets/mondatx.dat.json");
            nlohmann::json j = jsonDoc.json;

            Eigen::MatrixXd matrix;
            j.get_to(matrix);

            this->itemResponseMatrix = matrix(Eigen::seq(0, matrix.rows() - 1), Eigen::seq(0, matrix.cols() - 3));
            this->rawScores = matrix(Eigen::seq(0, matrix.rows() - 1), Eigen::seq(matrix.cols() - 2, matrix.cols() - 1));
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes

#endif