#ifndef ANALYSES_UNIVARIATE_SCORE_DATASET_HPP
#define ANALYSES_UNIVARIATE_SCORE_DATASET_HPP

#include <string>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Datasets {
    struct BivariateScoreDataset {
      std::string datasetName;
      std::string rowVariableName;
      std::string columnVariableName;
      std::string rowVariableId;
      std::string columnVariableId;
      double minimumScoreX;
      double maximumScoreX;
      double scoreIncrementX;
      double minimumScoreY;
      double maximumScoreY;
      double scoreIncrementY;
      Eigen::MatrixXd scores;
    };
  }
}

#endif