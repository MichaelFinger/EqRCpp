#ifndef ANALYSES_UNIVARIATE_SCORE_DATASET_HPP
#define ANALYSES_UNIVARIATE_SCORE_DATASET_HPP

#include <string>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Datasets {
    struct UnivariateScoreDataset {
      std::string datasetName;
      std::string variableName;
      std::string id;
      double minimumScore;
      double maximumScore;
      double scoreIncrement;
      Eigen::VectorXd scores;
      Eigen::VectorXd scoreFreqDist;
    };
  }
}

#endif