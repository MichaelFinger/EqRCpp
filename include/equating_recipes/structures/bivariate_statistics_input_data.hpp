#ifndef STRUCTURES_BIVARIATE_STATISTICS_INPUT_DATA_HPP
#define STRUCTURES_BIVARIATE_STATISTICS_INPUT_DATA_HPP

#include <string>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateStatisticsInputData {
        std::string rowVariableName;
        std::string columnVariableName;

        Eigen::VectorXd rowScores;
        double rowMinimumScore;
        double rowMaximumScore;
        double rowScoreIncrement;

        Eigen::VectorXd columnScores;
        double columnMinimumScore;
        double columnMaximumScore;
        double columnScoreIncrement;

        std::string rowScoreId = "X";
        std::string columnScoreId = "Y";
    };
  }
}

#endif