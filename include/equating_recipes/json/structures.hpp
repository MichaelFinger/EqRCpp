#ifndef JSON_STRUCTURES_HPP
#define JSON_STRUCTURES_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <nlohmann/json.hpp>
#include <equating_recipes/structures/all_structures.hpp>

namespace EquatingRecipes {
  namespace Structures {
    void to_json(nlohmann::json& j, const Eigen::VectorXi& rec) {
      std::vector<int> values;

      if (rec.size() >= 1) {
        values.resize(rec.size());

        for (size_t index = 0; index < rec.size(); index++) {
          values[index] = rec(index);
        }
      }

      j = values;
    }

    void to_json(nlohmann::json& j, const Eigen::VectorXd& rec) {
      std::vector<double> values;

      if (rec.size() >= 1) {
        values.resize(rec.size());

        for (size_t index = 0; index < rec.size(); index++) {
          values[index] = rec(index);
        }
      }

      j = values;
    }

    void to_json(nlohmann::json& j, const Eigen::MatrixXi& rec) {
     std::vector<std::vector<int>> data;
      data.resize(rec.rows());

      for (size_t rowIndex = 0; rowIndex < rec.rows(); rowIndex++) {
        std::vector<int> rowValues;
        rowValues.resize(rec.cols());

        for (size_t columnIndex = 0; columnIndex < rec.cols(); columnIndex++) {
          rowValues[columnIndex] = rec(rowIndex, columnIndex);
        }

        data[rowIndex] = rowValues;
      }

      j = data;
    }

    void to_json(nlohmann::json& j, const Eigen::MatrixXd& rec) {
      std::vector<std::vector<double>> data;
      data.resize(rec.rows());

      for (size_t rowIndex = 0; rowIndex < rec.rows(); rowIndex++) {
        std::vector<double> rowValues;
        rowValues.resize(rec.cols());

        for (size_t columnIndex = 0; columnIndex < rec.cols(); columnIndex++) {
          rowValues[columnIndex] = rec(rowIndex, columnIndex);
        }

        data[rowIndex] = rowValues;
      }

      j = data;
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::UnivariateStatistics& rec) {
      j = nlohmann::json {{"id", rec.id},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"minimumScore", rec.minimumScore},
                          {"maximumScore", rec.maximumScore},
                          {"frequencyDisributiontMinimumScoreObserved", rec.freqDistMinimumScore},
                          {"frequencyDisributiontMaximumScoreObserved", rec.freqDistMaximumScore},
                          {"scoreIncrement", rec.scoreIncrement},
                          {"numberOfScores", rec.numberOfScores},
                          {"frequencyDistribution", rec.freqDist},
                          {"frequencyDistributionDouble", rec.freqDistDouble},
                          {"cumulativeFrequencyDistribution", rec.cumulativeFreqDist},
                          {"relativeFrequencyDistribution", rec.relativeFreqDist},
                          {"cumulativeRelativeFrequencyDistribution", rec.cumulativeRelativeFreqDist},
                          {"percentileRankDistribution", rec.percentileRankDist},
                          {"moments", rec.momentValues}};
    }
  } // namespace Structures
} // namespace EquatingRecipes

#endif