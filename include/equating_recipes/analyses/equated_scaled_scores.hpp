#ifndef ANALYSES_EQUATED_SCALED_SCORES_HPP
#define ANALYSES_EQUATED_SCALED_SCORES_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/implementation/utilities.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct EquatedScaledScores {
      struct InputData {
        std::string title;
        std::string datasetName;
        double lowestObservableEquatedRawScore;
        double highestObservableEquatedRawScore;
        double scoreIncrementEquatedRawScore;
        double lowestObservableScaledScore;
        double highestObservableScaledScore;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;

        size_t roundToNumberOfDecimalPlaces = 1;
      };

      struct OutputData {
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::EquatedScaledScores::InputData& inputData,
                                EquatingRecipes::Analyses::EquatedScaledScores::OutputData& outputData) {
        outputData.pData = inputData.pData;

        EquatingRecipes::Implementation::Utilities::runEquatedScaledScores(outputData.pData,
                                                                           inputData.equatedRawScoreResults,
                                                                           inputData.lowestObservableEquatedRawScore,
                                                                           inputData.highestObservableEquatedRawScore,
                                                                           inputData.scoreIncrementEquatedRawScore,
                                                                           inputData.rawToScaledScoreTable,
                                                                           inputData.roundToNumberOfDecimalPlaces,
                                                                           inputData.lowestObservableScaledScore,
                                                                           inputData.highestObservableScaledScore,
                                                                           outputData.equatedScaledScoreResults);

        nlohmann::json results = nlohmann::json::object();
        results["DatasetName"] = inputData.datasetName;
        results["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults;

        nlohmann::json j = {{"analysis_title", inputData.title},
                            {"analysis_type", "equated_scaled_scores"},
                            {"analysis_results", results}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif