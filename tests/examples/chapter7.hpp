#ifndef TESTS_EXAMPLES_CHAPTER_7_HPP
#define TESTS_EXAMPLES_CHAPTER_7_HPP

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>

#include <equating_recipes/analyses/linear_equating_random_groups.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/analytic_standard_errors.hpp>

#include "datasets/actmathfreq.hpp"
#include "datasets/yctmath.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter7 {
        void operator()() {
          /* Random Groups Design: 
              Kolen and Brennan (2004): Chapter 2 example: 
              Standard Errors of Equipercentile equating (see pp. 81) */
          EquatingRecipes::Tests::Examples::Datasets::ACTMathFreq actMathFreq;

          EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataX;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataY;
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;

          inputDataX.datasetName = "ACT Math";
          inputDataX.id = "X";
          inputDataX.maximumScore = 40;
          inputDataX.minimumScore = 0;
          inputDataX.scoreIncrement = 1;
          inputDataX.scoreFrequencies = actMathFreq.freqX;
          inputDataX.title = "ACT Math X---Univariate Statistics";
          inputDataX.variableName = "X";

          nlohmann::json univariateStatisticsXJson = univariateStatistics(inputDataX,
                                                                          univariateStatisticsX);

          inputDataY.datasetName = "ACT Math";
          inputDataY.id = "Y";
          inputDataY.maximumScore = 40;
          inputDataY.minimumScore = 0;
          inputDataY.scoreIncrement = 1;
          inputDataY.scoreFrequencies = actMathFreq.freqY;
          inputDataY.title = "ACT Math Y---Univariate Statistics";
          inputDataY.variableName = "Y";

          nlohmann::json univariateStatisticsYJson = univariateStatistics(inputDataY,
                                                                          univariateStatisticsY);

          EquatingRecipes::Analyses::AnalyticStandardErrors analyticStandardErrors;
          EquatingRecipes::Analyses::AnalyticStandardErrors::InputData inputData;

          inputData.title = "ACT Math - Equating Standard Errors";
          inputData.datasetName = "ACT Math";
          inputData.univariateStatisticsX = univariateStatisticsX;
          inputData.univariateStatisticsY = univariateStatisticsY;

          Eigen::VectorXd stdErrors;

          nlohmann::json analyticStandardErrorsJson = analyticStandardErrors(inputData,
                                                                             stdErrors);

          nlohmann::json j = {univariateStatisticsXJson,
                              univariateStatisticsYJson,
                              analyticStandardErrorsJson};

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(j);
          jsonDoc.toTextFile("chapter7.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif