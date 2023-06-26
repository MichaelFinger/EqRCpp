#ifndef TESTS_EXAMPLES_CHAPTER_2_HPP
#define TESTS_EXAMPLES_CHAPTER_2_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

#include "datasets/actmathfreq.hpp"
#include "datasets/mondatx.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter2 {
        void operator()() {
          EquatingRecipes::Tests::Examples::Datasets::ACTMathFreq actMathFreq;

          /* Random Groups Design: Kolen and Brennan (2004)
           Chapter 2 example (see pp. 50-52) */

          EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;
          EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;

          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataACTMath;

          inputDataACTMath.id = "X";
          inputDataACTMath.variableName = "ACTMathScore";
          inputDataACTMath.minimumScore = 0;
          inputDataACTMath.maximumScore = 40;
          inputDataACTMath.scoreIncrement = 1;
          inputDataACTMath.scoreFrequencies = actMathFreq.freqX;

          EquatingRecipes::Analyses::UnivariateStatistics::OutputData univariateStatisticsACTMathOutputData;

          nlohmann::json univariateStatisticsACTMathJson = univariateStatistics(inputDataACTMath,
                                                                                univariateStatisticsACTMathOutputData);

          /* Common-items Nonequivalent Groups Design: 
          Kolen and Brennan (2004) Chapter 4 example (see page 123)*/
          EquatingRecipes::Tests::Examples::Datasets::MondatX mondatX;

          EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataXY;

          inputDataXY.datasetName = "MondatX";
          inputDataXY.rowVariableName = "RawScoreForm1";
          inputDataXY.columnVariableName = "RawScoreForm2";
          inputDataXY.rowScores = mondatX.rawScores.col(0);
          inputDataXY.rowMinimumScore = 0;
          inputDataXY.rowMaximumScore = 36;
          inputDataXY.rowScoreIncrement = 1;
          inputDataXY.columnScores = mondatX.rawScores.col(1);
          inputDataXY.columnMinimumScore = 0;
          inputDataXY.columnMaximumScore = 12;
          inputDataXY.columnScoreIncrement = 1;
          inputDataXY.rowScoreId = "X";
          inputDataXY.columnScoreId = "Y";

          EquatingRecipes::Analyses::BivariateStatistics::OutputData outputDataXY;
          nlohmann::json bivariateStatisticsMondatXJson = bivariateStatistics(inputDataXY,
                                                                              outputDataXY);

          nlohmann::json j = nlohmann::json::array();
          j.push_back(univariateStatisticsACTMathJson);
          j.push_back(bivariateStatisticsMondatXJson);

          EquatingRecipes::JSON::JsonDocument doc;
          doc.setJson(j);
          doc.toTextFile("chapter2.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes
#endif