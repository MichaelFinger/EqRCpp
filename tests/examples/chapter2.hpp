#ifndef TESTS_EXAMPLES_CHAPTER_2_HPP
#define TESTS_EXAMPLES_CHAPTER_2_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/wrappers/utilities.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
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
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputData;

          inputData.datasetName = "ACTMathFreq";
          inputData.id = "X";
          inputData.variableName = "ACTMathScore";
          inputData.minimumScore = 0;
          inputData.maximumScore = 40;
          inputData.scoreIncrement = 1;
          inputData.scoreFrequencies = actMathFreq.freqX;

          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsACTMath;
          nlohmann::json univariateStatisticsACTMathJson = univariateStatistics(inputData, univariateStatisticsACTMath);

          /* Common-items Nonequivalent Groups Design: 
          Kolen and Brennan (2004) Chapter 4 example (see page 123)*/
          EquatingRecipes::Tests::Examples::Datasets::MondatX mondatX;

          EquatingRecipes::Structures::BivariateStatistics bivariateStatistics = EquatingRecipes::Utilities::bivariateFromScores(mondatX.rawScores,
                                                                                                                                 0,
                                                                                                                                 36,
                                                                                                                                 1,
                                                                                                                                 0,
                                                                                                                                 12,
                                                                                                                                 1,
                                                                                                                                 "X",
                                                                                                                                 "V",
                                                                                                                                 "MondatX",
                                                                                                                                 "RawScoreForm1",
                                                                                                                                 "RawScoreForm2");

          nlohmann::json j = {univariateStatisticsACTMathJson,
                              bivariateStatistics};

          EquatingRecipes::JSON::JsonDocument doc;
          doc.setJson(j);
          doc.toTextFile("chapter2.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes
#endif