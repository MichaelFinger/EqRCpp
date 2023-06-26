#ifndef TESTS_EXAMPLES_CHAPTER_4_HPP
#define TESTS_EXAMPLES_CHAPTER_4_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

#include <equating_recipes/analyses/bivariate_statistics.hpp>
// #include <equating_recipes/analyses/common_item_nonequivalent_groups.hpp>

#include "datasets/mondatx.hpp"
#include "datasets/mondaty.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter4 {
        void operator()() {
          /* Common-item Nonequivalent Groups Design: 
            Kolen and Brennan (2004) Chapter 4 example:
            Linear equating  pp. 121-124 */

          // convertFtoW("mondatx.dat",2,fieldsACT,"mondatx-temp");
          // ReadRawGet_BSTATS("mondatx-temp",1,2,0,36,1,0,12,1,'X','V',&xv);

          // convertFtoW("mondaty.dat",2,fieldsACT,"mondaty-temp");
          // ReadRawGet_BSTATS("mondaty-temp",1,2,0,36,1,0,12,1,'Y','V',&yv);

          // Wrapper_CN('C','L','N',-1,1,0,0,&xv,&yv,0,&pdCLN,&rCLN);
          // Print_CN(outf,"Chapter 4: proportional wts",&pdCLN,&rCLN);

          /* Common-items Nonequivalent Groups Design: 
          Kolen and Brennan (2004) Chapter 4 example (see page 123)*/
          
          // EquatingRecipes::Tests::Examples::Datasets::MondatX mondatX;
          // EquatingRecipes::Tests::Examples::Datasets::MondatY mondatY;

          // EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;
          // EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataXV;
          // EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataYV;

          // EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
          // EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;

          // inputDataXV.datasetName = "MondatX";
          // inputDataXV.rowVariableName = "X";
          // inputDataXV.columnVariableName = "V";
          // inputDataXV.rowScores = mondatX.rawScores.col(0);
          // inputDataXV.rowMinimumScore = 0;
          // inputDataXV.rowMaximumScore = 36;
          // inputDataXV.rowScoreIncrement = 1;
          // inputDataXV.columnScores = mondatX.rawScores.col(1);
          // inputDataXV.columnMinimumScore = 0;
          // inputDataXV.columnMaximumScore = 12;
          // inputDataXV.columnScoreIncrement = 1;
          // inputDataXV.rowScoreId = "X";
          // inputDataXV.columnScoreId = "V";

          // inputDataYV.datasetName = "MondatY";
          // inputDataYV.rowVariableName = "Y";
          // inputDataYV.columnVariableName = "V";
          // inputDataYV.rowScores = mondatY.rawScores.col(0);
          // inputDataYV.rowMinimumScore = 0;
          // inputDataYV.rowMaximumScore = 36;
          // inputDataYV.rowScoreIncrement = 1;
          // inputDataYV.columnScores = mondatY.rawScores.col(1);
          // inputDataYV.columnMinimumScore = 0;
          // inputDataYV.columnMaximumScore = 12;
          // inputDataYV.columnScoreIncrement = 1;
          // inputDataYV.rowScoreId = "Y";
          // inputDataYV.columnScoreId = "V";

          // nlohmann::json bivariateStatisticsXVJson = bivariateStatistics(inputDataXV,
          //                                                                bivariateStatisticsXV);

          // nlohmann::json bivariateStatisticsYVJson = bivariateStatistics(inputDataYV,
          //                                                                bivariateStatisticsYV);

          // EquatingRecipes::Analyses::CommonItemNonequivalentGroups::InputData inputData;
          // EquatingRecipes::Analyses::CommonItemNonequivalentGroups::OutputData outputData;

          // inputData.datasetName = "mondat";
          // inputData.xVariableName = "X";
          // inputData.yVariableName = "Y";
          // inputData.bivariateStatisticsXV = bivariateStatisticsXV;
          // inputData.bivariateStatisticsYV = bivariateStatisticsYV;
          // inputData.design = EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS;
          // inputData.method = EquatingRecipes::Structures::Method::LINEAR;
          // inputData.smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;

          // EquatingRecipes::Analyses::CommonItemNonequivalentGroups cgEquatingNoSmoothing;
          // nlohmann::json commonItemNonequivalentGroupsJson = cgEquatingNoSmoothing(inputData,
          //                                                                          outputData);

          // nlohmann::json j = commonItemNonequivalentGroupsJson;

          // EquatingRecipes::JSON::JsonDocument jsonDoc;
          // jsonDoc.setJson(j);
          // jsonDoc.toTextFile("chapter4.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes
#endif