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
#include <equating_recipes/wrappers/cg_no_smoothing.hpp>
#include <equating_recipes/wrappers/utilities.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

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
          EquatingRecipes::Tests::Examples::Datasets::MondatX mondatX;
          EquatingRecipes::Tests::Examples::Datasets::MondatY mondatY;

          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV = EquatingRecipes::Utilities::bivariateFromScores(mondatX.rawScores,
                                                                                                                                   0,
                                                                                                                                   36,
                                                                                                                                   1,
                                                                                                                                   0,
                                                                                                                                   12,
                                                                                                                                   1,
                                                                                                                                   "X",
                                                                                                                                   "V");

          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV = EquatingRecipes::Utilities::bivariateFromScores(mondatY.rawScores,
                                                                                                                                   0,
                                                                                                                                   36,
                                                                                                                                   1,
                                                                                                                                   0,
                                                                                                                                   12,
                                                                                                                                   1,
                                                                                                                                   "Y",
                                                                                                                                   "V");

          EquatingRecipes::CGEquatingNoSmoothing cgEquatingNoSmoothing;
          EquatingRecipes::Structures::PData pData;
          EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;

          cgEquatingNoSmoothing.run(EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS,
                                    EquatingRecipes::Structures::Method::LINEAR,
                                    EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                    -1.0,
                                    true,
                                    0.0,
                                    0.0,
                                    bivariateStatisticsXV,
                                    bivariateStatisticsYV,
                                    0,
                                    pData,
                                    equatedRawScoreResults);

          nlohmann::json j = nlohmann::json::object();

          j["PData"] = pData;
          j["EquatedRawScoreResults"] = equatedRawScoreResults;

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(j);
          jsonDoc.toTextFile("chapter4.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes
#endif