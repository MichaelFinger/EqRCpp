#ifndef TESTS_EXAMPLES_CHAPTER_2_HPP
#define TESTS_EXAMPLES_CHAPTER_2_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

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

          // ReadFdGet_USTATS("actmathfreq.dat",1,2,0,40,1,'X',&x);
          // Print_USTATS(outf,"ACT Math X",&x);

          EquatingRecipes::Structures::UnivariateStatistics univariateStatistics = EquatingRecipes::Utilities::univariateFromScoreFrequencies(actMathFreq.freqX,
                                                                                                                                              0,
                                                                                                                                              40,
                                                                                                                                              1,
                                                                                                                                              "X");

          nlohmann::json j = univariateStatistics;

          EquatingRecipes::JSON::JsonDocument doc;
          doc.setJson(j);
          doc.toTextFile("chapter2_ACTMath.json");

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
                                                                                                                                 "V");
          j = bivariateStatistics;

          doc.setJson(j);
          doc.toTextFile("chapter2_MONDAT_X.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes
#endif