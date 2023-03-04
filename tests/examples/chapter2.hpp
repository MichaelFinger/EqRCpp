#ifndef TESTS_CHAPTER_2_HPP
#define TESTS_CHAPTER_2_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include "fixtures/actmathfreq.hpp"
#include "fixtures/mondatx.hpp"

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

namespace EquatingRecipes {
  namespace Tests {
    struct Chapter2 {
      void run() {
        EquatingRecipes::JSON::JsonDocument jsonDoc;
        jsonDoc.fromTextFile("../docs/EquatingRecipesExamples/actmathfreq.dat.json");

        EquatingRecipes::Tests::Fixtures::ACTMathFreq actMathFreq;
        actMathFreq.configure(jsonDoc.json);

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
        jsonDoc.fromTextFile("../docs/EquatingRecipesExamples/mondatx.dat.json");

        EquatingRecipes::Tests::Fixtures::MondatX mondatX;
        mondatX.configure(jsonDoc.json);

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
        doc.toTextFile("chpater2_MONDAT_X");
      }
    };
  } // namespace Tests
} // namespace EquatingRecipes
#endif