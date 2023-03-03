#ifndef TESTS_UNIVARIATE_STATISTICS_HPP
#define TESTS_UNIVARIATE_STATISTICS_HPP

#include <filesystem>
#include <iostream>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include "fixtures/lsat6.hpp"
#include "fixtures/actmathfreq.hpp"

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

namespace Tests {
  struct UnivariateStatistics {
    void runLSAT6() {
      Eigen::VectorXd lsat6FreqDist = EquatingRecipes::Tests::Fixtures::LSAT6::rawScoreFrequencyDistribution();

      std::for_each(lsat6FreqDist.begin(),
                    lsat6FreqDist.end(),
                    [&](const double& scoreFreq) {
                      std::cout << scoreFreq << "\n";
                    });

      Eigen::VectorXd lsat6RelFreqDist = EquatingRecipes::Tests::Fixtures::LSAT6::rawScoreRelativeFrequencyDistribution();

      std::for_each(lsat6RelFreqDist.begin(),
                    lsat6RelFreqDist.end(),
                    [&](const double& scoreRelFreq) {
                      std::cout << scoreRelFreq << "\n";
                    });

      Eigen::VectorXi lsat6CumFreqDist(lsat6FreqDist.size());
      lsat6CumFreqDist(0) = lsat6FreqDist(0);
      for (size_t index = 1; index < lsat6FreqDist.size(); index++) {
        lsat6CumFreqDist(index) = lsat6CumFreqDist(index - 1) + lsat6FreqDist(index);
      }

      std::for_each(lsat6CumFreqDist.begin(),
                    lsat6CumFreqDist.end(),
                    [&](const int& lsat6CumFreq) {
                      std::cout << lsat6CumFreq << "\n";
                    });

      Eigen::VectorXd lsat6CumRelFreqDist = lsat6CumFreqDist.cast<double>() / static_cast<double>(lsat6FreqDist.sum());

      std::for_each(lsat6CumRelFreqDist.begin(),
                    lsat6CumRelFreqDist.end(),
                    [&](const double& lsat6CumRelFreq) {
                      std::cout << lsat6CumRelFreq << "\n";
                    });

      double minimumScore = 0;
      double maxixmumScore = 5;
      double scoreIncrement = 1;

      double maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(maxixmumScore,
                                                                                 minimumScore,
                                                                                 scoreIncrement);

      EquatingRecipes::Structures::UnivariateStatistics univariateStatistics =
          EquatingRecipes::Utilities::univariateFromScoreFrequencies(lsat6FreqDist,
                                                                    0,
                                                                    5,
                                                                    1,
                                                                    "X");

      nlohmann::json j = univariateStatistics;

      EquatingRecipes::JSON::JsonDocument jsonDoc;
      jsonDoc.setJson(j);
      jsonDoc.toTextFile("tmp.json");
    }

    void runACTMath() {
      std::filesystem::path datasetFolder = std::filesystem::relative("../../ER for distribution 9-10-09/Examples/");

      EquatingRecipes::Tests::Fixtures::ACTMathFreq actMathFreq;
      // actMathFreq.import(datasetFolder.string());
      actMathFreq.import("/Users/michaelfinger/Developer/EquatingRecipes/ER for distribution 9-10-09/Examples/");

      std::cout << actMathFreq.data.rawScores << "\n";
      std::cout << actMathFreq.data.freqX << "\n";
      std::cout << actMathFreq.data.freqY << "\n";
    }
  };
} // namespace Tests
#endif