#include <iostream>

#include "chapter2.hpp"
#include "chapter3.hpp"
// #include "chapter4.hpp"
// #include "chapter5.hpp"
// #include "chapter6.hpp"
// #include "chapter7.hpp"
// #include "chapter8.hpp"
// #include "irt_scale_transformation.hpp"

#include <equating_recipes/analyses/analyses.hpp>

#include <equating_recipes/implementation/observed_score_distribution.hpp>

int main(int argc, char const* argv[]) {
  // EquatingRecipes::Tests::Examples::Chapter2 ch2;
  // ch2();

  // EquatingRecipes::Tests::Examples::Chapter3 ch3;
  // ch3();

  // EquatingRecipes::Tests::Examples::Chapter4 ch4;
  // ch4();

  // EquatingRecipes::Tests::Examples::Chapter5 ch5;
  // ch5();

  // EquatingRecipes::Tests::Examples::Chapter6 ch6;
  // ch6();

  // EquatingRecipes::Tests::Examples::Chapter7 ch7;
  // ch7();

  // EquatingRecipes::Tests::Examples::Chapter8 ch8;
  // ch8();

  // EquatingRecipes::Tests::Examples::IRTScaleTransformation irtScaleTransformation;
  // irtScaleTransformation();

  EquatingRecipes::Implementation::ObservedScoreDistribution observedScoreDistribution;
  std::cout << observedScoreDistribution.probResp3P(1.5,
                                                    1.0,
                                                    0.2,
                                                    -1,
                                                    1.702,
                                                    EquatingRecipes::Implementation::ObservedScoreDistribution::Density::LOGISTIC,
                                                    false)
            << "\n";

  Eigen::VectorXd thresholds(5);
  thresholds(0) = -1;
  thresholds(1) = 0.05;
  thresholds(2) = 0.7;
  thresholds(3) = 1.1;
  thresholds(4) = 2.01;

  std::cout << observedScoreDistribution.probRespGR(1.5,
                                                    thresholds,
                                                    -1,
                                                    1.702,
                                                    EquatingRecipes::Implementation::ObservedScoreDistribution::Density::LOGISTIC,
                                                    false)
            << "\n";

  return 0;
}
