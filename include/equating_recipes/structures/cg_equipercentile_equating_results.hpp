#ifndef STRUCTURES_CG_EQUIPERCENTILE_EQUATING_RESULTS_HPP
#define STRUCTURES_CG_EQUIPERCENTILE_EQUATING_RESULTS_HPP

#include <optional>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct CGEquipercentileEquatingResults {
      Eigen::VectorXd syntheticPopulationRelativeFreqDistX;              //   fxs[]      rel freq dist for x in syn pop
      Eigen::VectorXd syntheticPopulationRelativeFreqDistY;              //   gys[]      rel freq dist for y in syn pop
      Eigen::VectorXd equatedRawScores;                                  //   eraw[]     equipercentile equated raw scores
      std::optional<double> slope;                                       //   a        slope for Braun-Holland Method
      std::optional<double> intercept;                                   //   b        intercept for Braun-Holland Method
      std::optional<Eigen::VectorXd> braunHollandEquatedRawScores;       //   erawBH[]   Braun-Holland linear equated raw scores
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif