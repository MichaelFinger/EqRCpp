/* 
  From Source: IRTeq.h 
  Original Struct: RawFitDist
  Description: Structure that stores IRT fitted distributions
*/

#ifndef IMPLEMENTATION_IRT_FITTED_DISTRIBUTION_HPP
#define IMPLEMENTATION_IRT_FITTED_DISTRIBUTION_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTFittedDistribution {
      size_t numberOfRawScoreCategories;                                      /* Number of raw score categories */
      Eigen::VectorXd rawScores;                                              /* Raw scores; zero-offset */
      Eigen::VectorXd fittedDistributionNewGroup;                             /* Fitted distribution for the new group; zero-offset */
      Eigen::VectorXd fittedDistributionOldGroup;                             /* Fitted distribution for the old group; zero-offset */
      Eigen::VectorXd fittedDistributionSyntheticGroup;                       /* Fitted distribution for the synthetic group; zero-offset */
      Eigen::VectorXd momentsFittedDistributionNewGroup;                      /* moments for fitted distribtion for new group */
      Eigen::VectorXd momentsFittedDistributionOldGroup;                      /* moments for fitted distribtion for old group */
      Eigen::VectorXd momentsFittedDistributionSyntheticGroup;                /* moments for fitted distribtion for synthetic group */      
    };
  }
}

#endif