/* 
  From Source: IRTst.h
  Original Struct: IRTstControl
  Description: 
*/

#ifndef STRUCTURES_IRT_SCALE_TRANSFORMATION_CONTROL_HPP
#define STRUCTURES_IRT_SCALE_TRANSFORMATION_CONTROL_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTScaleTransformationControl {
      size_t numberOfItemsNewForm; /* Numbers of items on the new (X) form */
      size_t numberOfItemsOldForm; /* Numbers of items on the old (Y) form */
      size_t numberOfCommonItems; /* Number of the common items */
      size_t numberOfThetasNewForm; /* Numbers of ability points on the new (X) scale */
      size_t numberOfThetasOldForm; /* Numbers of ability points on the old (Y) scale */
      double minimumRawScoreNewForm;
      double maximumRawScoreNewForm;
      double rawScoreIncrementNewForm;
      double minimumRawScoreOldForm;
      double maximumRawScoreOldForm;
      double rawScoreIncrementOldForm;
      Eigen::VectorXd thetaValuesNewForm;
      Eigen::VectorXd thetaWeightsNewForm;
      Eigen::VectorXd thetaValuesOldForm;
      Eigen::VectorXd thetaWeightsOldForm;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif