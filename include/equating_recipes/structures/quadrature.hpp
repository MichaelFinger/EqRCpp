#ifndef STRUCTURES_QUADRATURE_HPP
#define STRUCTURES_QUADRATURE_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct Quadrature {
      Eigen::VectorXd thetaValues;
      Eigen::VectorXd thetaWeights;
    };
  }
}

#endif