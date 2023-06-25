#ifndef STRUCTURES_QUADRATURE_HPP
#define STRUCTURES_QUADRATURE_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct Quadrature {
      Eigen::VectorXd thetaValues;
      Eigen::VectorXd thetaWeights;

      bool empty() {
        bool isEmpty = (thetaValues.size() == 0 || thetaWeights.size() == 0);

        return isEmpty;
      }
    };
  }
}

#endif