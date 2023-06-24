#ifndef IMPLEMENTATION_QUADRATURE_HPP
#define IMPLEMENTATION_QUADRATURE_HPP

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