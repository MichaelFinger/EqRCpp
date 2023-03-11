#ifndef LBFGS_OPTIMIZER_HPP
#define LBFGS_OPTIMIZER_HPP

#include <memory>
#include <vector>

#include <nlopt.hpp>

#include <equating_recipes/wrappers/optimization_function.hpp>

namespace EquatingRecipes {
  class LBFGSOptimizer {
  public:
    static double wrap(const std::vector<double>& x, std::vector<double>& grad, void* data) {
      return (*reinterpret_cast<EquatingRecipes::OptimizationFunction*>(data))(x, grad);
    }

    double optimize(const std::vector<double>& x,
                           std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction,
                           const double& ftol) {
      nlopt::opt opt(nlopt::LD_LBFGS, x.size());

      opt.set_min_objective(EquatingRecipes::LBFGSOptimizer::wrap, &(*optimizationFunction));

      opt.set_ftol_rel(ftol);

      double minf;

      opt.optimize(x);

      return minf;
    }
  };
} // namespace EquatingRecipes

#endif