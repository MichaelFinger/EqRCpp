#ifndef LBFGS_OPTIMIZER_HPP
#define LBFGS_OPTIMIZER_HPP

#include <memory>
#include <optional>
#include <vector>

#include <nlopt.hpp>

#include <equating_recipes/optimization_function.hpp>

namespace EquatingRecipes {
  class LBFGSOptimizer {
  public:
    static double wrap(const std::vector<double>& x, std::vector<double>& grad, void* data) {
      return (*reinterpret_cast<EquatingRecipes::OptimizationFunction*>(data))(x, grad);
    }

    double optimize(const std::vector<double>& x,
                    std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction,
                    const unsigned long& maximumNumberOfIterations,
                    const std::optional<double>& maximumAbsoluteChangeInFunctionValue,
                    const std::optional<double>& maximumRelativeChangeInFunctionValue,
                    const std::optional<double>& maximumAbsoluteChangeInParameterValues,
                    const std::optional<double>& maximumRelativeChangeInParameterValues) {
      nlopt::opt opt(nlopt::LD_LBFGS, x.size());

      opt.set_min_objective(EquatingRecipes::LBFGSOptimizer::wrap, &(*optimizationFunction));

      opt.set_maxeval(maximumNumberOfIterations);

      if (maximumAbsoluteChangeInFunctionValue.has_value()) {
        opt.set_ftol_abs(maximumAbsoluteChangeInFunctionValue.value());
      }

      if (maximumRelativeChangeInFunctionValue.has_value()) {
        opt.set_ftol_rel(maximumRelativeChangeInFunctionValue.value());
      }

      if (maximumAbsoluteChangeInParameterValues.has_value()) {
        opt.set_xtol_abs(maximumAbsoluteChangeInParameterValues.value());
      }

      if (maximumRelativeChangeInParameterValues.has_value()) {
        opt.set_xtol_rel(maximumRelativeChangeInParameterValues.value());
      }

      double minf;

      try {
         opt.optimize(x);
      } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
      }

      return minf;
    }
  };
} // namespace EquatingRecipes

#endif