#ifndef LBFGS_OPTIMIZER_HPP
#define LBFGS_OPTIMIZER_HPP

#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include <nlopt.hpp>

#include <equating_recipes/optimization_function.hpp>
#include <equating_recipes/structures/optimization_results.hpp>

#include <iostream>

namespace EquatingRecipes {
  class LBFGSOptimizer {
  public:
    static double wrap(const std::vector<double>& x, std::vector<double>& grad, void* data) {
      double func = (*reinterpret_cast<EquatingRecipes::OptimizationFunction*>(data))(x, grad);

      std::cout.precision(10);
      std::cout << "func = " << func
                << ", slope = " << x[0]
                << ", intercept = " << x[1]
                << " grad = (" << grad[0] << ", " << grad[1] << ")\n";
      std::cout.precision();

      return func;
    }

    EquatingRecipes::Structures::OptimizationResults optimize(std::vector<double>& x,
                                                              std::shared_ptr<EquatingRecipes::OptimizationFunction> optimizationFunction,
                                                              const int& maximumNumberOfIterations,
                                                              const std::optional<double>& maximumAbsoluteChangeInFunctionValue,
                                                              const std::optional<double>& maximumRelativeChangeInFunctionValue,
                                                              const std::optional<double>& maximumAbsoluteChangeInParameterValues,
                                                              const std::optional<double>& maximumRelativeChangeInParameterValues) {
      nlopt::opt opt(nlopt::LD_LBFGS, x.size());
    
      opt.set_min_objective(EquatingRecipes::LBFGSOptimizer::wrap, &(*optimizationFunction));

      opt.set_maxeval(maximumNumberOfIterations);

      if (maximumAbsoluteChangeInFunctionValue.has_value()) {
        double value = maximumAbsoluteChangeInFunctionValue.value();
        opt.set_ftol_abs(value);
      }

      if (maximumRelativeChangeInFunctionValue.has_value()) {
        double value = maximumRelativeChangeInFunctionValue.value();
        opt.set_ftol_rel(value);
      }

      opt.set_x_weights(1.0);

      if (maximumAbsoluteChangeInParameterValues.has_value()) {
        double value = maximumAbsoluteChangeInParameterValues.value();
        opt.set_xtol_abs(value);
      }

      if (maximumRelativeChangeInParameterValues.has_value()) {
        double value = maximumRelativeChangeInParameterValues.value();
        opt.set_xtol_rel(value);
      }

      EquatingRecipes::Structures::OptimizationResults results;

      std::for_each(x.begin(),
                    x.end(),
                    [&](const double& value) {
                      results.parameterStartingValues.push_back(value);
                    });

      double minf;

      // try {
        nlopt::result result = opt.optimize(x, minf);

        results.functionValue = minf;

        std::for_each(x.begin(),
                      x.end(),
                      [&](const double& value) {
                        results.parameterEstimates.push_back(value);
                      });

        std::vector<double> grad;
        grad.resize(x.size());
        optimizationFunction->operator()(results.parameterEstimates,
                                         grad);

        std::for_each(grad.begin(),
                      grad.end(),
                      [&](const double& value) {
                        results.gradientVector.push_back(value);
                      });

        std::vector<double> xtolAbs = opt.get_xtol_abs();
        auto iter = std::max_element(xtolAbs.begin(), xtolAbs.end());

        results.numberOfIterations = opt.get_numevals();
        results.maximumAbsoluteChangeInFunctionValue = opt.get_ftol_abs();
        results.maximumRelativeChangeInFunctionValue = opt.get_ftol_rel();
        results.maximumAbsoluteChangeInParameterValues = *iter;
        results.maximumRelativeChangeInParameterValues = opt.get_xtol_rel();

        results.resultCode = getOptimizationResultCode(result);
      // } catch (const std::exception& e) {
      //   std::cerr << e.what() << '\n';
      // }

      return results;
    }

  private:
    EquatingRecipes::Structures::OptimizationResultCode getOptimizationResultCode(const nlopt::result& nloptResult) {
      EquatingRecipes::Structures::OptimizationResultCode resultCode;

      switch (nloptResult) {
        case nlopt::result::FAILURE:
        case nlopt::result::INVALID_ARGS:
        case nlopt::result::OUT_OF_MEMORY:
        case nlopt::result::ROUNDOFF_LIMITED:
        case nlopt::result::FORCED_STOP:
        case nlopt::result::NUM_FAILURES:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::FAILED;
          break;
        case nlopt::result::SUCCESS:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::SUCCESS;
          break;
        case nlopt::result::STOPVAL_REACHED:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::FUNCTION_STOP_VALUE;
          break;
        case nlopt::result::FTOL_REACHED:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::FUNCTION_TOLERANCE;
          break;
        case nlopt::result::XTOL_REACHED:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::PARAMETER_TOLERANCE;
          break;
        case nlopt::result::MAXEVAL_REACHED:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::MAXIMUM_NUMBER_OF_ITERATIONS;
          break;
        case nlopt::result::MAXTIME_REACHED:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::MAXIMUM_TIME;
          break;
        case nlopt::result::NUM_RESULTS:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::SUCCESS;
        default:
          resultCode = EquatingRecipes::Structures::OptimizationResultCode::FAILED;
          break;
      }

      return resultCode;
    }
  };
} // namespace EquatingRecipes

#endif