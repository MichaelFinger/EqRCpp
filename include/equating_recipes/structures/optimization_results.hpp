#ifndef STRUCTURES_OPTIMIZATION_RESULTS_HPP
#define STRUCTURES_OPTIMIZATION_RESULTS_HPP

#include <optional>
#include <vector>

namespace EquatingRecipes {
  namespace Structures {
    enum class OptimizationResultCode {
      SUCCESS,
      FAILED,
      FUNCTION_STOP_VALUE,
      FUNCTION_TOLERANCE,
      PARAMETER_TOLERANCE,
      MAXIMUM_NUMBER_OF_ITERATIONS,
      MAXIMUM_TIME
    };

    struct OptimizationResults {
      std::vector<double> parameterStartingValues;
      std::vector<double> parameterEstimates;
      double functionValue;
      std::vector<double> gradientVector;
      OptimizationResultCode resultCode;
      unsigned long numberOfIterations;
      double maximumAbsoluteChangeInFunctionValue;
      double maximumRelativeChangeInFunctionValue;
      double maximumAbsoluteChangeInParameterValues;
      double maximumRelativeChangeInParameterValues;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif