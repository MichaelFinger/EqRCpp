#ifndef OPTIMIZATION_FUNCTION_HPP
#define OPTIMIZATION_FUNCTION_HPP

#include <vector>
#include <Eigen/Core>

#include <equating_recipes/irt_model_functions.hpp>
#include <equating_recipes/optimization_function.hpp>
#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>

namespace EquatingRecipes {
  class OptimizationFunction {
  public:
    void configure(const EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      this->oldThetaValues = irtScaleTransformationData.quadratureOldForm.thetaValues;
      this->oldThetaWeights = irtScaleTransformationData.quadratureOldForm.thetaWeights;
      this->newThetaValues = irtScaleTransformationData.quadratureNewForm.thetaValues;
      this->newThetaWeights = irtScaleTransformationData.quadratureNewForm.thetaWeights;

      switch (irtScaleTransformationData.irtScaleTranformationMethod) {
        case EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA:
          this->symmetry = irtScaleTransformationData.haebaraSymmetryOption;
          this->functionStandardization = irtScaleTransformationData.haebaraFunctionStandardization;
          break;
        case EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD:
          this->symmetry = irtScaleTransformationData.stockingLordSymmetryOption;
          this->functionStandardization = irtScaleTransformationData.stockingLordFunctionStandardization;
           break;
        default:
          break;
      }

      this->commonItems = irtScaleTransformationData.commonItems;
    }

    double operator()(const std::vector<double>& x,
                      std::vector<double>& grad) {
      double funcValue = functionValue(x);

      if (!grad.empty()) {
        functionGradient(x, grad);
      }

      return funcValue;
    }

  protected:
    virtual double functionValue(const std::vector<double>& x) = 0;

    virtual void functionGradient(const std::vector<double>& x,
                                  std::vector<double>& grad) = 0;

    Eigen::VectorXd oldThetaValues;
    Eigen::VectorXd oldThetaWeights;
    Eigen::VectorXd newThetaValues;
    Eigen::VectorXd newThetaWeights;
    std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;
    EquatingRecipes::Structures::Symmetry symmetry;
    bool functionStandardization;

    EquatingRecipes::IRTModelFunctions irtModelFunctions;
  };
} // namespace EquatingRecipes

#endif