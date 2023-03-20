#ifndef IRT_SCALE_TRANSFORMATION_HPP
#define IRT_SCALE_TRANSFORMATION_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/form_type.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_input.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/haebara_function.hpp>
#include <equating_recipes/irt_model_functions.hpp>
#include <equating_recipes/lbfgs_optimizer.hpp>
#include <equating_recipes/optimization_function.hpp>
#include <equating_recipes/stocking_lord_function.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class IRTScaleTransformation {
  public:
    void runIRTScaleTransformation(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      bool missingQuadrature = irtScaleTransformationData.quadratureNewForm.thetaValues.size() == 0 ||
                               irtScaleTransformationData.quadratureNewForm.thetaWeights.size() == 0 ||
                               irtScaleTransformationData.quadratureOldForm.thetaValues.size() == 0 ||
                               irtScaleTransformationData.quadratureOldForm.thetaWeights.size() == 0;

      if (missingQuadrature) {
        throw std::runtime_error("Haebara or Stocking-Lord results cannot be computed because at least\none of the ability distributions is not present");
      }

      StMeanSigma(irtScaleTransformationData);
      StMeanMean(irtScaleTransformationData);

      if (irtScaleTransformationData.runHaebara) {
        StHaebara(irtScaleTransformationData);
      }

      if (irtScaleTransformationData.runStockingLord) {
        StStockingLord(irtScaleTransformationData);
      }

      getScaleTransformResults(irtScaleTransformationData);
    }

  private:
    /*--------------------------------------------------------------------------------------
      Functionality:
        Use the values of the slope and intercept of the new-to-old transformation
        to convert both the parameter estimates of the new form items and
        the ability points for the new group distribution.
        
        Input:
        
          ItemOutF  The name of a file in which the transformed output for items on
                    the new form is saved; The format of output is, in essence, the
                    same as that of input, so the output file can be read by the function
                    ItemInfoRead without any syntax error.
          DistOutF  The name of a file in which the transformed ability points for the
                    new group distribution are saved along with the original weights;
                    As with ItemOutF, the output is saved using the format of input for
                    the ability distribution, so the file DistOutF can be read by the
                    function ThetaInfoRead without problems.
          slope     The value of A (slope) for a chosen linking method
          intercept The value of B (intercept) for a chosen linking method
          NewItem   A pointer to an array of the ItemSpec structure, which is for the new form
          Handle    A pointer to a variable of the IRTstControl structure
        
        Output:
          ItemOutF  The transformed output file for the new form
          DistOutF  The transformed output file for the new group's ability distribution

          * Note: The original input remains intact for both the items and distribution.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    --------------------------------------------------------------------------------------*/
    void
    getScaleTransformResults(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      double slope = std::numeric_limits<double>::quiet_NaN();
      double intercept = std::numeric_limits<double>::quiet_NaN();

      switch (irtScaleTransformationData.irtScaleTranformationMethod) {
        case EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA:
          slope = irtScaleTransformationData.haebaraSlope.value_or(std::numeric_limits<double>::quiet_NaN());
          intercept = irtScaleTransformationData.haebaraIntercept.value_or(std::numeric_limits<double>::quiet_NaN());
          break;

        case EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_MEAN:
          slope = irtScaleTransformationData.meanMeanSlope;
          intercept = irtScaleTransformationData.meanMeanIntercept;
          break;

        case EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_SIGMA:
          slope = irtScaleTransformationData.meanSigmaSlope;
          intercept = irtScaleTransformationData.meanSigmaIntercept;
          break;

        case EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD:
          slope = irtScaleTransformationData.stockingLordSlope.value_or(std::numeric_limits<double>::quiet_NaN());
          intercept = irtScaleTransformationData.stockingLordIntercept.value_or(std::numeric_limits<double>::quiet_NaN());
          break;

        default:
          break;
      }

      if (slope == std::numeric_limits<double>::quiet_NaN() ||
          intercept == std::numeric_limits<double>::quiet_NaN()) {
        return;
      }

      std::for_each(irtScaleTransformationData.newItems.begin(),
                    irtScaleTransformationData.newItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& newItem) {
                      EquatingRecipes::Structures::IRTScaleTransformationItemResults newItemResult;
                      newItemResult.configure(newItem);

                      switch (newItem.irtModel) {
                        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
                          newItemResult.transformedA(1) = newItem.a(1) / slope;
                          newItemResult.transformedB(1) = newItem.b(1) * slope + intercept;
                          newItemResult.transformedC(1) = newItem.c(1);

                          break;
                        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
                          newItemResult.transformedA(1) = newItem.a(1) / slope;
                          newItemResult.transformedB(Eigen::seq(1, newItem.numberOfCategories - 1)) =
                              newItem.b(Eigen::seq(1, newItem.numberOfCategories - 1)) * slope +
                              Eigen::VectorXd::Constant(newItem.numberOfCategories - 2, intercept);

                          break;
                        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
                          newItemResult.transformedA(1) = newItem.a(1) / slope;
                          newItemResult.transformedB(0) = newItem.b(0) * slope + intercept;
                          newItemResult.transformedB(Eigen::seq(1, newItem.numberOfCategories)) =
                              newItem.d(Eigen::seq(1, newItem.numberOfCategories - 1)) * slope;

                          break;
                        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
                          newItemResult.transformedA(Eigen::seq(0, newItem.numberOfCategories - 1)) =
                              newItem.a(Eigen::seq(0, newItem.numberOfCategories - 1)) / slope;

                          newItemResult.transformedC(Eigen::seq(0, newItem.numberOfCategories - 1)) =
                              newItem.c(Eigen::seq(0, newItem.numberOfCategories - 1)) -
                              (intercept / slope) * newItem.a(Eigen::seq(0, newItem.numberOfCategories - 1));
                          break;
                        default:
                          break;
                      }

                      irtScaleTransformationData.itemResultsNewForm.push_back(newItemResult);
                    });

      irtScaleTransformationData.transformedQuadratureNewForm.thetaValues = (irtScaleTransformationData.quadratureNewForm.thetaValues * slope) +
                                                                            Eigen::VectorXd::Constant(irtScaleTransformationData.quadratureNewForm.thetaValues.size(), intercept);
      irtScaleTransformationData.transformedQuadratureNewForm.thetaWeights = irtScaleTransformationData.quadratureNewForm.thetaWeights;
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Estimate the slope (A) and intercept (B) of the new-to-old transformation
        by using the Haebara method.
        
        Input:
      Handle: A pointer to control solutions
      ComItem: A pointer to an array of the CommonItemSpec structure
      SYM: Symmetric option (old_scale, new_scale, or symmetric)
      FuncStd: Function standardization option (on or off)
      S0: A starting value for the slope
      I0: A starting value of the intercept
        Output:    
      *slope: estimate of A
      *intercept: estimate of B

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void StHaebara(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      double ftol = 0.0000000001;

      std::vector<double> x {irtScaleTransformationData.haebaraSlopeStartingValue.value(),
                             irtScaleTransformationData.haebaraInterceptStartingValue.value()};

      std::shared_ptr<EquatingRecipes::HaebaraFunction> optimizationFunction =
          std::make_shared<EquatingRecipes::HaebaraFunction>();

      optimizationFunction->configure(irtScaleTransformationData);

      EquatingRecipes::LBFGSOptimizer optimizer;
      double functionValue = optimizer.optimize(x,
                                                optimizationFunction,
                                                ftol);

      double slope = x[0];
      double intercept = x[1];

      irtScaleTransformationData.haebaraSlope = slope;
      irtScaleTransformationData.haebaraIntercept = intercept;
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Estimate the slope (A) and intercept (B) of the new-to-old transformation
        by using the Mean/Mean method.
        
        Input:
      Handle: A pointer to an array of the CommonItemSpec structure
      ComItem: A pointer to the CommonItemSpec structure
        Output:    
      *slope: estimate of A
      *intercept: estimate of B

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void StMeanMean(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      int j, k, a_num = 0, b_num = 0, index_a = 0, index_b = 0;
      double new_mu_a = 0.0, old_mu_a = 0.0;
      double new_mu_b = 0.0, old_mu_b = 0.0;
      Eigen::VectorXd na_vec;
      Eigen::VectorXd oa_vec;
      Eigen::VectorXd nb_vec;
      Eigen::VectorXd ob_vec;

      /* counting valid a- and b-parameters by model */

      for (size_t itemIndex = 0; itemIndex < irtScaleTransformationData.commonItems.size(); itemIndex++) {
        EquatingRecipes::Structures::CommonItemSpecification commonItem = irtScaleTransformationData.commonItems[itemIndex];
        EquatingRecipes::Structures::IRTModel irtModel = commonItem.irtModel;

        if (irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC ||
            irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE ||
            irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
          a_num++;
          b_num += (commonItem.numberOfCategories - 1);
        } else {
          a_num += commonItem.numberOfCategories;
          b_num += commonItem.numberOfCategories;
        }
      }

      na_vec.resize(a_num + 1);
      oa_vec.resize(a_num + 1);
      nb_vec.resize(b_num + 1);
      ob_vec.resize(b_num + 1);

      /* copying valid a- and b-parameters by model into nb_vec and ob_vec */
      for (size_t itemIndex = 0; itemIndex < irtScaleTransformationData.commonItems.size(); itemIndex++) {
        EquatingRecipes::Structures::CommonItemSpecification commonItem = irtScaleTransformationData.commonItems[itemIndex];
        EquatingRecipes::Structures::IRTModel irtModel = commonItem.irtModel;

        if (irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC ||
            irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE ||
            irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
          index_a++;
          na_vec(index_a) = commonItem.newA(1);
          oa_vec(index_a) = commonItem.oldA(1);
          for (size_t categoryIndex = 1; categoryIndex < commonItem.numberOfCategories; categoryIndex++) {
            index_b++;
            nb_vec(index_b) = commonItem.newB(categoryIndex);
            ob_vec(index_b) = commonItem.oldB(categoryIndex);
          }
        } else {
          for (size_t categoryIndex = 1; categoryIndex < commonItem.numberOfCategories; categoryIndex++) {
            index_a++;
            index_b++;
            na_vec(index_a) = commonItem.newA(categoryIndex);
            oa_vec(index_a) = commonItem.oldA(categoryIndex);
            nb_vec(index_b) = -1.0 * (commonItem.newC(categoryIndex) / commonItem.newA(categoryIndex));
            ob_vec(index_b) = -1.0 * (commonItem.oldC(categoryIndex) / commonItem.oldA(categoryIndex));
          }
        }
      }

      new_mu_a += na_vec.sum();
      old_mu_a += oa_vec.sum();

      new_mu_b += nb_vec.sum();
      old_mu_b += ob_vec.sum();

      double slope;
      double intercept;

      if (a_num != 0 && b_num != 0) {
        new_mu_a /= a_num;
        old_mu_a /= a_num;
        new_mu_b /= b_num;
        old_mu_b /= b_num;

        slope = new_mu_a / old_mu_a;
        intercept = old_mu_b - slope * new_mu_b;
      } else {
        slope = 1.0;
        intercept = 0.0;
      }

      irtScaleTransformationData.meanMeanSlope = slope;
      irtScaleTransformationData.meanMeanIntercept = intercept;
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Estimate the slope (A) and intercept (B) of the new-to-old transformation
        by using the Mean/Sigma method.
        
        Input:
      Handle: A pointer to an array of the CommonItemSpec structure
      ComItem: A pointer to the CommonItemSpec structure
        Output:    
      *slope: estimate of A
      *intercept: estimate of B

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void StMeanSigma(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      int j, k, b_num = 0, index = 0;
      double new_mu_b = 0.0, old_mu_b = 0.0;
      double new_si_b = 0.0, old_si_b = 0.0;
      Eigen::VectorXd nb_vec;
      Eigen::VectorXd ob_vec;

      /* counting valid b-parameters by model */
      for (size_t itemIndex = 0; itemIndex < irtScaleTransformationData.commonItems.size(); itemIndex++) {
        EquatingRecipes::Structures::CommonItemSpecification commonItem = irtScaleTransformationData.commonItems[itemIndex];
        EquatingRecipes::Structures::IRTModel irtModel = commonItem.irtModel;

        if (irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC ||
            irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE ||
            irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
          b_num += commonItem.numberOfCategories - 1;
        } else {
          b_num += commonItem.numberOfCategories;
        }
      }

      nb_vec.resize(b_num + 1);
      ob_vec.resize(b_num + 1);

      /* copying valid b-parameters by model into nb_vec and ob_vec */
      for (size_t itemIndex = 0; itemIndex < irtScaleTransformationData.commonItems.size(); itemIndex++) {
        EquatingRecipes::Structures::CommonItemSpecification commonItem = irtScaleTransformationData.commonItems[itemIndex];
        EquatingRecipes::Structures::IRTModel irtModel = commonItem.irtModel;

        if (irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC ||
            irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE ||
            irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
          for (size_t categoryIndex = 1; categoryIndex < commonItem.numberOfCategories; categoryIndex++) {
            index++;
            nb_vec(index) = commonItem.newB(categoryIndex);
            ob_vec(index) = commonItem.oldB(categoryIndex);
          }
        } else {
          for (size_t categoryIndex = 0; categoryIndex < commonItem.numberOfCategories; categoryIndex++) {
            index++;
            nb_vec(index) = -1.0 * (commonItem.newC(categoryIndex) / commonItem.newA(categoryIndex));
            ob_vec(index) = -1.0 * (commonItem.oldC(categoryIndex) / commonItem.oldA(categoryIndex));
          }
        }
      }

      new_mu_b += nb_vec.sum();
      old_mu_b += ob_vec.sum();
      new_si_b += nb_vec.cwiseProduct(nb_vec).sum();
      old_si_b += ob_vec.cwiseProduct(ob_vec).sum();

      double slope;
      double intercept;

      if (b_num != 0) {
        new_mu_b /= b_num;
        old_mu_b /= b_num;
        new_si_b = std::sqrt(new_si_b / b_num - new_mu_b * new_mu_b);
        old_si_b = std::sqrt(old_si_b / b_num - old_mu_b * old_mu_b);

        slope = old_si_b / new_si_b;
        intercept = old_mu_b - slope * new_mu_b;
      } else {
        slope = 1.0;
        intercept = 0.0;
      }

      irtScaleTransformationData.meanSigmaSlope = slope;
      irtScaleTransformationData.meanSigmaIntercept = intercept;
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Estimate the slope (A) and intercept (B) of the new-to-old transformation
        by using the Stocking-Lord method.
        
        Input:
      Handle: A pointer to control solutions
      ComItem: A pointer to an array of the CommonItemSpec structure
      SYM: Symmetric option (old_scale, new_scale, symmetric)
      FuncStd: Function standardization option (on or off)
      S0: A starting value for the slope
      I0: A starting value of the intercept
        Output:    
      *slope: estimate of A
      *intercept: estimate of B

      Author: Seonghoon Kim
      Date of last revision 9/25/08
------------------------------------------------------------------------------*/
    void StStockingLord(EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
      std::vector<double> x {irtScaleTransformationData.stockingLordSlopeStartingValue.value(),
                             irtScaleTransformationData.stockingLordInterceptStartingValue.value()};
      double ftol = 0.0000000001;

      std::shared_ptr<EquatingRecipes::StockingLordFunction> optimizationFunction =
          std::make_shared<EquatingRecipes::StockingLordFunction>();

      optimizationFunction->configure(irtScaleTransformationData);

      EquatingRecipes::LBFGSOptimizer optimizer;
      double functionValue = optimizer.optimize(x,
                                                optimizationFunction,
                                                ftol);

      double slope = x[0];
      double intercept = x[1];

      irtScaleTransformationData.stockingLordSlope = slope;
      irtScaleTransformationData.stockingLordIntercept = intercept;
    }
  };
} // namespace EquatingRecipes

#endif