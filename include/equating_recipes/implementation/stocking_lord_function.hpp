#ifndef EQUATING_RECIPES_STOCKING_LORD_FUNCTION_HPP
#define EQUATING_RECIPES_STOCKING_LORD_FUNCTION_HPP

#include <cmath>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/implementation/optimization_function.hpp>

namespace EquatingRecipes {
  namespace Implementation {
    class StockingLordFunction : public EquatingRecipes::Implementation::OptimizationFunction {
    protected:
      /*------------------------------------------------------------------------------
      Functionality:
        calculate the value of the Stocking-Lord criterion function.
        
        Input:
      x: the point of slope and intercept
        Output:    
      the value of the criterion function

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double functionValue(const std::vector<double>& x) override {
        // int i, j, k;
        // double w1_sum = 0.0, w2_sum = 0.0;
        // double sym_f1, sym_f2;
        // double theta, th_weight;
        // double p_origi, p_trans, Ujk;
        // double func1, func2;
        double func1_sum = 0.0;
        double func2_sum = 0.0;

        /* F1: old scale */
        for (size_t quadratureIndex = 0; quadratureIndex < oldThetaValues.size(); quadratureIndex++) {
          double theta = oldThetaValues(quadratureIndex);
          double th_weight = oldThetaWeights(quadratureIndex);

          double func1 = 0.0;
          for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
            EquatingRecipes::Structures::CommonItemSpecification commonItem = commonItems[itemIndex];

            for (size_t categoryIndex = 1; categoryIndex <= commonItem.numberOfCategories; categoryIndex++) {
              double uJK = commonItem.scoringFunctionValues(categoryIndex);
              double p_origi = irtModelFunctions.probOld(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         true,
                                                         x[0],
                                                         x[1]);
              double p_trans = irtModelFunctions.probOld(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         false,
                                                         x[0],
                                                         x[1]);
              func1 += uJK * (p_origi - p_trans);
            }
          }

          func1_sum += std::pow(func1, 2) * th_weight;
        }

        double w1_sum = oldThetaWeights.sum();

        /* F2: new scale */

        for (size_t quadratureIndex = 0; quadratureIndex < oldThetaValues.size(); quadratureIndex++) {
          double theta = newThetaValues(quadratureIndex);
          double th_weight = newThetaWeights(quadratureIndex);

          double func2 = 0.0;

          for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
            EquatingRecipes::Structures::CommonItemSpecification commonItem = commonItems[itemIndex];

            for (size_t categoryIndex = 1; categoryIndex <= commonItem.numberOfCategories; categoryIndex++) {
              double uJK = commonItem.scoringFunctionValues(categoryIndex);
              double p_orig = irtModelFunctions.probNew(commonItem,
                                                        categoryIndex,
                                                        theta,
                                                        true,
                                                        x[0],
                                                        x[1]);

              double p_trans = irtModelFunctions.probNew(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         false,
                                                         x[0],
                                                         x[1]);

              func2 += uJK * (p_orig - p_trans);
            }
          }

          func2_sum += std::pow(func2, 2) * th_weight;
        }

        double w2_sum = newThetaWeights.sum();

        /* symmetric or non-symmetric setting*/
        double sym_f1;
        double sym_f2;

        if (this->symmetry == EquatingRecipes::Structures::Symmetry::SYMMETRIC) {
          sym_f1 = 1.0;
          sym_f2 = 1.0;
        } else {
          if (this->symmetry == EquatingRecipes::Structures::Symmetry::OLD_SCALE) {
            sym_f1 = 1.0;
            sym_f2 = 0.0;
          } else {
            sym_f1 = 0.0;
            sym_f2 = 1.0;
          }
        }

        /* function standardization */
        double functionValue;
        if (this->functionStandardization) {
          functionValue = (sym_f1 * func1_sum / w1_sum + sym_f2 * func2_sum / w2_sum);
        } else {
          functionValue = (sym_f1 * func1_sum + sym_f2 * func2_sum);
        }

        return functionValue;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        calculate the partial derivatives of the Stocking-Lord criterion function
        with respect to slope (S or A) and intercept (I or B).
        
        Input:
      x: the point of slope and intercept
        Output:    
      the gradient of the criterion function w.r.t. A and B.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      void functionGradient(const std::vector<double>& x,
                            std::vector<double>& grad) override {
        // int i, j, k;
        // double w1_sum = 0.0, w2_sum = 0.0;
        // double sym_f1, sym_f2;

        // double theta, th_weight;
        // double p_origi, p_trans, ps, pi, Ujk;
        // double func1, func2, ps_sum, pi_sum;
        // double ps_f1_sum = 0.0, pi_f1_sum = 0.0;
        // double ps_f2_sum = 0.0, pi_f2_sum = 0.0;

        double ps_f1_sum = 0.0;
        double pi_f1_sum = 0.0;

        /* F1: old scale */
        for (size_t quadratureIndex = 0; quadratureIndex < oldThetaValues.size(); quadratureIndex++) {
          double theta = oldThetaValues(quadratureIndex);
          double th_weight = oldThetaWeights(quadratureIndex);

          double func1 = 0.0;
          double ps_sum = 0.0;
          double pi_sum = 0.0;

          for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
            EquatingRecipes::Structures::CommonItemSpecification commonItem = commonItems[itemIndex];

            for (size_t categoryIndex = 1; categoryIndex <= commonItem.numberOfCategories; categoryIndex++) {
              double uJK = commonItem.scoringFunctionValues(categoryIndex);

              double p_origi = irtModelFunctions.probOld(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         true,
                                                         x[0],
                                                         x[1]);
              double p_trans = irtModelFunctions.probOld(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         false,
                                                         x[0],
                                                         x[1]);

              double ps = irtModelFunctions.itemResponseFunctionDerivativeOldOverS(commonItem,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   x[0],
                                                                                   x[1]);

              double pi = irtModelFunctions.itemResponseFunctionDerivativeOldOverI(commonItem,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   x[0],
                                                                                   x[1]);

              func1 += uJK * (p_origi - p_trans);
              ps_sum += uJK * ps;
              pi_sum += uJK * pi;
            }
          }

          ps_f1_sum += func1 * ps_sum * th_weight;
          pi_f1_sum += func1 * pi_sum * th_weight;
        }

        double ps_f2_sum = 0.0;
        double pi_f2_sum = 0.0;

        /* F2: new scale */
        for (size_t quadratureIndex = 0; quadratureIndex < newThetaValues.size(); quadratureIndex++) {
          double theta = newThetaValues(quadratureIndex);
          double th_weight = newThetaWeights(quadratureIndex);

          double func2 = 0.0;
          double ps_sum = 0.0;
          double pi_sum = 0.0;

          for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
            EquatingRecipes::Structures::CommonItemSpecification commonItem = commonItems[itemIndex];

            for (size_t categoryIndex = 1; categoryIndex <= commonItem.numberOfCategories; categoryIndex++) {
              double uJK = commonItem.scoringFunctionValues(categoryIndex);

              double p_origi = irtModelFunctions.probNew(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         true,
                                                         x[0],
                                                         x[1]);

              double p_trans = irtModelFunctions.probNew(commonItem,
                                                         categoryIndex,
                                                         theta,
                                                         false,
                                                         x[0],
                                                         x[1]);
              double ps = irtModelFunctions.itemResponseFunctionDerivativeNewOverS(commonItem,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   x[0],
                                                                                   x[1]);

              double pi = irtModelFunctions.itemResponseFunctionDerivativeNewOverI(commonItem,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   x[0],
                                                                                   x[1]);
              func2 += uJK * (p_origi - p_trans);
              ps_sum += uJK * ps;
              pi_sum += uJK * pi;
            }
          }

          ps_f2_sum += func2 * ps_sum * th_weight;
          pi_f2_sum += func2 * pi_sum * th_weight;
        }

        double sym_f1 = 0.0;
        double sym_f2 = 0.0;

        /* symmetric or non-symmetric setting*/
        if (this->symmetry == EquatingRecipes::Structures::Symmetry::SYMMETRIC) {
          sym_f1 = 1.0;
          sym_f2 = 1.0;
        } else {
          if (this->symmetry == EquatingRecipes::Structures::Symmetry::OLD_SCALE) {
            sym_f1 = 1.0;
            sym_f2 = 0.0;
          } else {
            sym_f1 = 0.0;
            sym_f2 = 1.0;
          }
        }

        /* function standardization */
        double w1_sum = oldThetaWeights.sum();
        double w2_sum = newThetaWeights.sum();

        if (this->functionStandardization) {
          grad[0] = -2.0 * (sym_f1 * ps_f1_sum / w1_sum + sym_f2 * ps_f2_sum / w2_sum);
          grad[1] = -2.0 * (sym_f1 * pi_f1_sum / w1_sum + sym_f2 * pi_f2_sum / w2_sum);
        } else {
          grad[0] = -2.0 * (sym_f1 * ps_f1_sum + sym_f2 * ps_f2_sum);
          grad[1] = -2.0 * (sym_f1 * pi_f1_sum + sym_f2 * pi_f2_sum);
        }
      }
    };
  } // namespace Implementation
} // namespace EquatingRecipes

#endif