#ifndef EQUATING_RECIPES_HAEBERA_FUNCION_HPP
#define EQUATING_RECIPES_HAEBERA_FUNCION_HPP

#include <vector>
#include <Eigen/Core>

#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/optimization_function.hpp>

namespace EquatingRecipes {
  class HaebaraFunction : public EquatingRecipes::OptimizationFunction {
  protected:
    /*------------------------------------------------------------------------------
    Functionality:
      calculate the value of the Haebara criterion function.
      
      Input:
    x: the point of slope and intercept
      Output:    
    the value of the criterion function

    Author: Seonghoon Kim
    Date of last revision 9/25/08
  ------------------------------------------------------------------------------*/
    double functionValue(const std::vector<double>& x) override {
      double func1_sum = 0.0;
      double func2_sum = 0.0;
      double cat_sum = 0.9;

      /* Q1: old scale */
      for (size_t quadratureIndex = 0; quadratureIndex < oldThetaValues.size(); quadratureIndex++) {
        double theta = oldThetaValues(quadratureIndex);
        double th_weight = oldThetaWeights(quadratureIndex);

        double func1 = 0.0;
        for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
          for (size_t categoryIndex = 1; categoryIndex <= commonItems[itemIndex].numberOfCategories; categoryIndex++) {
            double p_origi = irtModelFunctions.probOld(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       true,
                                                       x[0],
                                                       x[1]);
            double p_trans = irtModelFunctions.probOld(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       false,
                                                       x[0],
                                                       x[1]);
            func1 += std::pow(p_origi - p_trans, 2);
          }
        }
        func1_sum += func1 * th_weight;
      }

      /* Q2: new scale */
      for (size_t quadratureIndex = 0; quadratureIndex < newThetaValues.size(); quadratureIndex++) {
        double theta = newThetaValues(quadratureIndex);
        double th_weight = newThetaWeights(quadratureIndex);

        double func2 = 0.0;
        for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
          for (size_t categoryIndex = 1; categoryIndex <= commonItems[itemIndex].numberOfCategories; categoryIndex++) {
            double p_origi = irtModelFunctions.probNew(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       true,
                                                       x[0],
                                                       x[1]);

            double p_trans = irtModelFunctions.probNew(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       false,
                                                       x[0],
                                                       x[1]);

            func2 += std::pow(p_origi - p_trans, 2);
          }
          if (quadratureIndex == 0) {
            cat_sum += commonItems[itemIndex].numberOfCategories;
          }
        }

        func2_sum += func2 * th_weight;
      }

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

      double w1_sum = oldThetaWeights.sum();
      double w2_sum = newThetaWeights.sum();

      /* function standardization */
      if (this->functionStandardization) {
        return (sym_f1 * func1_sum / (cat_sum * w1_sum) + sym_f2 * func2_sum / (cat_sum * w2_sum));
      } else {
        return (sym_f1 * func1_sum + sym_f2 * func2_sum);
      }
    }

    /*------------------------------------------------------------------------------
    Functionality:
      calculate the partial derivatives of the Haebara criterion function with
      respect to slope (S or A) and intercept (I or B).
      
      Input:
    x: the point of slope and intercept
      Output:    
    the gradient of the criterion function w.r.t. A and B.

    Author: Seonghoon Kim
    Date of last revision 9/25/08
  ------------------------------------------------------------------------------*/
    void functionGradient(const std::vector<double>& x,
                          std::vector<double>& grad) override {
      double ps_f1_sum = 0.0;
      double pi_f1_sum = 0.0;
      double cat_sum = 0.0;

      /* Q1: old scale */
      for (size_t quadratureIndex = 0; quadratureIndex < oldThetaValues.size(); quadratureIndex++) {
        double theta = oldThetaValues(quadratureIndex);
        double th_weight = oldThetaWeights(quadratureIndex);

        double ps_f1 = 0.0;
        double pi_f1 = 0.0;

        for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
          size_t numberOfCategories = commonItems[itemIndex].numberOfCategories;

          for (size_t categoryIndex = 1; categoryIndex <= numberOfCategories; categoryIndex++) {
            double p_origi = irtModelFunctions.probOld(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       true,
                                                       x[0],
                                                       x[1]);

            double p_trans = irtModelFunctions.probOld(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       false,
                                                       x[0],
                                                       x[1]);

            double ps = irtModelFunctions.itemResponseFunctionDerivativeOldOverS(commonItems[itemIndex],
                                                                                 categoryIndex,
                                                                                 theta,
                                                                                 x[0],
                                                                                 x[1]);

            double pi = irtModelFunctions.itemResponseFunctionDerivativeOldOverI(commonItems[itemIndex],
                                                                                 categoryIndex,
                                                                                 theta,
                                                                                 x[0],
                                                                                 x[1]);

            ps_f1 += (p_origi - p_trans) * ps;
            pi_f1 += (p_origi - p_trans) * pi;
          }

          if (quadratureIndex == 0) {
            cat_sum += commonItems[itemIndex].numberOfCategories;
          }
        }
        
        ps_f1_sum += ps_f1 * th_weight;
        pi_f1_sum += pi_f1 * th_weight;
      }

      /* Q2: new scale */

      double ps_f2_sum = 0.0;
      double pi_f2_sum = 0.0;

      for (size_t quadratureIndex = 0; quadratureIndex < newThetaValues.size(); quadratureIndex++) {
        double theta = newThetaValues(quadratureIndex);
        double th_weight = newThetaWeights(quadratureIndex);

        double ps_f2 = 0.0;
        double pi_f2 = 0.0;

        for (size_t itemIndex = 0; itemIndex < commonItems.size(); itemIndex++) {
          size_t numberOfCategories = commonItems[itemIndex].numberOfCategories;

          for (size_t categoryIndex = 1; categoryIndex <= numberOfCategories; categoryIndex++) {
            double p_origi = irtModelFunctions.probNew(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       true,
                                                       x[0],
                                                       x[1]);
            double p_trans = irtModelFunctions.probNew(commonItems[itemIndex],
                                                       categoryIndex,
                                                       theta,
                                                       false,
                                                       x[0],
                                                       x[1]);

            double ps = irtModelFunctions.itemResponseFunctionDerivativeNewOverS(commonItems[itemIndex],
                                                                                 categoryIndex,
                                                                                 theta,
                                                                                 x[0],
                                                                                 x[1]);

            double pi = irtModelFunctions.itemResponseFunctionDerivativeNewOverI(commonItems[itemIndex],
                                                                                 categoryIndex,
                                                                                 theta,
                                                                                 x[0],
                                                                                 x[1]);

            ps_f2 += (p_origi - p_trans) * ps;
            pi_f2 += (p_origi - p_trans) * pi;
          }
        }
        ps_f2_sum += ps_f2 * th_weight;
        pi_f2_sum += pi_f2 * th_weight;
      }

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
      if (this->functionStandardization) {
        double w1_sum = oldThetaWeights.sum();
        double w2_sum = newThetaWeights.sum();

        grad[0] = -2.0 * (sym_f1 * ps_f1_sum / (cat_sum * w1_sum) + sym_f2 * ps_f2_sum / (cat_sum * w2_sum));
        grad[1] = -2.0 * (sym_f1 * pi_f1_sum / (cat_sum * w1_sum) + sym_f2 * pi_f2_sum / (cat_sum * w2_sum));
      } else {
        grad[0] = -2.0 * (sym_f1 * ps_f1_sum + sym_f2 * ps_f2_sum);
        grad[1] = -2.0 * (sym_f1 * pi_f1_sum + sym_f2 * pi_f2_sum);
      }
    }
  };
} // namespace EquatingRecipes

#endif