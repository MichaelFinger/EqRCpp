#ifndef IRT_MODEL_FUNCTIONS_HPP
#define IRT_MODEL_FUNCTIONS_HPP

#include <cmath>

#include <Eigen/Dense>

#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/form_type.hpp>

namespace EquatingRecipes {
  namespace Implementation {
    struct IRTModelFunctions {
      /* functions associated with different IRT models */

      /*------------------------------------------------------------------------------
    Functionality:
      Calculate the value of the category characteristic curve, P_jk, by model.

    Input:
      Item: One item of struct ItemSpec type 
      CatID: Response category ID in question (1 through CatNum)
      theta: ability value

    Output:
      Return the probability of an examinee having the value of theta
      responding in category ID
    Note:
      Uses arbitrarily the "old" scale with S = 1.0 and I = 0.0, since
      the functions in ProbDeriv.c are used.

    Author: Seonghoon Kim
    Date of last revision 9/25/08
  ------------------------------------------------------------------------------*/
      double itemResponseFunction(const EquatingRecipes::Structures::ItemSpecification& item,
                                  const size_t& categoryIndex,
                                  const double& theta) {
        double prob;

        switch (item.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            prob = itemResponseFunction3PL(categoryIndex,
                                           theta,
                                           item.scalingConstant,
                                           item.a(2),
                                           item.b(2),
                                           item.c(2),
                                           EquatingRecipes::Structures::FormType::OLD,
                                           1,
                                           0);

            break;
          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            prob = itemResponseFunctionLGR(item.numberOfCategories,
                                           categoryIndex,
                                           theta,
                                           item.scalingConstant,
                                           item.a(2),
                                           item.b,
                                           EquatingRecipes::Structures::FormType::OLD,
                                           1,
                                           0);
            break;
          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            prob = itemResponseFunctionGPC(item.numberOfCategories,
                                           categoryIndex,
                                           theta,
                                           item.scalingConstant,
                                           item.a(2),
                                           item.b,
                                           EquatingRecipes::Structures::FormType::OLD,
                                           1,
                                           0);
            break;
          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            prob = itemResponseFunctionNRM(item.numberOfCategories,
                                           categoryIndex,
                                           theta,
                                           item.a,
                                           item.c,
                                           EquatingRecipes::Structures::FormType::OLD,
                                           1,
                                           0);
            break;
          default:
            break;
        }

        return prob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        By model, calculate the first (partial) derivative of P_jk
        with respect to ability theta.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative(const EquatingRecipes::Structures::ItemSpecification& item,
                                            const size_t& categoryIndex,
                                            const double& theta) {
        double derivative;

        switch (item.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            derivative = itemResponseFunctionDerivative3PL(categoryIndex,
                                                           theta,
                                                           item.scalingConstant,
                                                           item.a(2),
                                                           item.b(2),
                                                           item.c(2));
            break;
          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            derivative = itemResponseFunctionDerivativeLGR(item.numberOfCategories,
                                                           categoryIndex,
                                                           theta,
                                                           item.scalingConstant,
                                                           item.a(2),
                                                           item.b);
            break;
          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            derivative = itemResponseFunctionDerivativeGPC(item.numberOfCategories,
                                                           categoryIndex,
                                                           theta,
                                                           item.scalingConstant,
                                                           item.a(2),
                                                           item.b);

            break;
          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            derivative = itemResponseFunctionDerivativeNRM(item.numberOfCategories,
                                                           categoryIndex,
                                                           theta,
                                                           item.a,
                                                           item.c);
            break;
          default:
            break;
        }

        return derivative;

        return 0;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Calculate the value of the characteristic curve on the old scale by model.

      For details about partial derivatives, refer to the following report:

      Kim, S., & Kolen, M.J. (2005). Methods for obtaining a common scale under
          unidimensional IRT models: A technical review and further extensions
          (Iowa Testing Programs Occasional Paper, No. 52). The University of
          Iowa.

      Input:
      Item: One item of struct CommonItemSpec type 
      CatID: Response category ID in question (1 through CatNum)
      theta: ability value
      original: on or off
        if on, then old scale's item parameters are used with
        S = 1.0 and I = 0.0. In this case, S and I are over-argumented.
        Otherwise, transformed item parameters (from new scale) are
        used through S and I.
      S: slope of the linear tranformation
      I: intercept of the linear transformation

      Output:
      Return either original (with old parameters and S = 1.0 and I = 0.0)
      or transformed (with new parameters and S and I) probability
      on the old scale

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double probOld(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                     const size_t& categoryIndex,
                     const double& theta,
                     const bool& original,
                     const double& S,
                     const double& I) {
        double prob;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            if (original) {
              prob = itemResponseFunction3PL(categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB(2),
                                             commonItem.oldC(2),
                                             EquatingRecipes::Structures::FormType::OLD,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunction3PL(categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB(2),
                                             commonItem.newC(2),
                                             EquatingRecipes::Structures::FormType::OLD,
                                             S,
                                             I);
            }
            break;
          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            if (original) {
              prob = itemResponseFunctionLGR(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionLGR(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             S,
                                             I);
            }
            break;
          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            if (original) {
              prob = itemResponseFunctionGPC(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionGPC(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             S,
                                             I);
            }
            break;
          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            if (original) {
              prob = itemResponseFunctionNRM(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.oldA,
                                             commonItem.oldC,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionNRM(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.newA,
                                             commonItem.newC,
                                             EquatingRecipes::Structures::FormType::OLD,
                                             S,
                                             I);
            }
            break;
          default:
            break;
        }
        return prob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Calculate the value of the characteristic curve on the new scale by model.

        Input:
      Item: One item of struct CommonItemSpec type 
      CatID: Response category ID in question (1 through CatNum)
      theta: ability value
      original: on or off
        if on, then new scale's item parameters are used with
        S = 1.0 and I = 0.0. In this case, S and I are over-argumented.
        Otherwise, transformed item parameters (from old scale) are
        used through S and I.
      S: slope of the linear tranformation
      I: intercept of the linear transformation

      Output:
      Return either original (with new parameters and S = 1.0 and I = 0.0)
      or transformed (with old parameters and S and I) probability
      on the new scale

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double probNew(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                     const size_t& categoryIndex,
                     const double& theta,
                     const bool& original,
                     const double& S,
                     const double& I) {
        double prob;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            if (original) {
              prob = itemResponseFunction3PL(categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB(2),
                                             commonItem.newC(2),
                                             EquatingRecipes::Structures::FormType::NEW,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunction3PL(categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB(2),
                                             commonItem.oldC(2),
                                             EquatingRecipes::Structures::FormType::NEW,
                                             S,
                                             I);
            }
            break;

          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            if (original) {
              prob = itemResponseFunctionLGR(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionLGR(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             S,
                                             I);
            }
            break;

          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            if (original) {
              prob = itemResponseFunctionGPC(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.newA(2),
                                             commonItem.newB,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionGPC(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.scalingConstant,
                                             commonItem.oldA(2),
                                             commonItem.oldB,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             S,
                                             I);
            }
            break;

          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            if (original) {
              prob = itemResponseFunctionNRM(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.newA,
                                             commonItem.newC,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             1,
                                             0);
            } else {
              prob = itemResponseFunctionNRM(commonItem.numberOfCategories,
                                             categoryIndex,
                                             theta,
                                             commonItem.oldA,
                                             commonItem.oldC,
                                             EquatingRecipes::Structures::FormType::NEW,
                                             S,
                                             I);
            }
            break;

          default:
            break;
        }
        return prob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        By model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to S.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeOldOverS(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                                                    const size_t& categoryIndex,
                                                    const double& theta,
                                                    const double& S,
                                                    const double& I) {
        double derivative;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            derivative = itemResponseFunctionDerivative3PLOldOverS(categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB(2),
                                                                   commonItem.newC(2),
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            derivative = itemResponseFunctionDerivativeLGROldOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            derivative = itemResponseFunctionDerivativeGPCOldOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            derivative = itemResponseFunctionDerivativeNRMOldOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.newA,
                                                                   commonItem.newC,
                                                                   S,
                                                                   I);

            break;
          default:
            break;
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        By model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to I.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeOldOverI(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                                                    const size_t& categoryIndex,
                                                    const double& theta,
                                                    const double& S,
                                                    const double& I) {
        double derivative;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            derivative = itemResponseFunctionDerivative3PLOldOverI(categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB(2),
                                                                   commonItem.newC(2),
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            derivative = itemResponseFunctionDerivativeLGROldOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            derivative = itemResponseFunctionDerivativeGPCOldOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.newA(2),
                                                                   commonItem.newB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            derivative = itemResponseFunctionDerivativeNRMOldOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.newA,
                                                                   commonItem.newC,
                                                                   S,
                                                                   I);

            break;
          default:
            break;
        }

        return derivative;
      }

      double itemResponseFunctionDerivativeNewOverS(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                                                    const size_t& categoryIndex,
                                                    const double& theta,
                                                    const double& S,
                                                    const double& I) {
        double derivative;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            derivative = itemResponseFunctionDerivative3PLNewOverS(categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB(2),
                                                                   commonItem.oldC(2),
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            derivative = itemResponseFunctionDerivativeLGRNewOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            derivative = itemResponseFunctionDerivativeGPCNewOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            derivative = itemResponseFunctionDerivativeNRMNewOverS(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.oldA,
                                                                   commonItem.oldC,
                                                                   S,
                                                                   I);

            break;
          default:
            break;
        }

        return derivative;
      }

      double itemResponseFunctionDerivativeNewOverI(const EquatingRecipes::Structures::CommonItemSpecification& commonItem,
                                                    const size_t& categoryIndex,
                                                    const double& theta,
                                                    const double& S,
                                                    const double& I) {
        double derivative;

        switch (commonItem.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            derivative = itemResponseFunctionDerivative3PLNewOverI(categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB(2),
                                                                   commonItem.oldC(2),
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            derivative = itemResponseFunctionDerivativeLGRNewOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            derivative = itemResponseFunctionDerivativeGPCNewOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.scalingConstant,
                                                                   commonItem.oldA(2),
                                                                   commonItem.oldB,
                                                                   S,
                                                                   I);

            break;

          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            derivative = itemResponseFunctionDerivativeNRMNewOverI(commonItem.numberOfCategories,
                                                                   categoryIndex,
                                                                   theta,
                                                                   commonItem.oldA,
                                                                   commonItem.oldC,
                                                                   S,
                                                                   I);

            break;
          default:
            break;
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the 3PL model, calculate the first (partial) derivative of P_j
        with respect to theta.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative3PL(const size_t& categoryIndex,
                                               const double& theta,
                                               const double& D,
                                               const double& a,
                                               const double& b,
                                               const double& c) {
        double derivative;
        double uj;

        uj = std::exp(D * a * (theta - b));

        if (categoryIndex == 1) {
          derivative = -1.0 * D * a * (1.0 - c) * uj / ((1.0 + uj) * (1.0 + uj));
        } else if (categoryIndex == 2) {
          derivative = D * a * (1.0 - c) * uj / ((1.0 + uj) * (1.0 + uj));
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GR model, calculate the first (partial) derivative of P_jk
        with respect to theta.
        
      Note:
        b[2] for the first item-step parameter
        b[CatNum] for the last item-step parameter

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeLGR(const size_t& numberOfCategories,
                                               const size_t& categoryIndex,
                                               const double& theta,
                                               const double& D,
                                               const double& a,
                                               const Eigen::VectorXd& b) {
        double cp_jk;
        double cp_jk1;

        if (categoryIndex == 1) {
          cp_jk = 1.0;
          cp_jk1 = 1.0 / (1.0 + std::exp(-1.0 * D * a * (theta - b(categoryIndex + 1))));
        } else {
          if (categoryIndex < numberOfCategories) {
            cp_jk = 1.0 / (1.0 + std::exp(-1.0 * D * a * (theta - b(categoryIndex))));
            cp_jk1 = 1.0 / (1.0 + std::exp(-1.0 * D * a * (theta - b(categoryIndex + 1))));
          } else {
            cp_jk = 1.0 / (1.0 + std::exp(-1.0 * D * a * (theta - b(categoryIndex))));
            cp_jk1 = 0.0;
          }
        }
        return (D * a * (cp_jk * (1.0 - cp_jk) - cp_jk1 * (1.0 - cp_jk1)));
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GPC model, calculate the first (partial) derivative of P_jk
        with respect to ability theta.
        
      Note:
        b[2] for the actual first item-step parameter
        b[CatNum] for the last item-step parameter
        is dealt with as a special case of the nominal response model

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeGPC(const size_t& numberOfCategories,
                                               const size_t& categoryIndex,
                                               const double& theta,
                                               const double& D,
                                               const double& a,
                                               const Eigen::VectorXd& b) {
        double aUSum = 0.0;
        double vj = 0.0;

        for (size_t index = 1; index <= numberOfCategories; index++) {
          double aSum = static_cast<double>(index + 1) * D * a;
          double bSum = -1.0 * D * a * b(Eigen::seq(2, index)).sum(); /* b[0] = 0 */
          double ujk = std::exp(aSum * theta + bSum);
          vj += ujk;
          aUSum += aSum * ujk;
        }

        double aSum = static_cast<double>(categoryIndex) * D * a;
        double bSum = -1.0 * D * a * b(Eigen::seq(2, categoryIndex)).sum(); /* b[0] = 0 */
        bSum *= -1.0 * D * a;
        double ujk = std::exp(aSum * theta + bSum);

        double derivative = aSum * ujk / vj - ujk * aSum / std::pow(vj, 2);

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the NR model, calculate the first (partial) derivative of P_jk
        with respect to ability theta.
        
      Note:
        a[1..CatNum] for the discrimination parameters
        c[1..CatNum] for the intercept parameters
        No scaling constant

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeNRM(const size_t& numberOfCategories,
                                               const size_t& categoryIndex,
                                               const double& theta,
                                               const Eigen::VectorXd& a,
                                               const Eigen::VectorXd& c) {
        double vj = 0.0;
        double aUSum = 0.0;

        for (size_t index = 1; index <= numberOfCategories; index++) {
          double ujk = std::exp(a(index) * theta + c(index));
          vj += ujk;
          aUSum += a(index) * ujk;
        }

        double ujk = std::exp(a(categoryIndex) * theta + c(categoryIndex));

        double derivative = a(categoryIndex) * ujk / vj - ujk * aUSum / std::pow(vj, 2);

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Calculate the value of the characteristic curve for the 3PL model.
        Input: 
      CatID: response category, 1 or 2 (1 for incorrect; 2 for correct)
      theta: ability value
      D: scaling constant (typically 1.7)
      a, b, c: item parameters
      scale: "old" or "new" ability scale, on which item and ability parameter
        estimates are placed on.
      S: slope of the linear tranformation
      I: intercept of the linear transformation

      Output:
      Return either original or transformed probability
      For the original probability, S = 1 and I = 0 with parameters on the
      reference scale.
      For the transformed probability, S and I with parameters on the
      transformed scale.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunction3PL(const size_t& categoryIndex,
                                     const double& theta,
                                     const double& D,
                                     const double& a,
                                     const double& b,
                                     const double& c,
                                     const EquatingRecipes::Structures::FormType& scale,
                                     const double& S,
                                     const double& I) {
        double aPar;
        double bPar;
        double cPar;

        if (scale == EquatingRecipes::Structures::FormType::OLD) {
          /* new-to-old scale transformation */
          aPar = a / S;
          bPar = S * b + I;
          cPar = c;
        } else {
          /* old-to-new scale transformation */
          aPar = S * a;
          bPar = (b - I) / S;
          cPar = c;
        }

        double deviate = D * aPar * (theta - bPar);

        double prob = cPar + (1.0 - cPar) / (1.0 + std::exp(-1.0 * deviate));

        if (categoryIndex == 1) {
          prob = 1.0 - prob;
        }

        return prob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the 3PL model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to S.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative3PLOldOverS(const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const double& nb,
                                                       const double& nc,
                                                       const double& S,
                                                       const double& I) {
        double as, cs;
        double ps, qs; /* Probability with transformed item parameters */
        double ps_over_S;

        /* Scale Transformation */
        as = na / S;
        cs = nc;

        ps = itemResponseFunction3PL(2,
                                     theta,
                                     D,
                                     na,
                                     nb,
                                     nc,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);
        qs = 1.0 - ps;
        ps_over_S = -D * as * ((theta - I) / S) * (ps - cs) * qs / (1.0 - cs);

        double derivative = ps_over_S;

        if (categoryIndex == 1) {
          derivative *= -1.0;
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the 3PL model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to I.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative3PLOldOverI(const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const double& nb,
                                                       const double& nc,
                                                       const double& S,
                                                       const double& I) {
        double as, cs;
        double ps, qs; /* Probability with transformed item parameters */
        double ps_over_I;

        /* Scale Transformation */
        as = na / S;
        cs = nc;

        ps = itemResponseFunction3PL(2,
                                     theta,
                                     D,
                                     na,
                                     nb,
                                     nc,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);
        qs = 1.0 - ps;
        ps_over_I = -D * as * (ps - cs) * qs / (1.0 - cs);

        if (categoryIndex == 1) {
          ps_over_I *= -1.0;
        }

        return ps_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the 3PL model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to S.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative3PLNewOverS(const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const double& ob,
                                                       const double& oc,
                                                       const double& S,
                                                       const double& I) {
        double cr;
        double pr, qr; /* Probability with transformed item parameter estimates */
        double pr_over_S;

        /* Scale Transformation */
        cr = oc;

        pr = itemResponseFunction3PL(2,
                                     theta,
                                     D,
                                     oa,
                                     ob,
                                     oc,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);
        qr = 1.0 - pr;

        pr_over_S = D * oa * theta * (pr - cr) * qr / (1.0 - cr);

        double derivative = pr_over_S;

        if (categoryIndex == 1) {
          derivative *= -1.0;
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the 3PL model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to I.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivative3PLNewOverI(const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const double& ob,
                                                       const double& oc,
                                                       const double& S,
                                                       const double& I) {
        double cr;
        double pr, qr; /* Probability with transformed item parameter estimates */
        double pr_over_I;

        /* Scale Transformation */
        cr = oc;

        pr = itemResponseFunction3PL(2,
                                     theta,
                                     D,
                                     oa,
                                     ob,
                                     oc,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);

        qr = 1.0 - pr;

        pr_over_I = D * oa * (pr - cr) * qr / (1.0 - cr);

        double derivative = pr_over_I;

        if (categoryIndex == 1) {
          derivative *= -1.0;
        }

        return derivative;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the logistic graded respnse model, calculate the cumulative category
        characteristic curve. Use the notation used in Kim and Kolen (2005).

        Input
      scale: "old" or "new"
        Output
      Return either original or transformed probability.
      For the original probability, S = 1 and I = 0 with parameters on the
      reference scale.
      For the transfomred probability, S and I with parameters on the
      transformed scale.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double cumulativeResponseProbabilityLGR(const size_t& categoryIndex,
                                              const double& theta,
                                              const double& D,
                                              const double& a,
                                              const double& b,
                                              const EquatingRecipes::Structures::FormType& scale,
                                              const double& S,
                                              const double& I) {
        double as, bs, ar, br;

        double aTransform;
        double bTransform;
        double cprob;

        if (categoryIndex == 0)
          cprob = 1.0;
        else {
          if (scale == EquatingRecipes::Structures::FormType::OLD) {
            /* new-to-old transformation */
            aTransform = a / S;
            bTransform = S * b + I;

          } else {
            /* old-to-new transformation */
            aTransform = S * a;
            bTransform = (b - I) / S;
          }

          cprob = 1.0 / (1.0 + exp(-1.0 * D * aTransform * (theta - bTransform)));
        }

        return cprob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the logistic graded response model, calculate the category
        characteristic curve. Use the notation used in Kim and Kolen (2005).

        Input
      CatNum: number of categories
      CatID: category response ID
      theta: ability value
      D: scaling constant
      a: discrimination parameter
      b[]: difficulty parameter array
          b[2] for the first difficulty parameter
          b[CatNum] for the last category
      scale: "old" or "new"
      S: slope of linear transformation
      I: intercept of linear transformation
        Output
      Return either original or transformed probability.
      For the original probability, S = 1 and I = 0 with parameters on the
      reference scale.
      For the transfomred probability, S and I with parameters on the
      transformed scale.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionLGR(const size_t& numberOfCategories,
                                     const size_t categoryIndex,
                                     const double& theta,
                                     const double& D,
                                     const double& a,
                                     const Eigen::VectorXd& b,
                                     const EquatingRecipes::Structures::FormType& scale,
                                     const double& S,
                                     const double& I) {
        double pre_cp;
        double pos_cp;

        if (scale == EquatingRecipes::Structures::FormType::OLD) {
          pre_cp = cumulativeResponseProbabilityLGR(categoryIndex,
                                                    theta,
                                                    D,
                                                    a,
                                                    b(categoryIndex),
                                                    scale,
                                                    S,
                                                    I);

          if (categoryIndex < numberOfCategories)
            pos_cp = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                      theta,
                                                      D,
                                                      a,
                                                      b(categoryIndex + 1),
                                                      scale,
                                                      S,
                                                      I);
          else {
            pos_cp = 0.0;
          }
        } else {
          pre_cp = cumulativeResponseProbabilityLGR(categoryIndex,
                                                    theta,
                                                    D,
                                                    a, b(categoryIndex),
                                                    EquatingRecipes::Structures::FormType::NEW,
                                                    S,
                                                    I);
          if (categoryIndex < numberOfCategories) {
            pos_cp = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                      theta,
                                                      D,
                                                      a,
                                                      b(categoryIndex + 1),
                                                      EquatingRecipes::Structures::FormType::NEW,
                                                      S,
                                                      I);
          } else {
            pos_cp = 0.0;
          }
        }

        double prob = pre_cp - pos_cp;

        return prob;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GR model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to S.
        
        Note:
      nb[2] for the first item-step parameter
      nb[CatNum] for the last item-step parameter

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeLGROldOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const Eigen::VectorXd& nb,
                                                       const double& S,
                                                       const double& I) {
        double as;
        double pre_cps;
        double pos_cps;
        double pre_cps_over_S;
        double pos_cps_over_S;
        double ps_over_S;

        as = na / S;

        if (categoryIndex == 1) {
          pos_cps = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                     theta,
                                                     D,
                                                     na,
                                                     nb(categoryIndex + 1),
                                                     EquatingRecipes::Structures::FormType::OLD,
                                                     S,
                                                     I);

          pos_cps_over_S = -1.0 * D * as * ((theta - I) / S) * pos_cps * (1.0 - pos_cps);
          ps_over_S = 0.0 - pos_cps_over_S;
        } else {
          if (categoryIndex < numberOfCategories) {
            pre_cps = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pos_cps = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex + 1),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pre_cps_over_S = -1.0 * D * as * ((theta - I) / S) * pre_cps * (1.0 - pre_cps);
            pos_cps_over_S = -1.0 * D * as * ((theta - I) / S) * pos_cps * (1.0 - pos_cps);
            ps_over_S = pre_cps_over_S - pos_cps_over_S;
          } else { /* CatId == CatNum */
            pre_cps = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pre_cps_over_S = -1.0 * D * as * ((theta - I) / S) * pre_cps * (1.0 - pre_cps);
            ps_over_S = pre_cps_over_S - 0.0;
          }
        }

        return ps_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GR model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to I.
        
        Note:
      nb[2] for the first item-step parameter
      nb[CatNum] for the last item-step parameter

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeLGROldOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const Eigen::VectorXd& nb,
                                                       const double& S,
                                                       const double& I) {
        double as;
        double pre_cps, pos_cps;
        double pre_cps_over_I;
        double pos_cps_over_I;
        double ps_over_I;

        as = na / S;

        if (categoryIndex == 1) {
          pos_cps = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                     theta,
                                                     D,
                                                     na,
                                                     nb(categoryIndex + 1),
                                                     EquatingRecipes::Structures::FormType::OLD,
                                                     S,
                                                     I);
          pos_cps_over_I = -1.0 * D * as * pos_cps * (1.0 - pos_cps);
          ps_over_I = 0.0 - pos_cps_over_I;
        } else {
          if (categoryIndex < numberOfCategories) {
            pre_cps = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pos_cps = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex + 1),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pre_cps_over_I = -1.0 * D * as * pre_cps * (1.0 - pre_cps);
            pos_cps_over_I = -1.0 * D * as * pos_cps * (1.0 - pos_cps);
            ps_over_I = pre_cps_over_I - pos_cps_over_I;
          } else {
            pre_cps = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       na,
                                                       nb(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::OLD,
                                                       S,
                                                       I);
            pre_cps_over_I = -1.0 * D * as * pre_cps * (1.0 - pre_cps);
            ps_over_I = pre_cps_over_I - 0.0;
          }
        }
        return ps_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GR model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to S.
        
        Note:
      ob[2] for the first item-step parameter
      ob[CatNum] for the last item-step parameter

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeLGRNewOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const Eigen::VectorXd& ob,
                                                       const double& S,
                                                       const double& I) {
        double pre_cpr;
        double pos_cpr;
        double pre_cpr_over_S;
        double pos_cpr_over_S;
        double pr_over_S;

        if (categoryIndex == 1) {
          pos_cpr = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                     theta,
                                                     D,
                                                     oa,
                                                     ob(categoryIndex + 1),
                                                     EquatingRecipes::Structures::FormType::NEW,
                                                     S,
                                                     I);
          pos_cpr_over_S = D * oa * theta * pos_cpr * (1.0 - pos_cpr);
          pr_over_S = 0.0 - pos_cpr_over_S;
        } else {
          if (categoryIndex < numberOfCategories) {
            pre_cpr = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);

            pos_cpr = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex + 1),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);
            pre_cpr_over_S = D * oa * theta * pre_cpr * (1.0 - pre_cpr);
            pos_cpr_over_S = D * oa * theta * pos_cpr * (1.0 - pos_cpr);
            pr_over_S = pre_cpr_over_S - pos_cpr_over_S;
          } else {
            pre_cpr = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);
            pre_cpr_over_S = D * oa * theta * pre_cpr * (1.0 - pre_cpr);
            pr_over_S = pre_cpr_over_S - 0.0;
          }
        }
        return pr_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GR model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to I.
        
        Note:
      ob[2] for the first item-step parameter
      ob[CatNum] for the last item-step parameter

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeLGRNewOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const Eigen::VectorXd& ob,
                                                       const double& S,
                                                       const double& I) {
        double pre_cpr, pos_cpr;
        double pre_cpr_over_I, pos_cpr_over_I;
        double pr_over_I;

        if (categoryIndex == 1) {
          pos_cpr = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                     theta,
                                                     D,
                                                     oa,
                                                     ob(categoryIndex + 1),
                                                     EquatingRecipes::Structures::FormType::NEW,
                                                     S,
                                                     I);
          pos_cpr_over_I = D * oa * pos_cpr * (1.0 - pos_cpr);
          pr_over_I = 0.0 - pos_cpr_over_I;
        } else {
          if (categoryIndex < numberOfCategories) {
            pre_cpr = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);
            pos_cpr = cumulativeResponseProbabilityLGR(categoryIndex + 1,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex + 1),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);
            pre_cpr_over_I = D * oa * pre_cpr * (1.0 - pre_cpr);
            pos_cpr_over_I = D * oa * pos_cpr * (1.0 - pos_cpr);
            pr_over_I = pre_cpr_over_I - pos_cpr_over_I;
          } else { /* CatID == CatNum */
            pre_cpr = cumulativeResponseProbabilityLGR(categoryIndex,
                                                       theta,
                                                       D,
                                                       oa,
                                                       ob(categoryIndex),
                                                       EquatingRecipes::Structures::FormType::NEW,
                                                       S,
                                                       I);
            pre_cpr_over_I = D * oa * pre_cpr * (1.0 - pre_cpr);
            pr_over_I = pre_cpr_over_I - 0.0;
          }
        }
        return pr_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the generalized partial credit model, calculate the category
        characteristic curve. Use the notation used in Kim and Kolen (2005).
        
        Note:
      is dealt with as a special case of the nominal response model.

        Input
      CatNum: number of categories
      CatID: category response ID
      theta: ability value
      D: scaling constant
      a: discrimination parameter
      b[]: difficulty parameter array
          It is assumed that b[1] = 0, but is not used.
          b[2] for the actual first difficulty (item-step) parameter
          b[CatNum] for the last category
      scale: "old" or "new"
      S: slope of linear transformation
      I: intercept of linear transformation
        Output
      Return either original or transformed probability.
      For the original probability, S = 1 and I = 0 with parameters on the
      reference scale.
      For the transfomred probability, S and I with parameters on the
      transformed scale.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionGPC(const size_t& numberOfCategories,
                                     const size_t categoryIndex,
                                     const double& theta,
                                     const double& D,
                                     const double& a,
                                     const Eigen::VectorXd& b,
                                     const EquatingRecipes::Structures::FormType& scale,
                                     const double& S,
                                     const double& I) {
        int k;
        int l;
        double vjs = 0.0;
        double vjr = 0.0;
        double a_sum;
        double b_sum;
        double as;
        double bs;
        double ar;
        double br;

        if (scale == EquatingRecipes::Structures::FormType::OLD) {
          for (k = 1; k <= numberOfCategories; k++) {
            a_sum = k * D * a;
            b_sum = 0.0;
            for (l = 2; l <= k; l++) {
              b_sum += b[l];
            }
            b_sum *= -1.0 * D * a;

            /* new-to-old transformation */
            as = a_sum / S;
            bs = b_sum - (I / S) * a_sum;
            vjs += std::exp(as * theta + bs);
          }
          a_sum = (categoryIndex)*D * a;
          b_sum = 0.0;
          for (l = 2; l <= categoryIndex; l++) {
            b_sum += b[l];
          }
          b_sum *= -D * a;

          /* new-to-old transformation */
          as = a_sum / S;
          bs = b_sum - (I / S) * a_sum;
          return std::exp(as * theta + bs) / vjs;
        } else {
          for (k = 1; k <= numberOfCategories; k++) {
            a_sum = k * D * a;
            b_sum = 0.0;
            for (l = 2; l <= k; l++) {
              b_sum += b[l];
            }
            b_sum *= -1.0 * D * a;

            /* old-to-new transformation */
            ar = S * a_sum;
            br = b_sum + I * a_sum;
            vjr += std::exp(ar * theta + br);
          }
          a_sum = categoryIndex * D * a;
          b_sum = 0.0;
          for (l = 2; l <= categoryIndex; l++) {
            b_sum += b[l];
          }
          b_sum *= -1.0 * D * a;

          /* old-to-new transformation */
          ar = S * a_sum;
          br = b_sum + I * a_sum;
          return std::exp(ar * theta + br) / vjr;
        }
      }

      /*------------------------------------------------------------------------------
    Functionality:
      Under the GPC model, calculate partial derivative of P* (from new-to-old
      transformation) with respect to S.
      
      Note:
      nb[2] for the actual first item-step parameter
      nb[CatNum] for the last item-step parameter
      is dealt with as a special case of the nominal response model

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeGPCOldOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const Eigen::VectorXd& nb,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double as, ps;
        double na_sum, as_ps_sum = 0.0;
        double ps_over_S;

        for (k = 1; k <= numberOfCategories; k++) {
          na_sum = k * D * na;
          as = na_sum / S; /* new-to-old transformation */
          ps = itemResponseFunctionGPC(numberOfCategories,
                                       k,
                                       theta,
                                       D,
                                       na,
                                       nb,
                                       EquatingRecipes::Structures::FormType::OLD,
                                       S,
                                       I);
          as_ps_sum += as * ps;
        }

        na_sum = categoryIndex * D * na; /* for the category in question */
        as = na_sum / S;                 /* new-to-old transformation */
        ps = itemResponseFunctionGPC(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     D,
                                     na,
                                     nb,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);
        ps_over_S = -ps * ((theta - I) / S) * (as - as_ps_sum);
        return ps_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GPC model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to I.
        
        Note:
      nb[2] for the actual first item-step parameter
      nb[CatNum] for the last item-step parameter
      is dealt with as a special case of the nominal response model

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeGPCOldOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& na,
                                                       const Eigen::VectorXd& nb,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double as;
        double ps;
        double na_sum = 0.0;
        double as_ps_sum = 0.0;
        double ps_over_I;

        for (k = 1; k <= numberOfCategories; k++) {
          na_sum = k * D * na;
          as = na_sum / S; /* new-to-old transformation */
          ps = itemResponseFunctionGPC(numberOfCategories,
                                       categoryIndex,
                                       theta,
                                       D,
                                       na,
                                       nb,
                                       EquatingRecipes::Structures::FormType::OLD,
                                       S,
                                       I);
          as_ps_sum += as * ps;
        }

        na_sum = categoryIndex * D * na; /* for the category in question */
        as = na_sum / S;                 /* new-to-old transformation */
        ps = itemResponseFunctionGPC(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     D,
                                     na,
                                     nb,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);

        ps_over_I = -1.0 * ps * (as - as_ps_sum);
        return ps_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GPC model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to S.
        
        Note:
      ob[2] for the actual first item-step parameter
      ob[CatNum] for the last item-step parameter
      is dealt with as a special case of the nominal response model

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeGPCNewOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const Eigen::VectorXd& ob,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double pr;
        double oa_sum, oa_pr_sum = 0.0;
        double pr_over_S;

        for (k = 1; k <= numberOfCategories; k++) {
          oa_sum = k * D * oa;
          pr = itemResponseFunctionGPC(numberOfCategories,
                                       k,
                                       theta,
                                       D,
                                       oa,
                                       ob,
                                       EquatingRecipes::Structures::FormType::NEW,
                                       S,
                                       I);
          oa_pr_sum += oa_sum * pr;
        }
        oa_sum = categoryIndex * D * oa; /* for the category in question */
        pr = itemResponseFunctionGPC(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     D,
                                     oa,
                                     ob,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);
        pr_over_S = pr * theta * (oa_sum - oa_pr_sum);
        return pr_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the GPC model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to I.
        
        Note:
      ob[2] for the actual first item-step parameter
      ob[CatNum] for the last item-step parameter
      is dealt with as a special case of the nominal response model

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeGPCNewOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const double& D,
                                                       const double& oa,
                                                       const Eigen::VectorXd& ob,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double pr;
        double oa_sum = 0.0;
        double oa_pr_sum = 0.0;
        double pr_over_I;

        for (k = 1; k <= numberOfCategories; k++) {
          oa_sum = k * D * oa;
          pr = itemResponseFunctionGPC(numberOfCategories,
                                       k,
                                       theta,
                                       D,
                                       oa,
                                       ob,
                                       EquatingRecipes::Structures::FormType::NEW,
                                       S,
                                       I);
          oa_pr_sum += oa_sum * pr;
        }
        oa_sum = categoryIndex * D * oa; /* for the category in question */
        pr = itemResponseFunctionGPC(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     D,
                                     oa,
                                     ob,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);
        pr_over_I = pr * (oa_sum - oa_pr_sum);
        return pr_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the nominal response model, calculate the category characteristic
        curve.
        
        Input
      CatNum: number of categories
      CatID: category response ID
      theta: ability value
      a[1..CatNum]: discrimination parameters
      c[1..CatNum]: intercept parameters
      scale: "old" or "new"
      S: slope of linear transformation
      I: intercept of linear transformation
      
      Note: No scaling constant
      
        Output
      Return either original or transformed probability.
      For the original probability, S = 1 and I = 0 with parameters on the
      reference scale.
      For the transfomred probability, S and I with parameters on the
      transformed scale.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionNRM(const size_t& numberOfCategories,
                                     const size_t& categoryIndex,
                                     const double& theta,
                                     const Eigen::VectorXd& a,
                                     const Eigen::VectorXd& c,
                                     const EquatingRecipes::Structures::FormType& scale,
                                     const double& S,
                                     const double& I) {
        int k;
        double vjs = 0.0;
        double vjr = 0.0;
        double as;
        double cs;
        double ar;
        double cr;

        if (scale == EquatingRecipes::Structures::FormType::OLD) {
          for (k = 1; k <= numberOfCategories; k++) {
            as = a(k) / S;
            cs = c(k) - (I / S) * a(k);
            vjs += exp(as * theta + cs);
          }
          as = a(categoryIndex) / S;
          cs = c(categoryIndex) - (I / S) * a(categoryIndex);
          return std::exp(as * theta + cs) / vjs;
        } else {
          for (k = 1; k <= numberOfCategories; k++) {
            ar = S * a(k);
            cr = c(k) + I * a(k);
            vjr += std::exp(ar * theta + cr);
          }
          ar = S * a(categoryIndex);
          cr = c(categoryIndex) + I * a(categoryIndex);
          return std::exp(ar * theta + cr) / vjr;
        }
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the NR model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to S.
        
        Note:
      na[1..CatNum] for the discrimination parameters
      nc[1..CatNum] for the intercept parameters
            No scaling constant

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeNRMOldOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const Eigen::VectorXd& na,
                                                       const Eigen::VectorXd& nc,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double as;
        double ps;
        double as_ps_sum = 0.0;
        double ps_over_S;

        for (k = 1; k <= numberOfCategories; k++) {
          as = na(k) / S;
          ps = itemResponseFunctionNRM(numberOfCategories,
                                       k,
                                       theta,
                                       na,
                                       nc,
                                       EquatingRecipes::Structures::FormType::OLD,
                                       S,
                                       I);
          as_ps_sum += as * ps;
        }
        as = na(categoryIndex) / S;
        ps = itemResponseFunctionNRM(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     na,
                                     nc,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);
        ps_over_S = -1.0 * ps * ((theta - I) / S) * (as - as_ps_sum);
        return ps_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the NR model, calculate partial derivative of P* (from new-to-old
        transformation) with respect to I.
        
        Note:
      na[1..CatNum] for the discrimination parameters
      nc[1..CatNum] for the intercept parameters
            No scaling constant

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeNRMOldOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const Eigen::VectorXd& na,
                                                       const Eigen::VectorXd& nc,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double as;
        double ps;
        double as_ps_sum = 0.0;
        double ps_over_I;

        for (k = 1; k <= numberOfCategories; k++) {
          as = na(k) / S;
          ps = itemResponseFunctionNRM(numberOfCategories,
                                       k,
                                       theta,
                                       na,
                                       nc,
                                       EquatingRecipes::Structures::FormType::OLD,
                                       S,
                                       I);
          as_ps_sum += as * ps;
        }

        as = na(categoryIndex) / S;
        ps = itemResponseFunctionNRM(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     na,
                                     nc,
                                     EquatingRecipes::Structures::FormType::OLD,
                                     S,
                                     I);
        ps_over_I = -1.0 * ps * (as - as_ps_sum);
        return ps_over_I;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the NR model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to S.
        
        Note:
      oa[1..CatNum] for the discrimination parameters
      oc[1..CatNum] for the intercept parameters
            No scaling constant

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeNRMNewOverS(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const Eigen::VectorXd& oa,
                                                       const Eigen::VectorXd& oc,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double pr;
        double oa_pr_sum = 0.0;
        double pr_over_S;

        for (k = 1; k <= numberOfCategories; k++) {
          pr = itemResponseFunctionNRM(numberOfCategories,
                                       k,
                                       theta,
                                       oa,
                                       oc,
                                       EquatingRecipes::Structures::FormType::NEW,
                                       S,
                                       I);
          oa_pr_sum += oa(k) * pr;
        }
        pr = itemResponseFunctionNRM(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     oa,
                                     oc,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);
        pr_over_S = pr * theta * (oa(categoryIndex) - oa_pr_sum);
        return pr_over_S;
      }

      /*------------------------------------------------------------------------------
      Functionality:
        Under the NR model, calculate partial derivative of P# (from old-to-new
        transformation) with respect to I.
        
        Note:
      oa[1..CatNum] for the discrimination parameters
      oc[1..CatNum] for the intercept parameters
            No scaling constant

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
      double itemResponseFunctionDerivativeNRMNewOverI(const size_t& numberOfCategories,
                                                       const size_t& categoryIndex,
                                                       const double& theta,
                                                       const Eigen::VectorXd& oa,
                                                       const Eigen::VectorXd& oc,
                                                       const double& S,
                                                       const double& I) {
        int k;
        double pr;
        double oa_pr_sum = 0.0;
        double pr_over_I;

        for (k = 1; k <= numberOfCategories; k++) {
          pr = itemResponseFunctionNRM(numberOfCategories,
                                       k,
                                       theta,
                                       oa,
                                       oc,
                                       EquatingRecipes::Structures::FormType::NEW,
                                       S,
                                       I);
          oa_pr_sum += oa(k) * pr;
        }

        pr = itemResponseFunctionNRM(numberOfCategories,
                                     categoryIndex,
                                     theta,
                                     oa,
                                     oc,
                                     EquatingRecipes::Structures::FormType::NEW,
                                     S,
                                     I);
        pr_over_I = pr * (oa(categoryIndex) - oa_pr_sum);
        return pr_over_I;
      }
    };
  } // namespace Implementation
} // namespace EquatingRecipes

#endif