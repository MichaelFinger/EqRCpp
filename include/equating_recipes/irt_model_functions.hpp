#ifndef IRT_MODEL_FUNCTIONS_HPP
#define IRT_MODEL_FUNCTIONS_HPP

#include <cmath>

#include <Eigen/Dense>

#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/form_type.hpp>

namespace EquatingRecipes {
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
                                                                item.scaleConstant,
                                                                item.a(1),
                                                                item.b(1),
                                                                item.c(1),
                                                                "old",
                                                                1,
                                                                0);

          break;
        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          prob = itemResponseFunctionLGR(item.numberOfCategories,
                                                                categoryIndex,
                                                                theta,
                                                                item.scaleConstant,
                                                                item.a(1),
                                                                item.b,
                                                                "old",
                                                                1,
                                                                0);
          break;
        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          prob = itemResponseFunctionGPC(item.numberOfCategories,
                                                                categoryIndex,
                                                                theta,
                                                                item.scaleConstant,
                                                                item.a(1),
                                                                item.b,
                                                                "old",
                                                                1,
                                                                0);
          break;
        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          prob = itemResponseFunctionNRM(item.numberOfCategories,
                                                                categoryIndex,
                                                                theta,
                                                                item.a,
                                                                item.c,
                                                                "old",
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
                                                         item.scaleConstant,
                                                         item.a(1),
                                                         item.b(1),
                                                         item.c(1));
          break;
        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          derivative = itemResponseFunctionDerivativeLGR(item.numberOfCategories,
                                                         categoryIndex,
                                                         theta,
                                                         item.scaleConstant,
                                                         item.a(1),
                                                         item.b);
          break;
        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          derivative = itemResponseFunctionDerivativeGPC(item.numberOfCategories,
                                                         categoryIndex,
                                                         theta,
                                                         item.scaleConstant,
                                                         item.a(1),
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

      if (categoryIndex == 0) {
        derivative = -1.0 * D * a * (1.0 - c) * uj / ((1.0 + uj) * (1.0 + uj));
      } else if (categoryIndex == 1) {
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

      if (categoryIndex == 0) {
        cp_jk = 1.0;
        cp_jk1 = 1.0 / (1.0 + std::exp(-1.0 * D * a * (theta - b(categoryIndex + 1))));
      } else {
        if (categoryIndex < numberOfCategories - 1) {
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

      for (size_t index = 0; index < numberOfCategories; index++) {
        double aSum = static_cast<double>(index + 1) * D * a;
        double bSum = -1.0 * D * a * b(Eigen::seq(1, b.size() - 1)).sum(); /* b[0] = 0 */
        double ujk = std::exp(aSum * theta + bSum);
        vj += ujk;
        aUSum += aSum * ujk;
      }

      double aSum = static_cast<double>(categoryIndex + 1) * D * a;
      double bSum = -1.0 * D * a * b(Eigen::seq(1, categoryIndex)).sum(); /* b[0] = 0 */
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

      for (size_t index = 0; index < numberOfCategories; index++) {
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
      // double as, bs, cs; /* parameters from new-to-old transformation */
      // double ar, br, cr; /* parameters from old-to-new transformation */
      // double devs, devr;

      double aPar;
      double bPar;
      double cPar;

      if (scale == EquatingRecipes::Structures::FormType::OLD) {
        /* new-to-old scale transformation */
        aPar = a/S;
        bPar = S*b + I;
        cPar = c;
      }
      else {
        /* old-to-new scale transformation */
        aPar = S*a;
        bPar = (b-I)/S;
        cPar = c;
      }

      double deviate = D*aPar*(theta-bPar);

      double prob =  cPar + (1.0 - cPar) / (1.0 + std::exp(-1.0 * deviate));

        if (categoryIndex == 0) {
          prob = 1.0 - prob;
        }

      return prob;
    }

    double itemResponseFunctionDerivative3PLOldOverS(const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const double& nb,
                                                     const double& nc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivative3PLOldOverI(const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const double& nb,
                                                     const double& nc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivative3PLNewOverS(const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const double& ob,
                                                     const double& oc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivative3PLNewOverI(const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const double& ob,
                                                     const double& oc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double cumulativeResponseProbabilityLGR(const size_t& categoryIndex,
                                            const double& theta,
                                            const double& D,
                                            const double& a,
                                            const double& b,
                                            const EquatingRecipes::Structures::FormType& scale,
                                            const double& S,
                                            const double& I) {
      return 0;
    }
    double itemResponseFunctionLGR(const size_t& numberOfCategories,
                                   const size_t categoryIndex,
                                   const double& theta,
                                   const double& D,
                                   const double& a,
                                   const Eigen::VectorXd& b,
                                   const std::string& scale,
                                   const double& S,
                                   const double& I) {
      return 0;
    }

    double itemResponseFunctionDerivativeLGROldOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const Eigen::VectorXd& nb,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeLGROldOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const Eigen::VectorXd& nb,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeLGRNewOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const Eigen::VectorXd& ob,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeLGRNewOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const Eigen::VectorXd& ob,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionGPC(const size_t& numberOfCategories,
                                   const size_t categoryIndex,
                                   const double& theta,
                                   const double& D,
                                   const double& a,
                                   const Eigen::VectorXd& b,
                                   const EquatingRecipes::Structures::FormType& scale,
                                   const double& S,
                                   const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeGPCOldOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const Eigen::VectorXd& nb,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeGPCOldOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& na,
                                                     const Eigen::VectorXd& nb,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeGPCNewOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const Eigen::VectorXd& ob,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeGPCNewOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const double& D,
                                                     const double& oa,
                                                     const Eigen::VectorXd& ob,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionNRM(const size_t& numberOfCategories,
                                   const size_t& categoryIndex,
                                   const double& theta,
                                   const Eigen::VectorXd& a,
                                   const Eigen::VectorXd& c,
                                   const EquatingRecipes::Structures::FormType& scale,
                                   const double& S,
                                   const double& I) {
      return 0;
    }

    double itemResponseFunctionDerivativeNRMOldOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const Eigen::VectorXd& na,
                                                     const Eigen::VectorXd& nc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeNRMOldOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const Eigen::VectorXd& na,
                                                     const Eigen::VectorXd& nc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeNRMNewOverS(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const Eigen::VectorXd& oa,
                                                     const Eigen::VectorXd& oc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
    double itemResponseFunctionDerivativeNRMNewOverI(const size_t& numberOfCategories,
                                                     const size_t& categoryIndex,
                                                     const double& theta,
                                                     const Eigen::VectorXd& oa,
                                                     const Eigen::VectorXd& oc,
                                                     const double& S,
                                                     const double& I) {
      return 0;
    }
  };
}

#endif