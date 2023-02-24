#ifndef IRT_SCALE_TRANSFORMATION_HPP
#define IRT_SCALE_TRANSFORMATION_HPP

#include <cmath>
#include <iostream>
#include <limits>
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
#include <equating_recipes/structures/irt_equating_results.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_input.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_control.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/irt_model_functions.hpp>

namespace EquatingRecipes {
  class IRTScaleTransformation {
  public:
    void runIRTScaleTransformation(FILE& outf,
                                   std::string& tt,
                                   std::string& ItemNewFile,
                                   std::string& ItemOldFile,
                                   std::string& ItemCommonFile,
                                   std::string& DistNewFile,
                                   std::string& DistOldFile,
                                   int HA,
                                   EquatingRecipes::Structures::Symmetry HAsym,
                                   bool HAfs,
                                   double HAs,
                                   double HAi,
                                   int SL,
                                   EquatingRecipes::Structures::Symmetry SLsym,
                                   bool SLfs,
                                   double SLs,
                                   double SLi,
                                   std::string ST,
                                   int PrintFiles) {}

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
    void ScaleTransform(const std::string& ItemOutF,
                        const std::string& DistOutF,
                        const double& slope,
                        const double& intercept,
                        const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                        const EquatingRecipes::Structures::IRTScaleTransformationControl& handle) {
    }

    std::vector<EquatingRecipes::Structures::ItemSpecification> importItemSpecifications(FILE& inf,
                                                                                         const std::string& oldOrnew,
                                                                                         const EquatingRecipes::Structures::IRTScaleTransformationControl& handle) {
      std::vector<EquatingRecipes::Structures::ItemSpecification> items;

      return items;
    }

    std::vector<EquatingRecipes::Structures::CommonItemSpecification> readCommonItems(FILE& inf,
                                                                                      const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                                                                                      const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                                                                                      const EquatingRecipes::Structures::IRTScaleTransformationControl& handle) {
      std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;

      return commonItems;
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          if (original) {
            prob = irtModelFunctions.itemResponseFunction3PL(categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB(1),
                                                             commonItem.oldC(1),
                                                             "old",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunction3PL(categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB(1),
                                                             commonItem.newC(1),
                                                             "old",
                                                             S,
                                                             I);
          }
          break;
        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionLGR(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB,
                                                             "old",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionLGR(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB,
                                                             "old",
                                                             S,
                                                             I);
          }
          break;
        case pc:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionGPC(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB,
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionGPC(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB,
                                                             "old",
                                                             S,
                                                             I);
          }
          break;
        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionNRM(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.oldA,
                                                             commonItem.oldC,
                                                             "old",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionNRM(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.newA,
                                                             commonItem.newC,
                                                             "old",
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          if (original) {
            prob = irtModelFunctions.itemResponseFunction3PL(categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB(1),
                                                             commonItem.newC(1),
                                                             "new",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunction3PL(categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB(1),
                                                             commonItem.oldC(1),
                                                             "new",
                                                             S,
                                                             I);
          }
          break;

        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionLGR(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB,
                                                             "new",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionLGR(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB,
                                                             "new",
                                                             S,
                                                             I);
          }
          break;

        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionGPC(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.newA(1),
                                                             commonItem.newB,
                                                             "new",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionGPC(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.scaleConstant,
                                                             commonItem.oldA(1),
                                                             commonItem.oldB,
                                                             "new",
                                                             S,
                                                             I);
          }
          break;

        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          if (original) {
            prob = irtModelFunctions.itemResponseFunctionNRM(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.newA,
                                                             commonItem.newC,
                                                             "new",
                                                             1,
                                                             0);
          } else {
            prob = irtModelFunctions.itemResponseFunctionNRM(commonItem.numberOfCategories,
                                                             categoryIndex,
                                                             theta,
                                                             commonItem.oldA,
                                                             commonItem.oldC,
                                                             "new",
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          derivative = irtModelFunctions.itemResponseFunctionDerivative3PLOldOverS(categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB(1),
                                                                                   commonItem.newC(1),
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeLGROldOverS(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeGPCOldOverS(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeNRMOldOverS(commonItem.numberOfCategories,
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          derivative = irtModelFunctions.itemResponseFunctionDerivative3PLOldOverI(categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB(1),
                                                                                   commonItem.newC(1),
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeLGROldOverI(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeGPCOldOverI(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.newA(1),
                                                                                   commonItem.newB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeNRMOldOverI(commonItem.numberOfCategories,
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          derivative = irtModelFunctions.itemResponseFunctionDerivative3PLNewOverS(categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldA(1),
                                                                                   commonItem.oldB(1),
                                                                                   commonItem.oldC(1),
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeLGRNewOverS(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldA(1),
                                                                                   commonItem.oldB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeGPCNewOverS(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldwA(1),
                                                                                   commonItem.oldwB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeNRMNewOverS(commonItem.numberOfCategories,
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

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      switch (commonItem.irtModel) {
        case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
          derivative = irtModelFunctions.itemResponseFunctionDerivative3PLNewOverI(categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldA(1),
                                                                                   commonItem.oldB(1),
                                                                                   commonItem.oldC(1),
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeLGRNewOverI(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldA(1),
                                                                                   commonItem.oldB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeGPCNewOverI(commonItem.numberOfCategories,
                                                                                   categoryIndex,
                                                                                   theta,
                                                                                   commonItem.scaleConstant,
                                                                                   commonItem.oldA(1),
                                                                                   commonItem.oldB,
                                                                                   S,
                                                                                   I);

          break;

        case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
          derivative = irtModelFunctions.itemResponseFunctionDerivativeNRMNewOverI(commonItem.numberOfCategories,
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

    void StHaebara(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                   const std::vector<EquatingRecipes::Structures::CommonItemSpecification>& commonItems,
                   const EquatingRecipes::Structures::Symmetry& symmetry,
                   const bool& funcStd,
                   const double& S0,
                   const double& I0,
                   Eigen::VectorXd& slope,
                   Eigen::VectorXd& intercept) {}

    double FuncHaebara(const Eigen::VectorXd& x) {
      return 0;
    }

    void GradHaebara(const Eigen::VectorXd& x,
                     Eigen::VectorXd& grad) {}

    void StMeanMean(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                    const std::vector<EquatingRecipes::Structures::CommonItemSpecification>& commonItems,
                    const double& slope,
                    const double& intercept) {}

    void StMeanSigma(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                     const std::vector<EquatingRecipes::Structures::CommonItemSpecification>& commonItems,
                     const double& slope,
                     const double& intercept);

    void StStockingLord(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                        const std::vector<EquatingRecipes::Structures::CommonItemSpecification>& commonItems,
                        const EquatingRecipes::Structures::Symmetry& symmetry,
                        const bool& funcStd,
                        const double& S0,
                        const double& I0,
                        const Eigen::VectorXd& slope,
                        const Eigen::VectorXd& intercept);
    double funcStockingLord(const Eigen::VectorXd& x);
    void gradStockingLord(const Eigen::VectorXd& x,
                          Eigen::VectorXd& grad);

    void readThetaInfo(FILE& inf,
                       const std::string& oldOrnew,
                       const EquatingRecipes::Structures::IRTScaleTransformationControl& handle) {}
  };
} // namespace EquatingRecipes

#endif