/* 
  From Source: IRTst.h
  Original Struct: ITEM_SPEC
  Description: 
*/

#ifndef STRUCTURES_ITEM_SPECIFICATION_HPP
#define STRUCTURES_ITEM_SPECIFICATION_HPP

#include <optional>
#include <Eigen/Core>

#include <equating_recipes/structures/irt_model.hpp>

namespace EquatingRecipes {
  namespace Structures {
    class ItemSpecification {
    public:
      int itemID;
      unsigned long numberOfCategories;
      double scalingConstant;
      Eigen::VectorXd scoringFunctionValues;
      Eigen::VectorXd a;
      Eigen::VectorXd b;
      Eigen::VectorXd c;
      Eigen::VectorXd d;
      EquatingRecipes::Structures::IRTModel irtModel;

      static ItemSpecification buildOneParameterLogistic(const int& itemID,
                                                         const double& itemDifficulty,
                                                         const double& scalingConstant,
                                                         const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                         const bool& locationParameterization) {
        ItemSpecification itemSpec = ItemSpecification::buildDichotomousParameterLogistic(itemID,
                                                                                          1.0,
                                                                                          itemDifficulty,
                                                                                          0.0,
                                                                                          scalingConstant,
                                                                                          scoringFunctionValues,
                                                                                          locationParameterization);
        return itemSpec;
      }

      static ItemSpecification buildTwoParameterLogistic(const int& itemID,
                                                         const double& itemDiscrimination,
                                                         const double& itemDifficulty,
                                                         const double& scalingConstant,
                                                         const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                         const bool& locationParameterization) {
        ItemSpecification itemSpec = ItemSpecification::buildDichotomousParameterLogistic(itemID,
                                                                                          itemDiscrimination,
                                                                                          itemDifficulty,
                                                                                          0.0,
                                                                                          scalingConstant,
                                                                                          scoringFunctionValues,
                                                                                          locationParameterization);

        return itemSpec;
      }

      static ItemSpecification buildThreeParameterLogistic(const int& itemID,
                                                           const double& itemDiscrimination,
                                                           const double& itemDifficulty,
                                                           const double& lowerAsymptote,
                                                           const double& scalingConstant,
                                                           const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                           const bool& locationParameterization) {
        ItemSpecification itemSpec = ItemSpecification::buildDichotomousParameterLogistic(itemID,
                                                                                          itemDiscrimination,
                                                                                          itemDifficulty,
                                                                                          lowerAsymptote,
                                                                                          scalingConstant,
                                                                                          scoringFunctionValues,
                                                                                          locationParameterization);
        return itemSpec;
      }

      static ItemSpecification buildGradedResponse(const int& itemID,
                                                   const double& itemDiscrimination,
                                                   const Eigen::VectorXd& thresholdParameters,
                                                   const double& scalingConstant,
                                                   const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                   const bool& locationParameterization) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE;
        itemSpec.numberOfCategories = thresholdParameters.size() + 1;
        itemSpec.scalingConstant = scalingConstant;
        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);
        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.b.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(2) = itemDiscrimination;

        if (locationParameterization) {
          itemSpec.b(Eigen::seq(2, itemSpec.numberOfCategories)) = thresholdParameters;
        } else {
          itemSpec.b(Eigen::seq(2, itemSpec.numberOfCategories)) = -1.0 * thresholdParameters / itemDiscrimination;
        }

        return itemSpec;
      }

      static ItemSpecification buildGeneralizedPartialCredit(const int& itemID,
                                                             const double& itemDiscrimination,
                                                             const Eigen::VectorXd& stepParameters,
                                                             const double& scalingConstant,
                                                             const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                             const bool& locationParameterization) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT;
        itemSpec.numberOfCategories = stepParameters.size() + 1;
        itemSpec.scalingConstant = scalingConstant;
        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);
        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.b.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(2) = itemDiscrimination;

        if (locationParameterization) {
          itemSpec.b(Eigen::seq(1, itemSpec.numberOfCategories)) = stepParameters;
        } else {
          itemSpec.b(Eigen::seq(1, itemSpec.numberOfCategories)) = -1.0 * stepParameters / itemDiscrimination;
        }

        return itemSpec;
      }

      static ItemSpecification buildGradedResponseRatingScale(const int& itemID,
                                                              const double& itemDiscrimination,
                                                              const double& itemLocation,
                                                              const Eigen::VectorXd& categoryParameters,
                                                              const double& scalingConstant,
                                                              const std::optional<Eigen::VectorXd>& scoringFunctionValues) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE;
        itemSpec.numberOfCategories = categoryParameters.size();
        itemSpec.scalingConstant = scalingConstant;
        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);

        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.b.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(2) = itemDiscrimination;

        itemSpec.b(Eigen::seq(2, itemSpec.numberOfCategories)) =
            Eigen::VectorXd::Constant(itemSpec.numberOfCategories, itemLocation) - categoryParameters;

        return itemSpec;
      }

      static ItemSpecification buildGeneralizedPartialCreditRatingScale(const int& itemID,
                                                                        const double& itemDiscrimination,
                                                                        const double& itemLocation,
                                                                        const Eigen::VectorXd& categoryParameters,
                                                                        const double& scalingConstant,
                                                                        const std::optional<Eigen::VectorXd>& scoringFunctionValues) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT;
        itemSpec.numberOfCategories = categoryParameters.size();
        itemSpec.scalingConstant = scalingConstant;
        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);

        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.b.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(2) = itemDiscrimination;
        itemSpec.b(0) = itemLocation;

        itemSpec.b(Eigen::seq(1, itemSpec.numberOfCategories)) =
            Eigen::VectorXd::Constant(itemSpec.numberOfCategories, itemSpec.b(0)) - categoryParameters;

        return itemSpec;
      }

      static ItemSpecification buildNominalResponse(const int& itemID,
                                                    const Eigen::VectorXd& a,
                                                    const Eigen::VectorXd& c,
                                                    const std::optional<Eigen::VectorXd>& scoringFunctionValues) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE;
        itemSpec.numberOfCategories = a.size();
        itemSpec.scalingConstant = 1.0;
        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);

        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.c.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(Eigen::seq(1, itemSpec.numberOfCategories)) = a;
        itemSpec.c(Eigen::seq(1, itemSpec.numberOfCategories)) = c;

        return itemSpec;
      }

    private:
      static ItemSpecification buildDichotomousParameterLogistic(const int& itemID,
                                                                 const double& itemDiscrimination,
                                                                 const double& itemDifficulty,
                                                                 const double& lowerAsymptote,
                                                                 const double& scalingConstant,
                                                                 const std::optional<Eigen::VectorXd>& scoringFunctionValues,
                                                                 const bool& locationParameterization) {
        ItemSpecification itemSpec;

        itemSpec.itemID = itemID;
        itemSpec.numberOfCategories = 2;
        itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC;
        itemSpec.a.resize(itemSpec.numberOfCategories + 1);
        itemSpec.b.resize(itemSpec.numberOfCategories + 1);
        itemSpec.c.resize(itemSpec.numberOfCategories + 1);

        itemSpec.a(itemSpec.numberOfCategories) = itemDiscrimination;

        if (locationParameterization) {
          itemSpec.b(itemSpec.numberOfCategories) = itemDifficulty;
        } else {
          itemSpec.b(itemSpec.numberOfCategories) = -1.0 * itemDifficulty / itemDiscrimination;
        }

        itemSpec.c(itemSpec.numberOfCategories) = lowerAsymptote;

        itemSpec.scalingConstant = scalingConstant;

        itemSpec.scoringFunctionValues = ItemSpecification::getScoringFunctionValues(itemSpec.numberOfCategories,
                                                                                     scoringFunctionValues);

        return itemSpec;
      }

      static Eigen::VectorXd getScoringFunctionValues(const size_t& numberOfCategories,
                                                      const std::optional<Eigen::VectorXd>& scoringFunctionValues) {
        Eigen::VectorXd result(numberOfCategories + 1);

        for (size_t index = 1; index <= numberOfCategories; index++) {
          if (scoringFunctionValues.has_value()) {
            result(index) = scoringFunctionValues.value()(index);
          } else {
            result(index) = static_cast<double>(index - 1);
          }
        }

        return result;
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif