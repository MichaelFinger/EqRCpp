#ifndef TESTS_EXAMPLES_IRT_SCALE_TRANSFORMATION_HPP
#define TESTS_EXAMPLES_IRT_SCALE_TRANSFORMATION_HPP

#include <string>
#include <boost/math/distributions/normal.hpp>

#include <datasets/dummyV3Items.hpp>
#include <datasets/item_specs_file.hpp>
#include <datasets/quadrature_file.hpp>
#include <equating_recipes/analyses/irt_scale_transformation.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <datasets/item_pair_file.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      class IRTScaleTransformation {
      public:
        void operator()() {
          std::string newItemsFilename = "LSAT6XItems";
          std::string oldItemsFilename = "LSAT6YItems";
          std::string itemPairsFilename = "LSAT6VItems";
          std::string newQuadDistFilename = "LSAT6Xdist";
          std::string oldQuadDistFilename = "LSAT6Ydist";

          EquatingRecipes::Tests::Examples::Datasets::ItemSpecsFile newItems;
          newItems.import(newItemsFilename);

          EquatingRecipes::Tests::Examples::Datasets::ItemSpecsFile oldItems;
          oldItems.import(oldItemsFilename);

          EquatingRecipes::Tests::Examples::Datasets::ItemPairFile itemPairs;
          itemPairs.import(itemPairsFilename);

          EquatingRecipes::Tests::Examples::Datasets::QuadratureFile oldDist;
          oldDist.import(oldQuadDistFilename);

          EquatingRecipes::Tests::Examples::Datasets::QuadratureFile newDist;
          newDist.import(newQuadDistFilename);

          std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems = createCommonItemSpecs(oldItems.itemSpecs,
                                                                                                                newItems.itemSpecs,
                                                                                                                itemPairs.itemPairSpecs);

          EquatingRecipes::Structures::IRTScaleTransformationData irtScaleTransformationData;

          irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_MEAN);
          irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_SIGMA);
          irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA);
          irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD);

          irtScaleTransformationData.commonItems = commonItems;

          irtScaleTransformationData.interceptStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 0;
          irtScaleTransformationData.slopeStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 1;

          irtScaleTransformationData.interceptStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 0;
          irtScaleTransformationData.slopeStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 1;
          
          irtScaleTransformationData.maximumNumberOfIterations[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 1000;
          irtScaleTransformationData.maximumNumberOfIterations[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 1000;
          
          irtScaleTransformationData.maximumRelativeChangeInParameterValues[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 1e-8;
          irtScaleTransformationData.maximumRelativeChangeInParameterValues[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 1e-8;

          irtScaleTransformationData.newItems = newItems.itemSpecs;
          irtScaleTransformationData.oldItems = oldItems.itemSpecs;

          irtScaleTransformationData.quadratureNewForm = newDist.quadrature;
          irtScaleTransformationData.quadratureOldForm = oldDist.quadrature;
          
          irtScaleTransformationData.standardizations[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = true;
          irtScaleTransformationData.symmetryOptions[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = EquatingRecipes::Structures::Symmetry::SYMMETRIC;
          
          irtScaleTransformationData.standardizations[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = true;
          irtScaleTransformationData.symmetryOptions[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = EquatingRecipes::Structures::Symmetry::OLD_SCALE;

          EquatingRecipes::Analyses::IRTScaleTransformation irtScaleTransformation;
          std::string title = "IRT Scale Transformation";
          nlohmann::json irtScaleTransformationResults = irtScaleTransformation(title, irtScaleTransformationData);

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(irtScaleTransformationResults);
          jsonDoc.toTextFile("irtScaleTransformation.json");
        }

      private:
        std::vector<EquatingRecipes::Structures::CommonItemSpecification> createCommonItemSpecs(const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldFormItemSpecs,
                                                                                                const std::vector<EquatingRecipes::Structures::ItemSpecification>& newFormItemSpecs,
                                                                                                const std::vector<std::pair<size_t, size_t>>& itemPairSpecs) {
          std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;

          std::for_each(itemPairSpecs.begin(),
                        itemPairSpecs.end(),
                        [&](const std::pair<size_t, size_t>& itemPairSpec) {
                          auto iterOld = std::find_if(oldFormItemSpecs.begin(),
                                                      oldFormItemSpecs.end(),
                                                      [&](const EquatingRecipes::Structures::ItemSpecification& itemSpec) {
                                                        return itemSpec.itemID == itemPairSpec.first;
                                                      });

                          auto iterNew = std::find_if(newFormItemSpecs.begin(),
                                                      newFormItemSpecs.end(),
                                                      [&](const EquatingRecipes::Structures::ItemSpecification& itemSpec) {
                                                        return itemSpec.itemID == itemPairSpec.first;
                                                      });

                          if (iterOld == oldFormItemSpecs.end() ||
                              iterNew == newFormItemSpecs.end()) {
                            return;
                          }

                          EquatingRecipes::Structures::CommonItemSpecification commonItem;

                          commonItem.irtModel = iterOld->irtModel;
                          commonItem.numberOfCategories = iterOld->numberOfCategories;
                          commonItem.scalingConstant = iterOld->scalingConstant;
                          commonItem.scoringFunctionValues = iterOld->scoringFunctionValues;

                          commonItem.oldA = iterOld->a;
                          commonItem.oldB = iterOld->b;
                          commonItem.oldC = iterOld->c;
                          commonItem.oldID = iterOld->itemID;

                          commonItem.newA = iterNew->a;
                          commonItem.newB = iterNew->b;
                          commonItem.newC = iterNew->c;
                          commonItem.newID = iterNew->itemID;

                          commonItems.push_back(commonItem);
                        });
          return commonItems;
        }

        EquatingRecipes::Structures::Quadrature getQuadrature() {
          EquatingRecipes::Structures::Quadrature quadrature;
          
          boost::math::normal_distribution<double> normalDist(0.0, 1.0);
          size_t numberOfQuadraturePoints = 31;
          double minTheta = -4;
          double maxTheta = 4;
          double increment = (maxTheta - minTheta) / static_cast<double>(numberOfQuadraturePoints - 1);

          quadrature.thetaValues.resize(numberOfQuadraturePoints);
          quadrature.thetaWeights.resize(numberOfQuadraturePoints);
          for (size_t quadIndex = 0; quadIndex < numberOfQuadraturePoints; quadIndex++) {
            quadrature.thetaValues(quadIndex) = minTheta + quadIndex * increment;
            quadrature.thetaWeights(quadIndex) = boost::math::pdf(normalDist, quadrature.thetaValues(quadIndex));
          }

          quadrature.thetaWeights = quadrature.thetaWeights.cwiseQuotient(Eigen::VectorXd::Constant(numberOfQuadraturePoints, quadrature.thetaWeights.sum()));

          return quadrature;
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif