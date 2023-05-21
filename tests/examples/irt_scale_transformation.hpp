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

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      class IRTScaleTransformation {
      public:
        void operator()() {
          EquatingRecipes::Tests::Examples::Datasets::ItemSpecsFile dummyXItems;
          dummyXItems.import("dummyXItems.txt");

          EquatingRecipes::Tests::Examples::Datasets::ItemSpecsFile dummyYItems;
          dummyYItems.import("dummyYItems.txt");

          EquatingRecipes::Tests::Examples::Datasets::DummyV3Items dummyV3Items;
          dummyV3Items.import();

          EquatingRecipes::Tests::Examples::Datasets::QuadratureFile dummyXDist;
          dummyXDist.import("dummyXdist.txt");

          EquatingRecipes::Tests::Examples::Datasets::QuadratureFile dummyYDist;
          dummyYDist.import("dummyYdist.txt");

          std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems = createCommonItemSpecs(dummyYItems.itemSpecs,
                                                                                                                dummyXItems.itemSpecs,
                                                                                                                dummyV3Items.itemPairSpecs);

          EquatingRecipes::Structures::IRTScaleTransformationData irtScaleTransformationData;

          // irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_MEAN);
          // irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_SIGMA);
          irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA);
          // irtScaleTransformationData.irtScaleTranformationMethods.insert(EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD);

          irtScaleTransformationData.commonItems = commonItems;

          irtScaleTransformationData.interceptStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 0;
          irtScaleTransformationData.slopeStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = 1;

          // irtScaleTransformationData.interceptStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 0;
          // irtScaleTransformationData.slopeStartingValue[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = 1;
          
          // irtScaleTransformationData.maximumNumberOfIterations = 100;
          // irtScaleTransformationData.maximumAbsoluteChangeInFunctionValue = 0.0001;
          irtScaleTransformationData.maximumRelativeChangeInFunctionValue = 0.0000000001;
          // irtScaleTransformationData.maximumAbsoluteChangeInParameterValues = 0.0001;
          // irtScaleTransformationData.maximumRelativeChangeInParameterValues = 1e-8;

          irtScaleTransformationData.newItems = dummyXItems.itemSpecs;
          irtScaleTransformationData.oldItems = dummyYItems.itemSpecs;

          irtScaleTransformationData.quadratureNewForm = dummyXDist.quadrature;
          irtScaleTransformationData.quadratureOldForm = dummyYDist.quadrature;
          
          irtScaleTransformationData.standardizations[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = true;
          irtScaleTransformationData.symmetryOptions[EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA] = EquatingRecipes::Structures::Symmetry::SYMMETRIC;
          
          // irtScaleTransformationData.standardizations[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = true;
          // irtScaleTransformationData.symmetryOptions[EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD] = EquatingRecipes::Structures::Symmetry::SYMMETRIC;

          EquatingRecipes::Analyses::IRTScaleTransformation irtScaleTransformation;
          std::string title = "IRT Scale Transformation";
          irtScaleTransformation(title, irtScaleTransformationData);
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