#ifndef TESTS_EXAMPLES_DATASETS_ITEM_SPECS_FILE_HPP
#define TESTS_EXAMPLES_DATASETS_ITEM_SPECS_FILE_HPP

#include <iostream>
#include <istream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

#include <boost/algorithm/string.hpp>

#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/item_specification.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        class ItemSpecsFile {
        public:
          size_t numberOfItems;
          std::vector<EquatingRecipes::Structures::ItemSpecification> itemSpecs;

          void import(const std::string& filename) {
            // itemSpecs = importItemSpecsFile(filename);
            itemSpecs = importItemSpecsFileAlt(filename);
          }

        private:
          std::vector<EquatingRecipes::Structures::ItemSpecification> importItemSpecsFile(const std::string& filename) {
            std::vector<std::string> linesRead;

            std::ifstream ifs(filename,
                              std::ios::in);

            for (std::string lineRead; std::getline(ifs, lineRead);) {
              linesRead.push_back(lineRead);
            }

            ifs.close();

            size_t lineIndex = 0;

            numberOfItems = static_cast<size_t>(std::stoi(linesRead[0]));

            for (size_t itemIndex = 0; itemIndex < numberOfItems; itemIndex++) {
              lineIndex++;

              size_t fieldIndex = 0;

              std::vector<std::string> values = splitString(linesRead[lineIndex]);

              EquatingRecipes::Structures::ItemSpecification itemSpec;
              itemSpec.itemID = std::stoi(values[fieldIndex]);

              fieldIndex++;

              if (values[fieldIndex] == "L3") {
                itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC;
              } else if (values[fieldIndex] == "GR") {
                itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE;
              } else if (values[fieldIndex] == "PC") {
                itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT;
              } else if (values[fieldIndex] == "NR") {
                itemSpec.irtModel = EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE;
              }

              fieldIndex++;

              itemSpec.numberOfCategories = static_cast<unsigned long>(std::stoi(values[fieldIndex]));
              fieldIndex++;

              itemSpec.scoringFunctionValues.resize(itemSpec.numberOfCategories);

              for (size_t respIndex = 0; respIndex < itemSpec.numberOfCategories; respIndex++) {
                itemSpec.scoringFunctionValues(respIndex) = std::stod(values[fieldIndex]);

                fieldIndex++;
              }

              itemSpec.scalingConstant = std::stod(values[fieldIndex]);
              fieldIndex++;

              itemSpec.a.setZero(itemSpec.numberOfCategories + 1);
              itemSpec.b.setZero(itemSpec.numberOfCategories + 1);
              itemSpec.c.setZero(itemSpec.numberOfCategories + 1);
              itemSpec.d.setZero(itemSpec.numberOfCategories + 1);

              if (itemSpec.irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                itemSpec.a(2) = std::stod(values[fieldIndex]);
                fieldIndex++;

                itemSpec.b(2) = std::stod(values[fieldIndex]);
                fieldIndex++;

                itemSpec.c(2) = std::stod(values[fieldIndex]);
                fieldIndex++;
              } else if (itemSpec.irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE) {
                itemSpec.a(2) = std::stod(values[fieldIndex]);
                fieldIndex++;

                for (size_t respIndex = 2; respIndex <= itemSpec.numberOfCategories; respIndex++) {
                  itemSpec.b(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }
              } else if (itemSpec.irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
                itemSpec.a(2) = std::stod(values[fieldIndex]);
                fieldIndex++;

                itemSpec.b(0) = std::stod(values[fieldIndex]); // location
                fieldIndex++;

                for (size_t respIndex = 1; respIndex <= itemSpec.numberOfCategories; respIndex++) {
                  itemSpec.d(respIndex) = std::stod(values[fieldIndex]);
                  itemSpec.b(respIndex) = itemSpec.b(0) - itemSpec.d(respIndex);

                  fieldIndex++;
                }
              } else if (itemSpec.irtModel == EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE) {
                itemSpec.scalingConstant = 1.0;

                for (size_t respIndex = 1; respIndex <= itemSpec.numberOfCategories; respIndex++) {
                  itemSpec.a(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }

                for (size_t respIndex = 1; respIndex <= itemSpec.numberOfCategories; respIndex++) {
                  itemSpec.c(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }
              }

              itemSpecs.push_back(itemSpec);
            }

            return itemSpecs;
          }

          std::vector<EquatingRecipes::Structures::ItemSpecification> importItemSpecsFileAlt(const std::string& filename) {
            std::vector<std::string> linesRead;

            std::ifstream ifs(filename,
                              std::ios::in);

            for (std::string lineRead; std::getline(ifs, lineRead);) {
              linesRead.push_back(lineRead);
            }

            ifs.close();

            size_t lineIndex = 0;

            numberOfItems = static_cast<size_t>(std::stoi(linesRead[0]));

            for (size_t itemIndex = 0; itemIndex < numberOfItems; itemIndex++) {
              lineIndex++;

              size_t fieldIndex = 0;

              std::vector<std::string> values = splitString(linesRead[lineIndex]);

              int itemID = std::stoi(values[fieldIndex]);
              fieldIndex++;

              EquatingRecipes::Structures::IRTModel irtModel;

              if (values[fieldIndex] == "L3") {
                irtModel = EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC;
              } else if (values[fieldIndex] == "GR") {
                irtModel = EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE;
              } else if (values[fieldIndex] == "PC") {
                irtModel = EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT;
              } else if (values[fieldIndex] == "NR") {
                irtModel = EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE;
              }
              fieldIndex++;

              unsigned long numberOfCategories = static_cast<unsigned long>(std::stoi(values[fieldIndex]));
              fieldIndex++;

              Eigen::VectorXd scoringFunctionValues(numberOfCategories + 1);
              for (size_t respIndex = 1; respIndex <= numberOfCategories; respIndex++) {
                scoringFunctionValues(respIndex) = std::stod(values[fieldIndex]);
                fieldIndex++;
              }

              double scalingConstant = std::stod(values[fieldIndex]);
              fieldIndex++;

              EquatingRecipes::Structures::ItemSpecification itemSpec;

              if (irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                double a = std::stod(values[fieldIndex]);
                fieldIndex++;

                double b = std::stod(values[fieldIndex]);
                fieldIndex++;

                double c = std::stod(values[fieldIndex]);
                fieldIndex++;

                itemSpec = EquatingRecipes::Structures::ItemSpecification::buildThreeParameterLogistic(itemID,
                                                                                                       a,
                                                                                                       b,
                                                                                                       c,
                                                                                                       scalingConstant,
                                                                                                       scoringFunctionValues);
              } else if (irtModel == EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE) {
                double a = std::stod(values[fieldIndex]);
                fieldIndex++;

                Eigen::VectorXd b(numberOfCategories - 1);
                for (size_t respIndex = 0; respIndex < numberOfCategories - 1; respIndex++) {
                  b(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }

                itemSpec = EquatingRecipes::Structures::ItemSpecification::buildGradedResponse(itemID,
                                                                                               a,
                                                                                               b,
                                                                                               scalingConstant,
                                                                                               scoringFunctionValues);
              } else if (irtModel == EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT) {
                double a = std::stod(values[fieldIndex]);
                fieldIndex++;

                double b = std::stod(values[fieldIndex]);
                fieldIndex++;

                Eigen::VectorXd d(numberOfCategories);
                for (size_t respIndex = 0; respIndex < numberOfCategories - 1; respIndex++) {
                  d(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }

                itemSpec = EquatingRecipes::Structures::ItemSpecification::buildGeneralizedPartialCreditRatingScale(itemID,
                                                                                                                    a,
                                                                                                                    b,
                                                                                                                    d,
                                                                                                                    scalingConstant,
                                                                                                                    scoringFunctionValues);
              } else if (irtModel == EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE) {
                scalingConstant = 1.0;
                Eigen::VectorXd a(numberOfCategories);
                Eigen::VectorXd c(numberOfCategories);

                for (size_t respIndex = 0; respIndex < numberOfCategories; respIndex++) {
                  a(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }

                for (size_t respIndex = 0; respIndex < numberOfCategories; respIndex++) {
                  c(respIndex) = std::stod(values[fieldIndex]);
                  fieldIndex++;
                }

                itemSpec = EquatingRecipes::Structures::ItemSpecification::buildNominalResponse(itemID,
                                                                                                a,
                                                                                                c,
                                                                                                scoringFunctionValues);
              }

              itemSpecs.push_back(itemSpec);
            }

            return itemSpecs;
          }

          std::vector<std::string> splitString(const std::string& lineRead) {
            std::vector<std::string> values;

            boost::split(values, lineRead, boost::is_any_of(" "), boost::algorithm::token_compress_on);

            return values;
          }
        };
      } // namespace Datasets
    }   // namespace Examples
  }     // namespace Tests
} // namespace EquatingRecipes

#endif