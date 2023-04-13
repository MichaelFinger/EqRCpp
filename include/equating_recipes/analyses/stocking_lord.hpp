#ifndef ANALYSES_STOCKING_LORD_HPP
#define ANALYSES_STOCKING_LORD_HPP

#include <string>
#include <Eigen/Core>
#include <equating_recipes/irt_scale_transformation.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct StockingLord {
      struct InputData {
        std::string title;
        std::string datasetName;
        std::string variableName;

        Eigen::VectorXd scoreFrequencies;
        double minimumScore;
        double maximumScore;
        double scoreIncrement;
        std::string id = "X";
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::StockingLord::InputData& inputData,
                                EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
        irtScaleTransformationData.irtScaleTranformationMethod = EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD;

        nlohmann::json j = {{"analysis_title", inputData.title},
                            {"analysis_type", "stocking_lord"},
                            {"analysis_results", irtScaleTransformationData}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif