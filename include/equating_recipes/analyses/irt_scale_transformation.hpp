#ifndef ANALYSES_STOCKING_LORD_HPP
#define ANALYSES_STOCKING_LORD_HPP

#include <map>
#include <optional>
#include <set>
#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/irt_scale_transformation.hpp>
#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <equating_recipes/structures/irt_scale_transformation_item_results.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/json/structures.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct IRTScaleTransformation {
      nlohmann::json operator()(std::string& title,
                                EquatingRecipes::Structures::IRTScaleTransformationData& irtScaleTransformationData) {
        EquatingRecipes::IRTScaleTransformation irtScaleTransformation;

        irtScaleTransformation.run(irtScaleTransformationData);

        nlohmann::json j = nlohmann::json::object();
        j["analysis_title"] = title;
        j["analysis_type"] = "irt_scale_transformation";
        j["analysis_results"] = irtScaleTransformationData;

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif