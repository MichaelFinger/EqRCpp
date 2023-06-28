#ifndef ANALYSES_KERNEL_EQUATING_RG_HPP
#define ANALYSES_KERNEL_EQUATING_RG_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/implementation/kernel_equating.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace KernelEquating {
      struct RandomGroupsEquating {
        struct InputData {

        };

        struct OutputData {

        };

        nlohmann::json operator()(const InputData& inputData,
        OutputData& outputData) {
          nlohmann::json j;
          return j;
        }
      };
    }
  }
}

#endif