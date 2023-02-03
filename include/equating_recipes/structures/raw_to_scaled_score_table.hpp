/* 
  From Source: ERutilities.h 
  Original Struct: ESS_RESULTS
  Description: raw to scaled score table
*/

#ifndef STRUCTURES_RAW_TO_SCALED_SCORE_TABLE_HPP
#define STRUCTURES_RAW_TO_SCALED_SCORE_TABLE_HPP

#include <map>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>

#include <equating_recipes/structures/method.hpp>

namespace EquatingRecipes {
  namespace Structures {
    class RawToScaledScoreTable {
    public:
      struct Entry {
        size_t scoreLocation;
        double rawScore;
        double scaledScore;
      };

      Entry getEntry(const size_t& scoreLocation) const {
        Entry entry = lookup.at(scoreLocation);
        return entry;
      }

      Entry getEntry(const double& rawScore) const {
        auto iter = std::find_if(lookup.begin(),
                                 lookup.end(),
                                 [&](const std::pair<size_t, Entry>& entry) {
                                   return rawScore == entry.second.rawScore;
                                 });

        Entry entry = iter->second;

        return entry;
      }

      Entry getFirstEntry() const {
        std::pair<size_t, Entry> firstEntry = *(lookup.begin());

        Entry entry = firstEntry.second;

        return entry;
      }

      Entry getLastEntry() const {
        std::pair<size_t, Entry> lastEntry = *(lookup.rbegin());

        Entry entry = lastEntry.second;

        return entry;
      }

      std::string toString() const {
        std::string msg;

        return msg;
      }

    private:
      std::map<size_t, Entry> lookup;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif