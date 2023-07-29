#ifndef TESTS_EXAMPLES_DATASETS_DUMMY_V3_ITEMS_HPP
#define TESTS_EXAMPLES_DATASETS_DUMMY_V3_ITEMS_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      namespace Datasets {
        struct DummyV3Items {
          size_t numberOfItems;
          std::vector<std::pair<size_t, size_t>> itemPairSpecs;

          void import() {
            std::string filename = "dummyV3items.txt";

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

              std::pair<size_t, size_t> itemPairSpec;
              itemPairSpec.first = static_cast<size_t>(std::stoi(values[0]));
              itemPairSpec.second = static_cast<size_t>(std::stoi(values[1]));

              itemPairSpecs.push_back(itemPairSpec);
            }
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