#ifndef EQUATING_RECIPES_INPUT_DATASET_HPP
#define EQUATING_RECIPES_INPUT_DATASET_HPP

#include <fstream>
#include <ranges>
#include <string>
#include <boost/algorithm/string.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Fixtures {
      class EquatingRecipeInputDataset {
      public:
        void import(const std::string& foldername) {
          this->foldername = foldername;
          configure();
          std::ifstream ifs = open();
          readLines(ifs);
          close(ifs);
        }

      protected:
        std::string foldername;
        std::string filename;

        virtual void configure() = 0;

        void close(std::ifstream& ifs) {
          if (ifs.is_open()) {
            ifs.close();
          }
        }

        std::ifstream open() {
          std::ifstream ifs;

          std::string fullFilename = foldername + filename;

          ifs.open(fullFilename);

          return ifs;
        }

        void readLines(std::ifstream& ifs) {
          if (ifs.is_open()) {
            std::string lineRead;

            while (std::getline(ifs, lineRead)) {
              processLine(lineRead);
            }
          }
        }

        size_t getNumberOfLines(size_t startRow = 0) {
          std::string lineRead;
          size_t numberOfLines = 0;

          std::ifstream ifs = open();

          if (ifs.is_open()) {
            while (std::getline(ifs, lineRead)) {
              if (!lineRead.empty()) {
                numberOfLines++;
              }
            }
          }

          close(ifs);

          return numberOfLines;
        }

        bool lineIsEmpty(const std::string& line) {
          std::string testValue = line;

          boost::trim(testValue);
          
          bool isEmpty = testValue.empty();

          return isEmpty;
        }

        /*
        string str1("hello abc-*-ABC-*-aBc goodbye");

    typedef vector< iterator_range<string::iterator> > find_vector_type;
    
    find_vector_type FindVec; // #1: Search for separators
    ifind_all( FindVec, str1, "abc" ); // FindVec == { [abc],[ABC],[aBc] }

    typedef vector< string > split_vector_type;
    
    split_vector_type SplitVec; // #2: Search for tokens
    split( SplitVec, str1, is_any_of("-*"), token_compress_on ); // SplitVec == { "hello abc","ABC","aBc goodbye" }
        */
        std::vector<std::string> splitLine(const std::string& line, const std::string& delimiter = ",") {
          std::vector<std::string> values;
          std::string str = line;
          std::string delim = delimiter;

          boost::split(values, str, boost::is_any_of(delim), boost::token_compress_on );

          return values;
        }

        virtual void processLine(const std::string& lineRead) = 0;
      };
    } // namespace Fixtures
  }   // namespace Tests
} // namespace EquatingRecipes

#endif