#ifndef JSON_DOCUMENT_HPP
#define JSON_DOCUMENT_HPP

#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

namespace EquatingRecipes {
  namespace JSON {
    struct JsonDocument {
      std::string dump() {
        return this->json.dump(2);
      }

      void setJson(const nlohmann::json& json) {
        this->json = json;
      }

      void toTextFile(const std::string& filename) {
        std::ofstream ofs;

        ofs.open(filename);

        ofs << json.dump(2);

        ofs.close();
      }

      void fromTextFile(const std::string& filename) {
        std::ifstream ifs(filename);

        this->json = nlohmann::json::parse(ifs);

        ifs.close();
      }

      nlohmann::json json;      
    };
  } // namespace JSON
} // namespace EquatingRecipes

#endif