#ifndef JSON_DOCUMENT_HPP
#define JSON_DOCUMENT_HPP

#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

namespace EquatingRecipes {
  namespace JSON {
    class JsonDocument {
    public:
      void setJson(const nlohmann::json& json) {
        this->json = json;
      }

      void toTextFile(const std::string& filename) {
        std::ofstream ofs;

        ofs.open(filename);

        ofs << json.dump(2);

        ofs.close();
      }

    private:
      nlohmann::json json;
    };
  } // namespace JSON
} // namespace EquatingRecipes

#endif