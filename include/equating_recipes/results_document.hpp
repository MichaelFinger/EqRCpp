#ifndef RESULTS_DOCUMENT_HPP
#define RESULTS_DOCUMENT_HPP

#include <map>
#include <set>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fmt/core.h>
#include <nlohmann/json.hpp>

namespace EquatingRecipes {
  class ResultsDocument {
  public:
    void setTitle(const std::string& title,
                  const std::string& subtitle = "") {
      doc["title"] = title;
      doc["subtitle"] = subtitle;
    }

    void setDescription(const std::string& description) {
      doc["decription"] = description;
    }

    void addSection(const std::string& sectionName,
                    bool overrideIfExists = false) {
      bool sectionExists = hasSection(sectionName);

      if (!sectionExists || (sectionExists && overrideIfExists)) {
        doc[sectionName] = sectionName;
      }
    }

    template<typename T>
    void addNamedValue(const std::string& name,
                       const T& value,
                       const std::string& sectionName = "") {
      if (sectionName != "") {
        doc.at(sectionName)[name] = value;
      } else {
        doc[name] = value;
      }
    }

    template<typename T>
    void addNamedArray(const std::string& name,
                       const std::vector<T>& value,
                       const std::string& sectionName = "") {
      if (sectionName != "") {
        doc.at(sectionName)[name] = value;
      } else {
        doc[name] = value;
      }
    }

    template<typename Derived>
    void addNamedEigen(const std::string& name,
                       const Eigen::EigenBase<Derived>& value,
                       const std::vector<std::string>& rowHeaders,
                       const std::vector<std::string>& columnHeaders,
                       const std::string& sectionName = "") {
      nlohmann::json eigenObject = nlohmann::json::object();

      for (size_t rowIndex = 0; rowIndex < value.rows(); rowIndex++) {
        nlohmann::json rowObject = nlohmann::json::object();

        for (size_t columnIndex = 0; columnIndex < value.rows(); columnIndex++) {
          rowObject[columnHeaders[columnIndex]] = value(rowIndex, columnIndex);
        }

        std::string rowHeader = rowHeaders.empty() ? fmt::Format("Row {}", rowIndex + 1) : rowHeaders[rowIndex];
        eigenObject[rowHeader] = rowObject;
      }

      if (sectionName != "") {
        doc.at(sectionName)[name] = eigenObject;
      } else {
        doc[name] = eigenObject;
      }
    }

  private:
    nlohmann::json doc;

    bool hasSection(const std::string& sectionName) {
      return doc.contains(sectionName);
    }
  };
} // namespace EquatingRecipes

#endif