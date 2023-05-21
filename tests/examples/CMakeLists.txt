add_executable(EqRCppExamples ${CMAKE_CURRENT_SOURCE_DIR}/runner.cpp)
target_compile_features(EqRCppExamples PRIVATE cxx_std_20)
target_include_directories(EqRCppExamples PRIVATE 
                           ${PROJECT_SOURCE_DIR}/include
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           ${Boost_INCLUDE_DIRS}
                           ${eigen_SOURCE_DIR}
                           ${fmtlib_SOURCE_DIR}/include
                           ${json_SOURCE_DIR}/single_include)
target_link_libraries(EqRCppExamples PRIVATE
                      ${Boost_LIBRARIES}
                      fmt::fmt
                      nlohmann_json::nlohmann_json
                      nlopt)