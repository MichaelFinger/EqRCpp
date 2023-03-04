cmake_minimum_required(VERSION 3.25.0)

project(EQR_CPP
        VERSION 0.0.1
        DESCRIPTION "C++ adaptation of Equating Recipes"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckCXXCompilerFlag)

include(FetchContent)

# Boost library
find_package(Boost REQUIRED)
# Adds Boost::boost

# Eigen math library
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
endif()
# Eigen is header-only library

# Formatting library
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 9.1.0
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_MakeAvailable(fmtlib)
# Adds fmt::fmt

# json
FetchContent_Declare(json 
  URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_MakeAvailable(json)
# adds nlohmann_json::nlohmann_json

# nlopt
set(NLOPT_CXX ON CACHE BOOL "enable cxx routines" FORCE)
set(NLOPT_FORTRAN OFF CACHE BOOL "enable fortran tests" FORCE)
set(NLOPT_PYTHON OFF CACHE BOOL "build python bindings" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "build octave bindings" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "build matlab bindings" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "build guile bindings" FORCE)
set(NLOPT_SWIG OFF CACHE BOOL "use SWIG to build bindings" FORCE)

FetchContent_Declare(
  nlopt
  GIT_REPOSITORY https://github.com/stevengj/nlopt.git
  GIT_TAG        v2.7.1
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_MakeAvailable(nlopt)
# adds nlopt

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Eigen_SOURCE_DIR})
include_directories(${fmtlib_SOURCE_DIR}/include)
include_directories(${json_SOURCE_DIR}/single_include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# add_executable(EqRCppExamples ${CMAKE_CURRENT_SOURCE_DIR}/tests/examples/runner.cpp)
# target_compile_features(EqRCppExamples PRIVATE cxx_std_20)
# target_include_directories(EqRCppExamples PRIVATE 
#                            ${CMAKE_CURRENT_SOURCE_DIR}
#                            ${CMAKE_CURRENT_SOURCE_DIR}/include
#                            ${Boost_INCLUDE_DIRS}
#                            ${eigen_SOURCE_DIR}
#                            ${fmtlib_SOURCE_DIR}/include
#                            ${json_SOURCE_DIR}/single_include)
# target_link_libraries(EqRCppExamples PRIVATE
#                       ${Boost_LIBRARIES}
#                       fmt::fmt
#                       nlohmann_json::nlohmann_json
#                       nlopt
#                       pthread
#                       dl)

add_subdirectory("tests/examples")