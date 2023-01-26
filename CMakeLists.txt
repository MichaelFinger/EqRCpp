cmake_minimum_required(VERSION 3.25.0)

project(EQR_CPP
        VERSION 0.0.1
        DESCRIPTION "C++ adaptation of Equating Recipes"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS OFF)

include(CheckCXXCompilerFlag)
# include(GNUInstallDirs)
        
# set(EQR_CPP_ARCHIVE_INSTALL_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR} CACHE PATH "" FORCE)
# set(EQR_CPP_LIBRARY_INSTALL_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR} CACHE PATH "" FORCE)
# set(EQR_CPP_RUNTIME_INSTALL_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR} CACHE PATH "" FORCE)
# set(EQR_CPP_INCLUDE_INSTALL_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "" FORCE)
# if(WIN32 AND NOT CYGWIN)
#   set(DEF_INSTALL_CMAKEDIR CMake)
# else()
#   set(DEF_INSTALL_CMAKEDIR share/cmake/${PROJECT_NAME})
# endif()
# set(EQR_CPP_INSTALL_CMAKEDIR ${DEF_INSTALL_CMAKEDIR} CACHE PATH "Installation directory for CMake files" FORCE)

# if (CMAKE_GENERATOR MATCHES "Visual Studio")
#   # If Microsoft SDK is installed create script run-msbuild.bat that
#   # calls SetEnv.cmd to set up build environment and runs msbuild.
#   # It is useful when building Visual Studio projects with the SDK
#   # toolchain rather than Visual Studio.
#   include(FindSetEnv)
  
#   if (WINSDK_SETENV)
#     set(MSBUILD_SETUP "call \"${WINSDK_SETENV}\"")
#   endif()
#   # Set FrameworkPathOverride to get rid of MSB3644 warnings.
#   set(netfxpath
#        "C:\\Program Files\\Reference Assemblies\\Microsoft\\Framework\\.NETFramework\\v4.0")
#   file(WRITE run-msbuild.bat "
#     ${MSBUILD_SETUP}
#     ${CMAKE_MAKE_PROGRAM} -p:FrameworkPathOverride=\"${netfxpath}\" %*")
# endif()

include(FetchContent)

# Boost library
find_package(Boost REQUIRED COMPONENTS date_time)
include_directories(${Boost_INCLUDE_DIRS})
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
# include_directories(${Eigen_SOURCE_DIR})
# Eigen is header-only library

# Formatting library
FetchContent_Declare(
  fmtlib
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 8.0.1
  DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_MakeAvailable(fmtlib)
# include_directories(${fmtlib_SOURCE_DIR}/include)
# Adds fmt::fmt

file(GLOB_RECURSE HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

add_library(EqRCpp SHARED ${SOURCE_LIST} ${HEADER_LIST})
target_compile_features(EqRCpp PRIVATE cxx_std_20)

target_include_directories(EqRCpp PUBLIC 
                           ${CMAKE_CURRENT_SOURCE_DIR}/include
                           ${Boost_INCLUDE_DIRS}
                           ${eigen_SOURCE_DIR}
                           ${fmtlib_SOURCE_DIR}/include)

target_link_libraries(EqRCpp PUBLIC
                      ${Boost_LIBRARIES}
                      fmt::fmt
                      pthread
                      dl)