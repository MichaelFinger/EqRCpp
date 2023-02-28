cmake_minimum_required(VERSION 3.25.0)

project(ER
        VERSION 1.0.0
        DESCRIPTION "Equating Recipes"
        LANGUAGES C)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(EquatingRecipes SHARED ${CMAKE_CURRENT_SOURCE_DIR}/RGandSG_NoSmooth.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/Analytic_SEs.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/BetaBinomial.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/Bootstrap.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/CG_EquiEquate.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/CG_NoSmooth.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/CLL_Equate.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/ERutilities.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/IRTst.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/kernel_Equate.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/LogLinear.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/matrix.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/NRutilities.c
                                   ${CMAKE_CURRENT_SOURCE_DIR}/IRTeq.c)



target_include_directories(EquatingRecipes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Source_and_Headers)

target_compile_options(EquatingRecipes PRIVATE "-Wno-return-type;-Wno-visibility;-Wno-comment;-Wno-incompatible-pointer-types;-Wno-parentheses;-Wno-absolute-value;-Wno-unused-value;-Wno-format;-Wno-unsequenced;-Wno-format-extra-args")
