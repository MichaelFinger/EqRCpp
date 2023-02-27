cmake_minimum_required(VERSION 3.25.0)

project(ER
        VERSION 1.0.0
        DESCRIPTION "Equating Recipes"
        LANGUAGES C)

set(CMAKE_C_STANDARD 99)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(EquatingRecipes SHARED NRutilities.c
                                   ERutilities.c
                                   matrix.c
                                   IRTst.c
                                   Analytic_SEs.c
                                   BetaBinomial.c
                                   Bootstrap.c
                                   CG_EquiEquate.c
                                   CG_NoSmooth.c
                                   CLL_Equate.c
                                   CubicSpline.c
                                   LogLinear.c
                                   RGandSG_NoSmooth.c
                                   kernel_Equate.c
                                   IRTeq.c)

target_include_directories(EquatingRecipes PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_options(EquatingRecipes PRIVATE "-Wno-return-type;-Wno-visibility;-Wno-comment;-Wno-incompatible-pointer-types;-Wno-parentheses;-Wno-absolute-value;-Wno-unused-value;-Wno-format;-Wno-unsequenced;-Wno-format-extra-args")



