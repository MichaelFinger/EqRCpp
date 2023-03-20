#include "chapter2.hpp"
#include "chapter3.hpp"
#include "chapter4.hpp"

// #include <equating_recipes/kernel_equating.hpp>

int main(int argc, char const *argv[]) {
  EquatingRecipes::Tests::Examples::Chapter2 ch2;
  ch2();

  EquatingRecipes::Tests::Examples::Chapter3 ch3;
  ch3();

  EquatingRecipes::Tests::Examples::Chapter4 ch4;
  ch4();

  return 0;
}
