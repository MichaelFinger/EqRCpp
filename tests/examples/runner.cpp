#include "chapter2.hpp"
#include "chapter3.hpp"
#include "chapter4.hpp"
#include "chapter5.hpp"
#include "chapter6.hpp"

int main(int argc, char const *argv[]) {
  EquatingRecipes::Tests::Examples::Chapter2 ch2;
  ch2();

  EquatingRecipes::Tests::Examples::Chapter3 ch3;
  ch3();

  EquatingRecipes::Tests::Examples::Chapter4 ch4;
  ch4();

  EquatingRecipes::Tests::Examples::Chapter5 ch5;
  ch5();

  EquatingRecipes::Tests::Examples::Chapter6 ch6;
  ch6();

  return 0;
}
