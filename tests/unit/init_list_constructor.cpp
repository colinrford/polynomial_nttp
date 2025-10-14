/*
 *  init_list_constructor.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp initializer_list constructor
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{0, 1};
  static_assert(linear.coefficients[0] == 0);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 2> quadratic{1, 0, 1};
  if constexpr (quadratic.degree == 2
             && quadratic.coefficients.size() == 3
             && quadratic.coefficients[1] == 0
             && quadratic.coefficients[0] == quadratic.coefficients[2])
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
