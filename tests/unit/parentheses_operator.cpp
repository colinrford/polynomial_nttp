/*
 *  parenthesis_operator.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp parentheses operator()
 *  test could be better... (TODO)
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{-1, 1};
  static_assert(linear(0) == -1 && linear(1) == 0);
  static_assert(linear(-1) == -2 && linear(2) == 1);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-1, 1, -3, 2};
  if constexpr (cubic(0) == cubic(1)
             && cubic(2) == 5
             && cubic(-2) == -31)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
