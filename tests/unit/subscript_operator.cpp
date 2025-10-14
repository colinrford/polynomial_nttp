/*
 *  subscript_operator.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp subscript operator[]
 *  this operator is NOT an assignment operator / usable for assignment
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{};
  static_assert(linear[0] == 0 && linear[1] == 0 && linear[1] == linear[2]);
  static_assert(linear.coefficients[0] == 0 && linear.coefficients[1] == 0);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-1, 1, -3, 2};
  if constexpr (cubic.degree == 3
             && cubic.coefficients.size() == 4
             && cubic[1] == 1
             && cubic[0] == (-cubic[3] + 1)
             && cubic[2] == -3)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
