/*
 *  unary_negate_operator.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp parentheses operator-(coefficient_t)
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr auto negated_linear = -linear;
  static_assert(negated_linear[0] == -1 && negated_linear[1] == -3);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto negated_cubic = -cubic;
  if constexpr (negated_cubic(0) == 2
             && negated_cubic(1) == -3
             && negated_cubic[2] == 5
             && std::same_as<decltype(linear), decltype(negated_linear)>
             && std::same_as<decltype(cubic), decltype(negated_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
