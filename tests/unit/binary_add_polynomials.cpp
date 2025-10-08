/*
 *  binary_add_polynomials.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp adding two polynomials via
 *  operator+(polynomial_nttp<R, M>, polynomial_nttp<R, N>)
 */

import std;
import polynomial_nttp;

using namespace math_nttp::polynomial;

int main() // not finished!!!!!!!
{
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr auto twice_linear = linear + linear;
  static_assert(twice_linear[0] == 2 && twice_linear[1] == 6);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto zero_cubic = cubic + (-cubic);
  constexpr auto cubic_plus_linear = cubic + linear;
  if constexpr (zero_cubic(0) == 0
             && (cubic + cubic)(1) == 6
             && cubic_plus_linear[2] == -5
             && std::same_as<decltype(linear), decltype(twice_linear)>
             && std::same_as<decltype(cubic_plus_linear), decltype(zero_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
