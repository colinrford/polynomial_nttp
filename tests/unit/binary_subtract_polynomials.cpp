/*
 *  binary_subtract_polynomials.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp subtracting one polynomial from another
 *  via operator-(polynomial_nttp<R, M>, polynomial_nttp<R, N>)
 */

import std;
import polynomial_nttp;

using namespace math_nttp::polynomial;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr auto zero_linear = linear - linear;
  static_assert(zero_linear[0] == 0 && zero_linear[0] == zero_linear[1]);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto cubic_minus_linear = cubic - linear;
  constexpr auto linear_minus_cubic = linear - cubic;
  if constexpr (cubic_minus_linear[0] == -3
             && linear_minus_cubic(-1) == 15
             && cubic_minus_linear(1) == -1
             && std::same_as<decltype(linear), decltype(zero_linear)>
             && std::same_as<decltype(cubic_minus_linear),
                             decltype(linear_minus_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
