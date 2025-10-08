/*
 *  leading_coefficient.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp obtaining leading coefficient
 */

import polynomial_nttp;

using namespace math_nttp::polynomial;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr polynomial_nttp<double, 10> ten_degree_poly{};
  static_assert(leading(linear) == 3 && leading(ten_degree_poly) == 0);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-1, 1, -1, 1};
  if constexpr (leading(linear) == 3
             && leading(ten_degree_poly) == 0
             && leading(cubic) == 1
             && leading(linear - 2 * cubic) == -2)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
