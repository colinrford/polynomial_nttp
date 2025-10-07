/*
 *  binary_add_polynomial_with_constant.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp adding a polynomial with a coefficient_t
 *  via operator+(R, polynomial_nttp<R, N>) and its cyclic permutation
 */

import std;
import polynomial_nttp;

using namespace math_nttp::polynomial;

int main() // not finished!!!!!!!
{
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr auto linear_plus_constant = linear + 1;
  static_assert(linear_plus_constant[0] == 2 && linear_plus_constant[1] == 3);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto constant_plus_cubic = 2 + cubic;
  if constexpr (constant_plus_cubic[0] == constant_plus_cubic(0)
             && constant_plus_cubic(1) == 5
             && (-2 + cubic + 1)[0] == -3
             && std::same_as<decltype(linear), decltype(linear_plus_constant)>
             && std::same_as<decltype(cubic), decltype(constant_plus_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
