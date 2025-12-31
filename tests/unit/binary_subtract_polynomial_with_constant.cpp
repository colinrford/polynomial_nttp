/*
 *  binary_subtract_polynomial_with_constant.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp subtracting a polynomial with a coefficient_t
 *  via operator-(R, polynomial_nttp<R, N>) and its cyclic permutation
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{0, 3};
  constexpr auto linear_minus_constant = linear - 1;
  static_assert(linear_minus_constant[0] == -1
             && linear_minus_constant[1] == 3);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{3, -2, -5, 7};
  constexpr auto constant_minus_cubic = 2 - cubic;
  if constexpr (constant_minus_cubic[0] == -1
             && constant_minus_cubic[0] == constant_minus_cubic(0)
             && constant_minus_cubic(1) == -1
             && (-2 - cubic - 1)[0] == -6
             && std::same_as<decltype(linear), decltype(linear_minus_constant)>
             && std::same_as<decltype(cubic), decltype(constant_minus_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
