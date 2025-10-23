/*
 *  binary_subtract_polynomials.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp multiplying a polynomial with a coefficient_t
 *  via operator*(R, polynomial_nttp<R, N>) and its cyclic permutation
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{-2, 1};
  constexpr auto thrice_linear = 3 * linear;
  static_assert(thrice_linear[0] == -6 && thrice_linear[1] == 3);
  static_assert(thrice_linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto constant_times_cubic = -1 * cubic;
  if constexpr (constant_times_cubic[0] == 2
             && constant_times_cubic(-1) == 17
             && (2 * cubic * -1)(1) == -6
             && std::same_as<decltype(linear), decltype(thrice_linear)>
             && std::same_as<decltype(cubic), decltype(constant_times_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
