/*
 *  binary_multiply_polynomials.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp multiplying two polynomials
 *  via operator*(polynomial_nttp<R, M>, polynomial_nttp<R, N>)
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{-1, 3};
  constexpr auto linear_squared = linear * linear;
  static_assert(linear_squared[0] == 1 && linear_squared[2] == 9);
  static_assert(linear_squared.coefficients.size() == 3);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto cubic_times_linear = cubic * linear;
  constexpr auto cubic_squared = cubic * cubic;
  if constexpr (cubic_times_linear[0] == 2
             && cubic_squared(0) == 4
             && cubic_times_linear(-3) == 1490
             && std::same_as<std::remove_cvref_t<decltype(cubic_times_linear)>,
                             polynomial_nttp<int, 4>>
             && std::same_as<std::remove_cvref_t<decltype(cubic_squared)>,
                             polynomial_nttp<int, 6>>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
