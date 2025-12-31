/*
 *  norm.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp norm of a polynomial
 */

import lam.polynomial_nttp;

using namespace lam;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{1, 3};
  constexpr polynomial_nttp<double, 10> ten_degree_poly{};
  static_assert(norm(linear) == 1 && norm(ten_degree_poly) == 10);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, norm(linear) + 2> cubic{-1, 1, -1, 1};
  if constexpr (norm(linear) == 1
             && norm(ten_degree_poly) == 10
             && norm(cubic) == 3
             && norm(polynomial_nttp<float, 5>{}) == 5)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
