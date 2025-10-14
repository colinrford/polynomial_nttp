/*
 *  make_monomial.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp for make_monomial
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  constexpr auto linear = make_monomial<int, 1>();
  constexpr auto ten_degree_poly = make_monomial<double, 10>();
  static_assert(norm(linear) == 1 && norm(ten_degree_poly) == 10);
  static_assert(leading(linear) == 1
             && leading(linear) == leading(ten_degree_poly));
  // even if asserts are off this ought to still work as a test
  constexpr auto cubic = make_monomial<int, 3>();
  constexpr auto linear_plus_cubic = linear + cubic;
  if constexpr (norm(linear) == 1
             && norm(ten_degree_poly) == 10
             && leading(cubic) == leading(linear_plus_cubic)
             && std::same_as<polynomial_nttp<int, 1>,
                             std::remove_cvref_t<decltype(linear)>>
             && std::same_as<polynomial_nttp<double, 10>,
                             std::remove_cvref_t<decltype(ten_degree_poly)>>
             && std::same_as<polynomial_nttp<int, 3>,
                             std::remove_cvref_t<decltype(linear_plus_cubic)>>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
