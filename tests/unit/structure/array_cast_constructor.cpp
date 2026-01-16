/*
 *  array_cast_constructor.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp std::array constructor
 *  i.e. build a polynomial_nttp from a std::array
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

int main()
{
  constexpr std::array<int, 2> two_element_array{0, 1};
  constexpr auto linear = polynomial_nttp<int, 1>(two_element_array);
  static_assert(linear.coefficients[0] == 0);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr std::array<int, 3> three_element_array{-1, 0, 1};
  constexpr auto quadratic = polynomial_nttp<int, 2>(three_element_array);
  if constexpr (quadratic.degree == 2 && quadratic.coefficients.size() == 3 && quadratic.coefficients[1] == 0 &&
                quadratic.coefficients[0] == -quadratic.coefficients[2] &&
                std::same_as<std::array<int, 3>::value_type, polynomial_nttp<int, 2>::coefficient_t> &&
                std::same_as<std::array<int, 2>::size_type, polynomial_nttp<int, 1>::index_t>)
  {
    return 0; // pass
  }
  else
  {
    return 1; // fail
  }
}
