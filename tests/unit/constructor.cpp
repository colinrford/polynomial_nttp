/*
 *  constructor.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp constructor
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

int main()
{
  static_assert(polynomial_nttp{}.degree == 0);
  static_assert(polynomial_nttp{}.coefficients.size() == 1);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<double, 3> cubic{};
  if constexpr (cubic.degree == 3
             && cubic.coefficients.size() == 4
             && cubic.coefficients[0] == cubic.coefficients[3]
             && std::same_as<double,
                             polynomial_nttp<double, 3>::coefficient_t>
             && std::same_as<std::size_t,
                             polynomial_nttp<int, 3>::index_t>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
