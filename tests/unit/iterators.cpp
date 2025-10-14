/*
 *  iterators.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp iterators
 *  I wonder if I should even be using iterators at this point
 *  and if I'm being fair, this test could be better... (you'll see)
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  constexpr polynomial_nttp<int, 1> linear{};
  static_assert(linear.coefficients[0] == 0 && linear.coefficients[1] == 0);
  static_assert(linear.coefficients.size() == 2);
  // even if asserts are off this ought to still work as a test
  constexpr int result_should_be_zero = [&]() {
    int result = 2;
    for (auto&& coeff : linear)
      --result;
    return result;
  }(); // range-for better work
  constexpr auto quadratic = [&]() {
    std::array<int, 3> quadratic_poly{};
    for (auto&& coeff : quadratic_poly)
      coeff = 1;
    return polynomial_nttp<int, 2>(quadratic_poly);
  }();
  if constexpr (result_should_be_zero == 0
             && quadratic.degree == 2
             && quadratic.coefficients.size() == 3
             && quadratic.coefficients[1] == 1
             && quadratic.coefficients[0] == quadratic.coefficients[2]
             && std::same_as<decltype(std::array<int, 3>{}.begin()),
                             decltype(polynomial_nttp<int, 2>{}.begin())>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
