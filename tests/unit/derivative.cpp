/*
 *  derivative.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp taking the formal derivative of a polynomial
 */

import std;
import polynomial_nttp;

using namespace math_nttp;

int main()
{
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  constexpr polynomial_nttp<int, 1> linear{-1, 3};
  constexpr auto constant = derivative(linear);
  static_assert(constant[0] == 3 && norm(constant) == 0);
  static_assert(constant.coefficients.size() == 1);
  // even if asserts are off this ought to still work as a test
  constexpr polynomial_nttp<int, 3> cubic{-2, 7, -5, 3};
  constexpr auto quadratic = derivative(cubic);
  constexpr bool quadratic_coefficients_are_correct = [&]() {
    if (quadratic[0] != 7)
      return false;
    else if (quadratic[1] != -10)
      return false;
    else if (quadratic[2] != 9)
      return false;
    else
      return true;
  }();

  constexpr auto another_cubic = derivative(linear * cubic);
  constexpr bool cubic_coefficients_are_correct = [&]() {
    if (another_cubic[0] != -13)
      return false;
    else if (another_cubic[1] != 52)
      return false;
    else if (another_cubic[2] != -54)
      return false;
    else if (another_cubic[3] != 36)
      return false;
    else
      return true;
  }();

  if constexpr (leading(derivative(constant)) == 0
             && norm(derivative(constant)) == 0
             && quadratic_coefficients_are_correct
             && quadratic(-1) == 26
             && cubic_coefficients_are_correct
             && std::same_as<std::remove_cvref_t<decltype(quadratic)>,
                             polynomial_nttp<int, 2>>
             && std::same_as<decltype(cubic), decltype(another_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
}
