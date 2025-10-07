/*
 *  antiderivative.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp taking the formal antiderivative of a
 *  polynomial
 */

import std;
import polynomial_nttp;

using namespace math_nttp::polynomial;

int main()
{
  constexpr auto comparitore = [](double a, double b) {
    constexpr auto abs_val = [](double val) { return val > 0. ? val : -val; };
    if (abs_val(a - b) < 1e-7)
      return true;
    else
      return false;
  };
  constexpr polynomial_nttp<double, 1> linear{-1., 2.};
  constexpr auto quadratic = antiderivative(linear);
  static_assert(comparitore(quadratic[0], 0.) && norm(quadratic) == 2);
  static_assert(quadratic.coefficients.size() == 3);
  // even if asserts are off this ought to still work as a test
  constexpr bool quadratic_coefficients_are_correct = [&]() {
    if (!comparitore(quadratic[0], 0.))
      return false;
    else if (!comparitore(quadratic[1], -1.))
      return false;
    else if (!comparitore(quadratic[2], 1.))
      return false;
    else
      return true;
  }();

  constexpr polynomial_nttp<double, 3> cubic{-2., 7., -5., 3.};
  constexpr auto quartic = antiderivative(cubic);
  constexpr bool quartic_coefficients_are_correct = [&]() {
    if (!comparitore(quartic[0], 0.))
      return false;
    else if (!comparitore(quartic[1], -2.))
      return false;
    else if (!comparitore(quartic[2], 7. / 2.))
      return false;
    else if (!comparitore(quartic[3], -5. / 3.))
      return false;
    else if (!comparitore(quartic[4], 3. / 4.))
      return false;
    else
      return true;
  }();

  constexpr auto another_cubic = antiderivative(linear * linear);
  constexpr bool cubic_coefficients_are_correct = [&]() {
    if (!comparitore(another_cubic[0], 0.))
      return false;
    else if (!comparitore(another_cubic[1], 1.))
      return false;
    else if (!comparitore(another_cubic[2], -2.))
      return false;
    else if (!comparitore(another_cubic[3], 4. / 3.))
      return false;
    else
      return true;
  }();

  if constexpr (quadratic_coefficients_are_correct
             && quartic_coefficients_are_correct
             && cubic_coefficients_are_correct
             && std::same_as<std::remove_cvref_t<decltype(quadratic)>,
                             polynomial_nttp<double, 2>>
             && std::same_as<std::remove_cvref_t<decltype(quartic)>,
                             polynomial_nttp<double, 4>>
             && std::same_as<decltype(cubic), decltype(another_cubic)>)
  {
    return 0; // pass
  } else
  {
    return 1; // fail
  }
}
