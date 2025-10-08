/*
 *  constant_cast_constructor.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp converting a constant to an N-degree
 *  polynomial
 *  it ought to always be a 0-degree polynomial, but the syntax
 *    does allow for creating an N-degree polynomial that is
 *    identically 0, or that has 0 as coefficient N
 */

import std;
import polynomial_nttp;

using namespace math_nttp::polynomial;

constexpr bool comparitore(double a, double b)
{
  constexpr auto abs_val = [](double val) { return val > 0. ? val : -val; };
  if (abs_val(a - b) < 1e-7)
    return true;
  else
    return false;
}

int main()
{
  constexpr polynomial_nttp<double, 0> constant(1.);
  constexpr polynomial_nttp another_constant(2.);
  constexpr polynomial_nttp<double, 100> wasteful_constant(3.);
  constexpr bool wasteful_constant_is_indeed_wasteful = [&]() {
    for (auto&& coefficient : wasteful_constant)
    {
      if (!comparitore(0., coefficient)
       && !comparitore(wasteful_constant[0], coefficient))
        return false;
    }
    return true;
  }();
  if constexpr (norm(constant) == norm(another_constant)
              && wasteful_constant_is_indeed_wasteful)
  {
    return 0;
  } else
  {
    return 1;
  }
}
