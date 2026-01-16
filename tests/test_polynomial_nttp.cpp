/*
 *  test_poly_nttp_mod.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  import polynomial_nttp and test all the functions found in
 *  polynomial_nttp.cppm
 *  mostly an artifact of early development at this point, the unit tests
 *  found in /tests/unit are more comprehensive than what is found here
 */

import std;
import lam.concepts;
import lam.polynomial_nttp;

using namespace lam;

consteval double line(double x)
{
  constexpr std::array<double, 2> coefficients{1., 2.};
  constexpr lam::polynomial_nttp<double, 1> p(coefficients);
  return p(x);
}

consteval double parabola(double x)
{
  constexpr std::array<double, 3> coefficients_p{0., 0., 3.};
  constexpr lam::polynomial_nttp<double, 2> p(coefficients_p);
  constexpr std::array<double, 2> coefficients_q{0., 1.};
  constexpr lam::polynomial_nttp<double, 1> q(coefficients_q);
  constexpr auto neg_q = -q;
  return (p + neg_q)(x);
}

consteval double parabola2(double x)
{
  constexpr std::array<double, 3> coefficients_p{1., 2., 3.};
  constexpr lam::polynomial_nttp<double, 2> p(coefficients_p);
  constexpr std::array<double, 2> coefficients_q{-1., 1.};
  constexpr lam::polynomial_nttp<double, 1> q(coefficients_q);
  return (q - p)(x);
}

consteval double cubic(double x)
{
  constexpr std::array<double, 3> coefficients_p{-1., 0., 2.};
  constexpr lam::polynomial_nttp<double, 2> p(coefficients_p);
  constexpr std::array<double, 2> coefficients_q{0., 3.};
  constexpr lam::polynomial_nttp<double, 1> q(coefficients_q);
  return (p * q)(x);
}

consteval auto degree_of_cubic()
{
  constexpr lam::polynomial_nttp<double, 3> p{0, 0, 0, 1};
  return norm(p);
}

consteval auto division_test()
{
  constexpr std::array<double, 9> coefficients_a{-2, 3, 0, -2, 1, -1, 2, 5, 1};
  constexpr lam::polynomial_nttp<double, 8> a(coefficients_a);
  constexpr std::array<double, 5> coefficients_b{1, 0, 2, 0, 1};
  constexpr lam::polynomial_nttp<double, 4> b(coefficients_b);
  constexpr auto a_divided_by_b = lam::division_prototype<double, 8, a, 4, b>();
  return a_divided_by_b;
}

consteval auto derivative_test()
{
  constexpr std::array<double, 3> coefficients_p{-1., 0., 1.};
  constexpr lam::polynomial_nttp<double, 2> p(coefficients_p);
  constexpr auto d_dx_p = lam::derivative(p);
  return d_dx_p[1];
}

consteval auto antiderivative_test()
{
  constexpr std::array<double, 3> coefficients_p{-1., 2., 3.};
  constexpr lam::polynomial_nttp<double, 2> p(coefficients_p);
  constexpr auto integral_of_p = lam::antiderivative(p);
  return integral_of_p[3];
}

int main()
{
  constexpr auto y_int = line(0.0);
  constexpr auto left = line(-1);
  constexpr auto right = line(1);
  constexpr auto a_zero = line(0);
  std::println("line(-1) = {}; line(1) = {}; line[0] = {}.", left, right, a_zero);
  std::println("the y intercept of the line y = 2 x + 1 equals {}", y_int);
  constexpr auto val = 1.;
  constexpr auto a = parabola(val);
  std::println("3x^2 - x evaluated at x = {} is {}", val, a);
  constexpr auto b = parabola2(val);
  std::println("-3x^2 - x - 2 evaluated at x = {} is {}", val, b);
  constexpr auto notta = cubic(a);
  constexpr lam::polynomial_nttp<double, 2> p{-1., 0., 2.};
  constexpr lam::polynomial_nttp<double, 1> q{0., 3.};
  std::println("2x^2 - 1 times 3x at x = {} is {} (deg = {})", a, notta, lam::norm(p * q));
  constexpr auto degree = degree_of_cubic();
  constexpr lam::polynomial_nttp<double, 3> x_cubed{0, 0, 0, 1};
  std::println("degree of {} x^{} is {}", x_cubed[3], x_cubed.degree, degree);
  constexpr auto a_divided_by_b = division_test();
  constexpr auto quotient = a_divided_by_b.first;
  constexpr auto remainder = a_divided_by_b.second;
  std::print("q(x) = ");
  for (std::size_t i = 0; i <= lam::norm(quotient); ++i)
  {
    if (i != lam::norm(quotient))
      std::print("{} x^{} + ", quotient[i], i);
    else
      std::print("{} x^{}", quotient[i], i);
  }
  std::print(", deg = {}\n", lam::norm(quotient));
  std::print("r(x) = ");
  for (std::size_t i = 0; i <= lam::norm(remainder); ++i)
  {
    if (i != lam::norm(remainder))
      std::print("{} x^{} + ", remainder[i], i);
    else
      std::print("{} x^{}", remainder[i], i);
  }
  std::print(", deg = {}\n", lam::norm(remainder));
  static_assert(lam::norm(lam::derivative(lam::polynomial_nttp<double, 0>{})) == 0);
  constexpr auto slope_of_x_squared = derivative_test();
  std::println("slope of x^2 - 1 is {}x", slope_of_x_squared);
  constexpr auto coefficient_of_ad = antiderivative_test();
  std::println("4th coefficient of antiderivative of 3x^2 + 2x - 1 is {}", coefficient_of_ad);

  static_assert(lam::concepts::experimental::field_element_c_weak<double>);
}
