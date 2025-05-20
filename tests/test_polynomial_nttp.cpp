#include <polynomial_nttp.cpp>
#include <print>

consteval double line(const double x)
{
  constexpr std::array<double, 2> coefficients{1.,2.};
  constexpr math::polynomial_nttp<double, 1> p(std::move(coefficients));
  return p(x);
}

consteval double parabola(const double x)
{
  constexpr std::array<double, 3> coefficients_p{0.,0.,3.};
  constexpr math::polynomial_nttp<double, 2> p(std::move(coefficients_p));
  constexpr std::array<double, 2> coefficients_q{0.,1.};
  constexpr math::polynomial_nttp<double, 1> q(std::move(coefficients_q));
  constexpr auto neg_q = -q;
  return (p + neg_q)(x);
}

consteval double parabola2(const double x)
{
  constexpr std::array<double, 3> coefficients_p{1.,2.,3.};
  constexpr math::polynomial_nttp<double, 2> p(std::move(coefficients_p));
  constexpr std::array<double, 2> coefficients_q{-1.,1.};
  constexpr math::polynomial_nttp<double, 1> q(std::move(coefficients_q));
  return (q - p)(x);
}

consteval double cubic(const double x)
{
  constexpr std::array<double, 3> coefficients_p{-1.,0.,2.};
  constexpr math::polynomial_nttp<double, 2> p(std::move(coefficients_p));
  constexpr std::array<double, 2> coefficients_q{0.,3.};
  constexpr math::polynomial_nttp<double, 1> q(std::move(coefficients_q));
  return (p * q)(x);
}

consteval auto degree_of_cubic()
{
  constexpr math::polynomial_nttp<double, 3> p{0, 0, 0, 1};
  return norm(p);
}

consteval auto division_test()
{
  constexpr std::array<double, 9> coefficients_a{-2, 3, 0, -2, 1, -1, 2, 5, 1};
  constexpr math::polynomial_nttp<double, 8> a(std::move(coefficients_a));
  constexpr std::array<double, 5> coefficients_b{1, 0, 2, 0, 1};
  constexpr math::polynomial_nttp<double, 4> b(std::move(coefficients_b));
  constexpr auto a_divided_by_b = math::
                                  division_prototype<double, 8, a, 4, b>();
  return a_divided_by_b;
}

consteval auto derivative_test()
{
  constexpr std::array<double, 3> coefficients_p{-1.,0.,1.};
  constexpr math::polynomial_nttp<double, 2> p(std::move(coefficients_p));
  constexpr auto d_dx_p = math::derivative(p);
  return d_dx_p[1];
}

consteval auto antiderivative_test()
{
  constexpr std::array<double, 3> coefficients_p{-1.,2.,3.};
  constexpr math::polynomial_nttp<double, 2> p(std::move(coefficients_p));
  constexpr auto integral_of_p = math::antiderivative(p);
  return integral_of_p[3];
}

int main()
{
  constexpr auto y_int = line(0.0);
  constexpr auto left = line(-1);
  constexpr auto right = line(1);
  constexpr auto a_zero = line(0);
  std::println("line(-1) = {}; line(1) = {}; line[0] = {}.", left,
                                                             right,
                                                             a_zero);
  std::println("the y intercept of the line y = 2 x + 1 equals {}", y_int);
  constexpr auto val = 1.;
  constexpr auto a = parabola(val);
  std::println("3x^2 - x evaluated at x = {} is {}", val, a);
  constexpr auto b = parabola2(val);
  std::println("-3x^2 - x - 2 evaluated at x = {} is {}", val, b);
  constexpr auto notta = cubic(a);
  constexpr math::polynomial_nttp<double, 2> p{-1.,0.,2.};
  constexpr math::polynomial_nttp<double, 1> q{0.,3.};
  std::println("2x^2 - 1 times 3x at x = {} is {} (deg = {})",
               a,
               notta,
               math::norm(p * q));
  constexpr auto degree = degree_of_cubic();
  constexpr math::polynomial_nttp<double, 3> x_cubed{0, 0, 0, 1};
  std::println("degree of {} x^{} is {}", x_cubed[3], x_cubed.degree, degree);
  constexpr auto a_divided_by_b = division_test();
  constexpr auto quotient = a_divided_by_b.first;
  constexpr auto remainder = a_divided_by_b.second;
  std::print("q(x) = ");
  for (std::size_t i = 0; i <= math::norm(quotient); ++i)
  {
    if (i != math::norm(quotient))
      std::print("{} x^{} + ", quotient[i], i);
    else
      std::print("{} x^{}", quotient[i], i);
  }
  std::print(", deg = {}\n", math::norm(quotient));
  std::print("r(x) = ");
  for (std::size_t i = 0; i <= math::norm(remainder); ++i)
  {
    if (i != math::norm(remainder))
      std::print("{} x^{} + ", remainder[i], i);
    else
      std::print("{} x^{}", remainder[i], i);
  }
  std::print(", deg = {}\n", math::norm(remainder));
  constexpr auto slope_of_x_squared = derivative_test();
  std::println("slope of x^2 - 1 is {}x", slope_of_x_squared);
  constexpr auto coefficient_of_ad = antiderivative_test();
  std::println("4th coefficient of antiderivative of 3x^2 + 2x - 1 is {}",
                coefficient_of_ad);
}
