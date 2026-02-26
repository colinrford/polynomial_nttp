/*
 *  lagrange_test.cpp
 *  unit test for lagrange interpolation
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

// Helper for finite field approximation (exact equality for Zq)
template<typename T>
bool check_poly_approx_equal(const polynomial_nttp<T, 0>& p, std::initializer_list<T> coeffs)
{
  if (p.degree != 0)
    return false;
  return is_approx_equal(p.coefficients[0], *coeffs.begin());
}

template<typename T, std::size_t N>
bool check_poly_approx_equal(const polynomial_nttp<T, N>& p, std::initializer_list<T> coeffs)
{
  if (coeffs.size() != N + 1)
    return false;
  auto it = coeffs.begin();
  for (std::size_t i = 0; i <= N; ++i)
    if (!is_approx_equal(p.coefficients[i], *it++))
      return false;
  return true;
}

int main()
{
  int failures = 0;

  // Test 1: Real (double) - Linear: y = x
  {
    std::array<double, 2> x = {0.0, 1.0};
    std::array<double, 2> y = {0.0, 1.0};
    auto poly = lagrange_interpolate<double, 1>(x, y);

    // Expected: 0 + 1*x
    if (!is_approx_equal(poly.coefficients[0], 0.0) || !is_approx_equal(poly.coefficients[1], 1.0))
    {
      std::println("ACCELERATE FAIL: double linear");
      failures++;
    }
  }

  // Test 2: Real (double) - Quadratic: y = x^2
  {
    std::array<double, 3> x = {-1.0, 0.0, 1.0};
    std::array<double, 3> y = {1.0, 0.0, 1.0};
    auto poly = lagrange_interpolate<double, 2>(x, y);

    // Expected: x^2 -> [0, 0, 1]
    if (!is_approx_equal(poly.coefficients[0], 0.0) || !is_approx_equal(poly.coefficients[1], 0.0) ||
        !is_approx_equal(poly.coefficients[2], 1.0))
    {
      std::println("ACCELERATE FAIL: double quadratic");
      failures++;
    }
  }

  // Test 3: Complex (std::complex<double>) - y = i*x
  {
    using C = std::complex<double>;
    std::array<C, 2> x = {C(0, 0), C(1, 0)};
    std::array<C, 2> y = {C(0, 0), C(0, 1)}; // 0 -> 0, 1 -> i
    auto poly = lagrange_interpolate<C, 1>(x, y);

    // Expected: 0 + i*x
    if (!is_approx_equal(poly.coefficients[0], C(0, 0)) || !is_approx_equal(poly.coefficients[1], C(0, 1)))
    {
      std::println("ACCELERATE FAIL: complex linear");
      failures++;
    }
  }

  // Test 4: Single Point (Constant)
  {
    std::array<double, 1> x = {10.0};
    std::array<double, 1> y = {42.0};
    auto poly = lagrange_interpolate<double, 0>(x, y);

    // Expected: 42
    if (!is_approx_equal(poly.coefficients[0], 42.0))
    {
      std::println("ACCELERATE FAIL: single point");
      failures++;
    }
  }

  // Test 5: Real (double) - Collinear points for Quadratic (Effectively Linear)
  {
    std::array<double, 3> x = {0.0, 1.0, 2.0};
    std::array<double, 3> y = {1.0, 2.0, 3.0};
    auto poly = lagrange_interpolate<double, 2>(x, y);

    // Expected: 1 + 1*x + 0*x^2 -> [1, 1, 0]
    if (!check_poly_approx_equal(poly, {1.0, 1.0, 0.0}))
    {
      std::println("FAIL: double effectively linear");
      failures++;
    }
  }

  // Test 6: Real (double) - Cubic: y = x^3 - 2x + 1
  {
    std::array<double, 4> x = {-2.0, -1.0, 0.0, 1.0};
    std::array<double, 4> y = {-3.0, 2.0, 1.0, 0.0};
    auto poly = lagrange_interpolate<double, 3>(x, y);

    // Expected: 1 - 2*x + 0*x^2 + 1*x^3 -> [1, -2, 0, 1]
    if (!check_poly_approx_equal(poly, {1.0, -2.0, 0.0, 1.0}))
    {
      std::println("FAIL: double cubic");
      failures++;
    }
  }

  // Test 7: Compile-time evaluation (constexpr)
  {
    constexpr std::array<double, 3> x_c = {-1.0, 0.0, 1.0};
    constexpr std::array<double, 3> y_c = {1.0, 0.0, 1.0};
    constexpr auto poly_c = lagrange_interpolate<double, 2>(x_c, y_c);

    static_assert(poly_c.degree == 2, "Degree must be 2");
    static_assert(poly_c.coefficients[0] == 0.0, "Coefficient 0 must be 0");
    static_assert(poly_c.coefficients[1] == 0.0, "Coefficient 1 must be 0");
    static_assert(poly_c.coefficients[2] == 1.0, "Coefficient 2 must be 1");
  }

  // Test 8: Comprehensive Interpolation Checks
  {
    std::array<double, 4> x = {1.0, 2.0, 3.0, 4.0};
    std::array<double, 4> y = {5.0, 2.0, 6.0, 1.0};
    auto poly = lagrange_interpolate<double, 3>(x, y);

    bool evaluation_check_passed = true;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
      if (!is_approx_equal(poly(x[i]), y[i], 1e-12))
      {
        evaluation_check_passed = false;
        std::println("FAIL: double comprehensive. Expected {} at x={}, got {}", y[i], x[i], poly(x[i]));
      }
    }

    if (!evaluation_check_passed)
      failures++;
  }

  return failures;
}
