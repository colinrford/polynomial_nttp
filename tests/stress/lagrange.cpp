/*
 *  lagrange.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Stress test for polynomial_nttp Lagrange interpolation stability.
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;



// ----------------------------------------------------------------------------
// Compile-Time Evaluation Check
// ----------------------------------------------------------------------------
consteval bool verify_lagrange_constexpr()
{
  std::array<double, 5> x = {1.0, 2.0, 3.0, 4.0, 5.0};
  std::array<double, 5> y = {5.0, 2.0, 6.0, 1.0, 3.0};
  auto poly = lagrange_interpolate<double, 4>(x, y);

  for (std::size_t i = 0; i < x.size(); ++i)
  {
    if (!is_relatively_approx_equal(y[i], poly(x[i]), 1e-12))
    {
      return false;
    }
  }
  return true;
}
static_assert(verify_lagrange_constexpr(), "Lagrange interpolation constexpr evaluation failed parity.");

int main()
{
  std::println("Stress Test: Lagrange Interpolation");
  std::println("===================================");
  int failures = 0;

  // 1. Extreme Precision Scaling Test
  // Interpolating large differences like 10^6 mixing with 10^-6
  {
    std::println("  Running Precision Extreme Scaling Test...");
    std::array<double, 4> x = {-10.0, -1.0, 1.0, 10.0};
    // Create an extreme y-scale
    std::array<double, 4> y = {1e6, -1e-6, 1e-6, -1e6};

    auto poly = lagrange_interpolate<double, 3>(x, y);

    for (std::size_t i = 0; i < x.size(); ++i)
    {
      if (!is_relatively_approx_equal(y[i], poly(x[i]), 1e-8))
      {
        std::println("    [FAIL] Extreme Scale: Expected {} at x={}, got {}", y[i], x[i], poly(x[i]));
        failures++;
      }
    }
  }

  // 2. High-Degree Runge Phenomenon Resistance Test (Chebyshev Nodes)
  // N=15
  {
    std::println("  Running High-Degree Target Test (N=15, Chebyshev Nodes)...");
    constexpr std::size_t N = 15;
    std::array<double, N> x{};
    std::array<double, N> y{};

    // Target mathematical function: y = sin(x) + e^(x/10) over interval [-1, 1]
    for (std::size_t i = 0; i < N; ++i)
    {
      // Generate Chebyshev nodes mapped to [-1, 1]
      x[i] = std::cos((2.0 * i + 1.0) / (2.0 * N) * std::numbers::pi);
      y[i] = std::sin(x[i]) + std::exp(x[i] / 10.0);
    }

    auto poly = lagrange_interpolate<double, N - 1>(x, y);

    // Verify that the polynomial structurally loops back to the exact evaluation points robustly
    for (std::size_t i = 0; i < x.size(); ++i)
    {
      if (!is_relatively_approx_equal(y[i], poly(x[i]), 1e-9))
      {
        std::println("    [FAIL] High-Degree Chebyshev: Expected {} at x={}, got {}", y[i], x[i], poly(x[i]));
        failures++;
      }
    }
  }

  if (failures == 0)
  {
    std::println("All Lagrange Stress Tests Passed.");
    return 0;
  }
  else
  {
    std::println("Lagrange Stress Tests FAILED. ({} errors)", failures);
    return 1;
  }
}
