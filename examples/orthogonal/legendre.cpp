/*
 *  legendre.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Legendre polynomial generation example
 *    Demonstrates generating Legendre polynomials P_n(x) at compile time
 *    using recursive constexpr functions.
 */

#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <print>
#include <vector>

import lam.polynomial_nttp;

#ifdef HAS_BOOST_MATH
#include <boost/math/special_functions/legendre.hpp>
#endif

namespace lam::orthogonal
{

// Generate the nth Legendre polynomial P_n(x)
// Recurrence:
//   P_0(x) = 1
//   P_1(x) = x
//   (n+1)P_{n+1}(x) = (2n+1)x P_n(x) - n P_{n-1}(x)
//   => P_{n+1}(x) = ((2n+1)/(n+1)) x P_n(x) - (n/(n+1)) P_{n-1}(x)
template<typename R, std::size_t N>
struct legendre_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0) // P_0(x) = 1
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1) // P_1(x) = x
      return lam::polynomial_nttp<R, 1>{{R(0), R(1)}};
    else
    {
      // Recursive case: calculating P_n where n = N
      // Recurrence uses P_{N-1} and P_{N-2}

      // We essentially want P_N, so let's map the recurrence indices:
      // P_N = ((2(N-1)+1)/N) * x * P_{N-1} - ((N-1)/N) * P_{N-2}

      constexpr auto P_nm1 = legendre_memo<R, N - 1>::value;
      constexpr auto P_nm2 = legendre_memo<R, N - 2>::value;

      constexpr auto x = lam::make_monomial<R, 1>();

      constexpr R two_n_minus_1 = R(2 * N - 1); // 2(N-1) + 1 = 2N - 1
      constexpr R n_minus_1 = R(N - 1);
      constexpr R n = R(N);
      constexpr R inv_n = R(1) / n;

      return (two_n_minus_1 * x * P_nm1 - n_minus_1 * P_nm2) * inv_n;
    }
  }();
};

template<typename R, std::size_t N>
constexpr auto legendre_n()
{
  return legendre_memo<R, N>::value;
}

} // namespace lam::orthogonal

// Runtime Reference Implementation
double reference_legendre(unsigned n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;

  double P_prev2 = 1.0;
  double P_prev1 = x;

  for (unsigned k = 1; k < n; ++k)
  {
    // (k+1)P_{k+1} = (2k+1)x P_k - k P_{k-1}
    double P_curr = ((2.0 * k + 1.0) * x * P_prev1 - k * P_prev2) / (k + 1.0);
    P_prev2 = P_prev1;
    P_prev1 = P_curr;
  }
  return P_prev1;
}

int main()
{
  using R = double;

  std::println("=== Legendre Polynomials P_n(x) ===\n");

  constexpr auto P2 = lam::orthogonal::legendre_n<R, 2>();
  constexpr auto P3 = lam::orthogonal::legendre_n<R, 3>();
  constexpr auto P4 = lam::orthogonal::legendre_n<R, 4>();

  std::println("P_2 = {}", P2); // 0.5(3x^2 - 1) = 1.5x^2 - 0.5
  std::println("P_3 = {}", P3); // 0.5(5x^3 - 3x) = 2.5x^3 - 1.5x
  std::println("P_4 = {}", P4); // 0.125(35x^4 - 30x^2 + 3) = 4.375x^4 - 3.75x^2 + 0.375

  // Compile-time verification
  static_assert(lam::is_approx_equal(P2.coefficients[2], R(1.5)));
  static_assert(lam::is_approx_equal(P2.coefficients[0], R(-0.5)));

  static_assert(lam::is_approx_equal(P3.coefficients[3], R(2.5)));
  static_assert(lam::is_approx_equal(P3.coefficients[1], R(-1.5)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Comparison against reference
#ifdef HAS_BOOST_MATH
  std::println("Reference: Boost.Math (boost::math::legendre_p)");
#else
  std::println("Reference: Local fallback (recurrence relation)");
#endif

  double eps = std::numeric_limits<double>::epsilon();
  std::array<double, 7> test_points = {-1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0};

  auto run_check = [&](unsigned n, auto p) {
    std::println("--- P_{} ---", n);
    for (double xval : test_points)
    {
      double lam_val = p(xval);
#ifdef HAS_BOOST_MATH
      double ref_val = boost::math::legendre_p(n, xval);
#else
      double ref_val = reference_legendre(n, xval);
#endif
      double diff = lam_val - ref_val;
      double ulps = 0.0;
      if (std::abs(ref_val) > std::numeric_limits<double>::min())
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);

      std::println("x={:>5.2f} lam={:>10.5f} ref={:>10.5f} ulps={:>5.1f}", xval, lam_val, ref_val, ulps);
    }
  };

  run_check(3, P3);
  run_check(4, P4);

  // Benchmarking
  constexpr int iterations = 10'000'000;
  using Clock = std::chrono::steady_clock;
  volatile double sink = 0.0;
  constexpr auto P20 = lam::orthogonal::legendre_n<R, 20>();

  std::println("\n=== Performance Benchmark (P_20, 10M evals) ===");

  auto benchmark = [&](auto&& name, auto&& func) {
    auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
    {
      double xx = -1.0 + (i % 200) * 0.01;
      sink = func(xx);
    }
    auto end = Clock::now();
    std::chrono::duration<double> diff = end - start;
    std::println("{:<15}: {:.3f} s  ({:.1f} M/s)", name, diff.count(), iterations / diff.count() / 1e6);
  };

  benchmark("lam::legendre", [&](double val) { return P20(val); });
  benchmark("reference", [&](double val) {
#ifdef HAS_BOOST_MATH
    return boost::math::legendre_p(20, val);
#else
      return reference_legendre(20, val);
#endif
  });

  // Static Table Verification
  std::println("\n=== Static Table Verification (Exact Values) ===");
  struct Check
  {
    double x;
    double expected;
  };

  // P_2(x) = 1.5x^2 - 0.5
  constexpr auto P2_checks = std::to_array<Check>({
    {0.0, -0.5}, {1.0, 1.0}, {-1.0, 1.0}, {0.5, -0.125} // 1.5(0.25) - 0.5 = 0.375 - 0.5 = -0.125
  });

  // P_3(x) = 2.5x^3 - 1.5x
  constexpr auto P3_checks = std::to_array<Check>({
    {0.0, 0.0}, {1.0, 1.0}, {-1.0, -1.0}, {0.5, -0.4375} // 2.5(0.125) - 1.5(0.5) = 0.3125 - 0.75 = -0.4375
  });

  auto verify_table = [&](auto name, auto poly, auto checks) {
    bool all_passed = true;
    for (auto [x, expected] : checks)
    {
      double val = poly(x);
      double err = std::abs(val - expected);
      if (err > 1e-14)
      {
        std::println("FAIL: {}({}) = {:.16f}, expected {:.16f}", name, x, val, expected);
        all_passed = false;
      }
    }
    if (all_passed)
      std::println("{} passed {} table checks.", name, checks.size());
  };

  verify_table("P_3", P3, P3_checks);

  // Symmetry Verification: P_n(-x) = (-1)^n P_n(x)
  std::println("\n=== Symmetry Verification ===");
  auto verify_symmetry = [&](auto name, unsigned n, auto poly) {
    double x = 0.42; // arbitrary point
    double val_pos = poly(x);
    double val_neg = poly(-x);
    double expected_sign = (n % 2 == 0) ? 1.0 : -1.0;
    double diff = std::abs(val_neg - (expected_sign * val_pos));
    if (diff < 1e-14)
      std::println("{} symmetry passed (diff={:.1e})", name, diff);
    else
      std::println("{} symmetry FAILED (diff={:.1e})", name, diff);
  };
  verify_symmetry("P_2", 2, P2);
  verify_symmetry("P_3", 3, P3);

  // Differential Equation Verification (The "Best" Check)
  // Legendre ODE: (1-x^2)P_n'' - 2xP_n' + n(n+1)P_n = 0
  std::println("\n=== Differential Equation Verification (Analytic Proof) ===");

  auto verify_ode = [&](auto name, unsigned n, auto poly) {
    using T = R;
    constexpr auto x = lam::make_monomial<T, 1>();
    constexpr auto one = lam::polynomial_nttp<T, 0>{T(1)};

    // Calculate terms
    auto dP = lam::derivative(poly);
    auto d2P = lam::derivative(dP);

    // (1 - x^2) * P''
    auto term1 = (one - x * x) * d2P;

    // -2x * P'
    auto term2 = (T(-2) * x) * dP;

    // + n(n+1) * P
    T eig = T(n * (n + 1));
    auto term3 = eig * poly;

    auto result = term1 + term2 + term3;

    // Check if result is identically zero
    bool all_zero = true;
    for (auto c : result.coefficients)
    {
      if (std::abs(c) > 1e-13)
        all_zero = false;
    }

    if (all_zero)
      std::println("{} satisfies Legendre ODE analytically.", name);
    else
    {
      std::println("{} FAILED Legendre ODE check.", name);
      for (auto c : result.coefficients)
        std::print("{} ", c);
      std::println("");
    }
  };

  verify_ode("P_2 (n=2)", 2, P2);
  verify_ode("P_3 (n=3)", 3, P3);
  verify_ode("P_4 (n=4)", 4, P4);
  verify_ode("P_20 (n=20)", 20, P20); // Check high degree too!

  // Check for spherical harmonics (Associated Legendre) with even m
  // P_l^m(x) = (-1)^m (1-x^2)^(m/2) d^m/dx^m P_l(x)
  // For m positive even integer, (-1)^m = 1.
  // Example: m=2. P_l^2(x) = (1-x^2) * P_l''(x).
  
  std::println("\n=== Spherical Harmonics Demo (Associated Legendre, even m) ===");
  
  // Case 1: l=2, m=2
  // P_2(x) = 0.5(3x^2 - 1)
  // P_2'(x) = 3x
  // P_2''(x) = 3
  // P_2^2(x) = (1-x^2) * 3 = 3 - 3x^2
  // Reuse existing P2 from above
  constexpr auto P2_prime = lam::derivative(P2);
  constexpr auto P2_prime2 = lam::derivative(P2_prime);
  
  constexpr auto x_poly = lam::make_monomial<R, 1>();
  constexpr auto one = lam::polynomial_nttp<R, 0>{1.0};
  constexpr auto one_minus_x2 = one - x_poly * x_poly;
  
  constexpr auto P2_2 = one_minus_x2 * P2_prime2;
  
  std::println("P_2^2(x) = {}", P2_2);
  // Expected: 3 - 3x^2
  static_assert(lam::is_approx_equal(P2_2.coefficients[0], 3.0));
  static_assert(lam::is_approx_equal(P2_2.coefficients[2], -3.0));
  std::println("P_2^2 verified (3 - 3x^2)");

  // Case 2: l=3, m=2
  // P_3(x) = 0.5(5x^3 - 3x)
  // P_3'(x) = 0.5(15x^2 - 3)
  // P_3''(x) = 0.5(30x) = 15x
  // P_3^2(x) = (1-x^2) * 15x = 15x - 15x^3
  // Reuse existing P3 from above
  constexpr auto P3_prime = lam::derivative(P3);
  constexpr auto P3_prime2 = lam::derivative(P3_prime);
  
  constexpr auto P3_2 = one_minus_x2 * P3_prime2;
  
  std::println("P_3^2(x) = {}", P3_2);
  // Expected: 15x - 15x^3
  static_assert(lam::is_approx_equal(P3_2.coefficients[1], 15.0));
  static_assert(lam::is_approx_equal(P3_2.coefficients[3], -15.0));
  std::println("P_3^2 verified (15x - 15x^3)");

  return 0;
}
