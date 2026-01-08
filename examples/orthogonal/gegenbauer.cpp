/*
 *  gegenbauer.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Gegenbauer (Ultraspherical) polynomial generation example
 *    Demonstrates generating Gegenbauer polynomials C_n^(alpha)(x) at compile time
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
#include <boost/math/special_functions/gegenbauer.hpp>
#endif

namespace lam::orthogonal
{

// Generate the nth Gegenbauer polynomial C_n^(alpha)(x)
// (Also known as Ultraspherical polynomials)
//
// Recurrence:
//   C_0^(alpha)(x) = 1
//   C_1^(alpha)(x) = 2 * alpha * x
//   n * C_n^(alpha)(x) = 2(n + alpha - 1)x C_{n-1}^(alpha)(x) - (n + 2*alpha - 2)C_{n-2}^(alpha)(x)
//
//   => C_n = (2(n+a-1)/n) * x * C_{n-1} - ((n+2a-2)/n) * C_{n-2}
//
// Special case: if alpha = 0, C_n^0(x) are identically 0 for n >= 1 (or defined via limits related to Chebyshev T_n).
// Boost handles alpha=0 by returning 0 for n>=1. We will assume alpha != 0 or handle it similarly.

// Note: Alpha must be a compile-time constant for our template approach,
// OR we can fix alpha for the type.
// Here we will use a template parameter for Alpha.
// Since Alpha is double, we can't pass it easily as a non-type template param in C++20/23 standardly
// unless it's a structural type or we wrappers.
// For simplicity in this example, we will fix Alpha or generate a specific Alpha.
// Let's implement a generator where Alpha is passed as a template argument (using a ratio or similar if needed),
// BUT typically these libraries want `double alpha`.
//
// HOWEVER, `lam::polynomial_nttp` is a compile-time structure.
// If we want `C_n^{alpha}` where alpha is known at compile time, we can bake it into coefficients.
// If alpha is runtime, we can't fully bake it into `polynomial_nttp` (which stores fixed coeffs).
//
// The prompt implies we want compile-time generation. Let's fix Alpha = 1 (Chebyshev U) or Alpha = 0.5 (Legendre)
// or just have a template `double Alpha` if supported (C++20 allows double NTTP).

template<typename R, std::size_t N, double Alpha>
struct gegenbauer_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0)
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1)
    {
      // C_1^a(x) = 2 * alpha * x
      constexpr R two_alpha = R(2.0 * Alpha);
      return lam::polynomial_nttp<R, 1>{{R(0), two_alpha}};
    }
    else
    {
      constexpr auto C_nm1 = gegenbauer_memo<R, N - 1, Alpha>::value;
      constexpr auto C_nm2 = gegenbauer_memo<R, N - 2, Alpha>::value;

      constexpr auto x = lam::make_monomial<R, 1>();

      constexpr R n = R(N);
      // n * C_n = 2(n + a - 1)x C_{n-1} - (n + 2a - 2)C_{n-2}

      constexpr R term1_coeff = R(2.0 * (n + Alpha - 1.0));
      constexpr R term2_coeff = R(n + 2.0 * Alpha - 2.0);

      constexpr R inv_n = R(1.0) / n;

      return (term1_coeff * x * C_nm1 - term2_coeff * C_nm2) * inv_n;
    }
  }();
};

template<typename R, std::size_t N, double Alpha>
constexpr auto gegenbauer_n()
{
  return gegenbauer_memo<R, N, Alpha>::value;
}

} // namespace lam::orthogonal

// Runtime Reference
double reference_gegenbauer(unsigned n, double alpha, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return 2.0 * alpha * x;

  double C_prev2 = 1.0;
  double C_prev1 = 2.0 * alpha * x;

  for (unsigned k = 1; k < n; ++k)
  {
    // C_{k+1} = [2(k + a)x C_k - (k + 2a - 1)C_{k-1}] / (k+1)
    // Note: The recurrence index mapping can be tricky.
    // Standard: (n+1) C_{n+1} = 2(n+a) x C_n - (n+2a-1) C_{n-1}
    // Let's match the one used above:
    // n C_n = 2(n+a-1)x C_{n-1} - (n+2a-2)C_{n-2}
    // Mapping k -> n-1 implies (k+1) C_{k+1} = 2(k+1+a-1) x C_k ... = 2(k+a) x C_k

    double term1 = 2.0 * (k + alpha) * x * C_prev1;
    double term2 = (k + 2.0 * alpha - 1.0) * C_prev2;
    double C_curr = (term1 - term2) / (k + 1.0);

    C_prev2 = C_prev1;
    C_prev1 = C_curr;
  }
  return C_prev1;
}

int main()
{
  using R = double;

  // Test with Alpha = 1.0 (Should match Chebyshev U_n, but scale might differ? No, U_n is C_n^1)
  // Test with Alpha = 0.5 (Should match Legendre P_n)

  constexpr double alpha1 = 1.0;
  constexpr double alpha2 = 2.5; // Arbitrary parameter

  std::println("=== Gegenbauer Polynomials C_n^a(x) ===\n");

  constexpr auto C3_a1 = lam::orthogonal::gegenbauer_n<R, 3, alpha1>(); // Expect U_3 = 8x^3 - 4x
  constexpr auto C3_a2 = lam::orthogonal::gegenbauer_n<R, 3, alpha2>();

  std::println("Alpha = 1.0 (Chebyshev U_n):");
  std::println("C_3^1 = {}", C3_a1);

  std::println("\nAlpha = 2.5:");
  std::println("C_3^2.5 = {}", C3_a2);

  // Compile-time verification for A=1 (U_3)
  // Compile-time verification for A=1 (U_3)
  static_assert(lam::is_approx_equal(C3_a1.coefficients[3], R(8)));
  static_assert(lam::is_approx_equal(C3_a1.coefficients[1], R(-4)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Comparison against reference
#ifdef HAS_BOOST_MATH
  std::println("Reference: Boost.Math (boost::math::gegenbauer)");
#else
  std::println("Reference: Local fallback (recurrence relation)");
#endif

  double eps = std::numeric_limits<double>::epsilon();
  std::array<double, 7> test_points = {-1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0};

  auto run_check = [&](unsigned n, double alpha, auto p) {
    std::println("--- C_{}^{} ---", n, alpha);
    for (double xval : test_points)
    {
      double lam_val = p(xval);
#ifdef HAS_BOOST_MATH
      double ref_val = boost::math::gegenbauer(n, alpha, xval);
#else
      double ref_val = reference_gegenbauer(n, alpha, xval);
#endif
      double diff = lam_val - ref_val;
      double ulps = 0.0;
      if (std::abs(ref_val) > std::numeric_limits<double>::min())
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);

      std::println("x={:>5.2f} lam={:>10.5f} ref={:>10.5f} ulps={:>5.1f}", xval, lam_val, ref_val, ulps);
    }
  };

  run_check(3, alpha1, C3_a1);
  run_check(3, alpha2, C3_a2);

  // Symmetry Verification: C_n^a(-x) = (-1)^n C_n^a(x)
  std::println("\n=== Symmetry Verification ===");
  auto verify_symmetry = [&](auto name, unsigned n, auto poly) {
    double x = 0.5;
    double val_pos = poly(x);
    double val_neg = poly(-x);
    double expected_sign = (n % 2 == 0) ? 1.0 : -1.0;
    double diff = std::abs(val_neg - (expected_sign * val_pos));
    if (diff < 1e-13)
      std::println("{} symmetry passed (diff={:.1e})", name, diff);
    else
      std::println("{} symmetry FAILED (diff={:.1e})", name, diff);
  };
  verify_symmetry("C_3^1", 3, C3_a1);
  verify_symmetry("C_3^2.5", 3, C3_a2);

  // Differential Equation Verification (Analytic Proof)
  // Gegenbauer ODE: (1-x^2)y'' - (2a + 1)xy' + n(n + 2a)y = 0
  std::println("\n=== Differential Equation Verification (Analytic Proof) ===");

  auto verify_ode = [&](auto name, unsigned n, double alpha, auto poly) {
    using T = R;
    constexpr auto x = lam::make_monomial<T, 1>();
    constexpr auto one = lam::polynomial_nttp<T, 0>{T(1)};

    auto dy = lam::derivative(poly);
    auto d2y = lam::derivative(dy);

    // (1 - x^2) y''
    auto term1 = (one - x * x) * d2y;

    // - (2a + 1) x y'
    T coeff_y_prime = T(2.0 * alpha + 1.0);
    auto term2 = (coeff_y_prime * x) * dy; // subtract this

    // + n(n + 2a) y
    T coeff_y = T(n * (n + 2.0 * alpha));
    auto term3 = coeff_y * poly;

    auto result = term1 - term2 + term3;

    bool all_zero = true;
    for (auto c : result.coefficients)
    {
      if (std::abs(c) > 1e-12)
        all_zero = false;
    }

    if (all_zero)
      std::println("{} satisfies Gegenbauer ODE analytically.", name);
    else
    {
      std::println("{} FAILED Gegenbauer ODE check.", name);
      for (auto c : result.coefficients)
        std::print("{} ", c);
      std::println("");
    }
  };

  verify_ode("C_3^1 (n=3, a=1)", 3, alpha1, C3_a1);
  verify_ode("C_3^2.5 (n=3, a=2.5)", 3, alpha2, C3_a2);
  // C20 from benchmark
  constexpr auto C20_check = lam::orthogonal::gegenbauer_n<R, 20, alpha2>();
  verify_ode("C_20^2.5 (n=20, a=2.5)", 20, alpha2, C20_check);

  // Benchmarking
  constexpr int iterations = 10'000'000;
  using Clock = std::chrono::steady_clock;
  constexpr auto C20 = lam::orthogonal::gegenbauer_n<R, 20, alpha2>();
  volatile double sink = 0.0;

  std::println("\n=== Performance Benchmark (C_20^2.5, 10M evals) ===");

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

  benchmark("lam::gegenb", [&](double val) { return C20(val); });
  benchmark("reference", [&](double val) {
#ifdef HAS_BOOST_MATH
    return boost::math::gegenbauer(20, alpha2, val);
#else
      return reference_gegenbauer(20, alpha2, val);
#endif
  });

  return 0;
}
