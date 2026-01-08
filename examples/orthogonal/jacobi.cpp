/*
 *  jacobi.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Jacobi polynomial generation example
 *    Demonstrates generating Jacobi polynomials P_n^(alpha, beta)(x) at compile time
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
#include <boost/math/special_functions/jacobi.hpp>
#endif

namespace lam::orthogonal
{

// Generate the nth Jacobi polynomial P_n^(alpha, beta)(x)
// Recurrence relation for n >= 2:
//
// 2n(n + a + b)(2n + a + b - 2) P_n
//    = (2n + a + b - 1){ (2n + a + b)(2n + a + b - 2)x + (a^2 - b^2) } P_{n-1}
//    - 2(n + a - 1)(n + b - 1)(2n + a + b) P_{n-2}
//
// Base cases:
// P_0 = 1
// P_1 = (a - b)/2 + (1 + (a + b)/2 ) x  = 0.5 * ( (a - b) + (a + b + 2) * x )

template<typename R, std::size_t N, double Alpha, double Beta>
struct jacobi_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0)
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1)
    {
      // P_1 = 0.5 * (a - b) + 0.5 * (a + b + 2) * x
      constexpr R a = R(Alpha);
      constexpr R b = R(Beta);
      constexpr R c0 = R(0.5) * (a - b);
      constexpr R c1 = R(0.5) * (a + b + R(2));
      return lam::polynomial_nttp<R, 1>{{c0, c1}};
    }
    else
    {
      constexpr auto P_nm1 = jacobi_memo<R, N - 1, Alpha, Beta>::value;
      constexpr auto P_nm2 = jacobi_memo<R, N - 2, Alpha, Beta>::value;

      constexpr auto x = lam::make_monomial<R, 1>();

      constexpr R n = R(N);
      constexpr R a = R(Alpha);
      constexpr R b = R(Beta);

      // Terms for recurrence
      constexpr R apb = a + b;
      constexpr R two_n = R(2) * n;
      constexpr R denom = two_n * (n + apb) * (two_n + apb - R(2)); // LHS coefficient

      constexpr R factor1 = two_n + apb - R(1);
      constexpr R term_x_coeff = factor1 * (two_n + apb) * (two_n + apb - R(2));
      constexpr R term_const_coeff = factor1 * (a * a - b * b);

      constexpr R term_sub = R(2) * (n + a - R(1)) * (n + b - R(1)) * (two_n + apb);

      constexpr R inv_denom = R(1) / denom;

      // P_n = ( (term_const + term_x * x) * P_{n-1} - term_sub * P_{n-2} ) / denom
      return ((term_const_coeff + term_x_coeff * x) * P_nm1 - term_sub * P_nm2) * inv_denom;
    }
  }();
};

template<typename R, std::size_t N, double Alpha, double Beta>
constexpr auto jacobi_n()
{
  return jacobi_memo<R, N, Alpha, Beta>::value;
}

} // namespace lam::orthogonal

// Runtime Reference
double reference_jacobi(unsigned n, double alpha, double beta, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);

  double P_prev2 = 1.0;
  double P_prev1 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);

  for (unsigned k = 2; k <= n; ++k)
  {
    double k_d = static_cast<double>(k);

    // Using the same recurrence structure but solving for P_k
    // 2k(k + a + b)(2k + a + b - 2) P_k ...

    double apb = alpha + beta;
    double two_k = 2.0 * k_d;
    double denom = two_k * (k_d + apb) * (two_k + apb - 2.0);

    double factor1 = two_k + apb - 1.0;
    double term_x_part = (two_k + apb) * (two_k + apb - 2.0) * x;
    double term_const_part = (alpha * alpha - beta * beta);

    double term_sub = 2.0 * (k_d + alpha - 1.0) * (k_d + beta - 1.0) * (two_k + apb);

    double P_curr = (factor1 * (term_const_part + term_x_part) * P_prev1 - term_sub * P_prev2) / denom;

    P_prev2 = P_prev1;
    P_prev1 = P_curr;
  }
  return P_prev1;
}

int main()
{
  using R = double;

  // Test case 1: Alpha=0, Beta=0 (Should match Legendre P_n)

  constexpr double a1 = 0.0;
  constexpr double b1 = 0.0;

  // Test case 2: Alpha=0.5, Beta=0.5
  constexpr double a2 = 0.5;
  constexpr double b2 = 0.5;

  std::println("=== Jacobi Polynomials P_n^(a,b)(x) ===\n");

  constexpr auto P3_leg = lam::orthogonal::jacobi_n<R, 3, a1, b1>();
  constexpr auto P3_jac = lam::orthogonal::jacobi_n<R, 3, a2, b2>();

  std::println("Alpha=0, Beta=0 (Legendre):");
  std::println("P_3^(0,0) = {}", P3_leg); // Should be 2.5x^3 - 1.5x

  std::println("\nAlpha=0.5, Beta=0.5:");
  std::println("P_3^(0.5,0.5) = {}", P3_jac);

  // Compile-time verification for Legendre case
  // Compile-time verification for Legendre case
  static_assert(lam::is_approx_equal(P3_leg.coefficients[3], R(2.5)));
  static_assert(lam::is_approx_equal(P3_leg.coefficients[1], R(-1.5)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Comparison against reference
#ifdef HAS_BOOST_MATH
  std::println("Reference: Boost.Math (boost::math::jacobi_p)");
#else
  std::println("Reference: Local fallback (recurrence relation)");
#endif

  double eps = std::numeric_limits<double>::epsilon();
  std::array<double, 7> test_points = {-1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0};

  auto run_check = [&](unsigned n, double alpha, double beta, auto p) {
    std::println("--- P_{}^({},{}) ---", n, alpha, beta);
    for (double xval : test_points)
    {
      double lam_val = p(xval);
#ifdef HAS_BOOST_MATH
      double ref_val = boost::math::jacobi(n, alpha, beta, xval);
#else
      double ref_val = reference_jacobi(n, alpha, beta, xval);
#endif
      double diff = lam_val - ref_val;
      double ulps = 0.0;
      if (std::abs(ref_val) > std::numeric_limits<double>::min())
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);

      std::println("x={:>5.2f} lam={:>10.5f} ref={:>10.5f} ulps={:>5.1f}", xval, lam_val, ref_val, ulps);
    }
  };

  run_check(3, a1, b1, P3_leg);
  run_check(3, a2, b2, P3_jac);

  // Symmetry Verification: P_n^(a,a)(-x) = (-1)^n P_n^(a,a)(x)
  // Only valid here because we used alpha=beta for both test cases!
  std::println("\n=== Symmetry Verification (alpha=beta cases) ===");
  auto verify_symmetry = [&](auto name, unsigned n, auto poly) {
    double x = 0.8;
    double val_pos = poly(x);
    double val_neg = poly(-x);
    double expected_sign = (n % 2 == 0) ? 1.0 : -1.0;
    double diff = std::abs(val_neg - (expected_sign * val_pos));
    if (diff < 1e-13)
      std::println("{} symmetry passed (diff={:.1e})", name, diff);
    else
      std::println("{} symmetry FAILED (diff={:.1e})", name, diff);
  };
  verify_symmetry("P_3^(0,0)", 3, P3_leg);
  verify_symmetry("P_3^(0.5,0.5)", 3, P3_jac);

  // Differential Equation Verification (Analytic Proof)
  // Jacobi ODE: (1-x^2)y'' + [b - a - (a + b + 2)x]y' + n(n + a + b + 1)y = 0
  std::println("\n=== Differential Equation Verification (Analytic Proof) ===");

  auto verify_ode = [&](auto name, unsigned n, double alpha, double beta, auto poly) {
    using T = R;
    constexpr auto x = lam::make_monomial<T, 1>();
    constexpr auto one = lam::polynomial_nttp<T, 0>{T(1)};

    auto dy = lam::derivative(poly);
    auto d2y = lam::derivative(dy);

    // (1 - x^2) y''
    auto term1 = (one - x * x) * d2y;

    // [b - a - (a + b + 2)x] y'
    T b_minus_a = T(beta - alpha);
    T coeff_x = T(alpha + beta + 2.0);
    // poly term for y' coeff: (b-a) - (a+b+2)x
    // Constructing this poly explicitly:
    // We can use operators. Constant (b-a) needs to be a poly?
    // Our operator+(T, poly) works.
    // (b-a) - coeff_x * x
    auto y_prime_coeff_poly = b_minus_a - coeff_x * x;

    auto term2 = y_prime_coeff_poly * dy;

    // + n(n + a + b + 1) y
    T coeff_y = T(n * (n + alpha + beta + 1.0));
    auto term3 = coeff_y * poly;

    auto result = term1 + term2 + term3;

    bool all_zero = true;
    for (auto c : result.coefficients)
    {
      if (std::abs(c) > 1e-6)
        all_zero = false;
    }

    if (all_zero)
      std::println("{} satisfies Jacobi ODE analytically.", name);
    else
    {
      std::println("{} FAILED Jacobi ODE check.", name);
      for (auto c : result.coefficients)
        std::print("{} ", c);
      std::println("");
    }
  };

  verify_ode("P_3^(0,0)", 3, a1, b1, P3_leg);
  verify_ode("P_3^(0.5,0.5)", 3, a2, b2, P3_jac);
  constexpr auto P20_check = lam::orthogonal::jacobi_n<R, 20, a2, b2>();
  verify_ode("P_20^(0.5,0.5)", 20, a2, b2, P20_check);

  // Benchmarking
  constexpr int iterations = 10'000'000;
  using Clock = std::chrono::steady_clock;
  constexpr auto P20 = lam::orthogonal::jacobi_n<R, 20, a2, b2>();
  volatile double sink = 0.0;

  std::println("\n=== Performance Benchmark (P_20^(0.5,0.5), 10M evals) ===");

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

  benchmark("lam::jacobi", [&](double val) { return P20(val); });
  benchmark("reference", [&](double val) {
#ifdef HAS_BOOST_MATH
    return boost::math::jacobi(20, a2, b2, val);
#else
      return reference_jacobi(20, a2, b2, val);
#endif
  });

  return 0;
}
