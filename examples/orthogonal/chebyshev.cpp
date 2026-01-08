/*
 *  chebyshev.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Chebyshev polynomial generation example
 *    Demonstrates generating Chebyshev polynomials of First and Second kind
 *    at compile time using recursive constexpr functions.
 */

#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <print>
#include <vector>

import lam.polynomial_nttp;

#ifdef HAS_BOOST_MATH
#include <boost/math/special_functions/chebyshev.hpp>
#endif

namespace lam::orthogonal
{

// Chebyshev polynomials of the First Kind: T_n(x)
// Recurrence:
//   T_0(x) = 1
//   T_1(x) = x
//   T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)
template<typename R, std::size_t N>
struct chebyshev_t_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0) // T_0(x) = 1
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1) // T_1(x) = x
      return lam::polynomial_nttp<R, 1>{{R(0), R(1)}};
    else
    {
      // Recursive case
      constexpr auto T_nm1 = chebyshev_t_memo<R, N - 1>::value;
      constexpr auto T_nm2 = chebyshev_t_memo<R, N - 2>::value;

      // x as a polynomial
      constexpr auto x = lam::make_monomial<R, 1>();

      // T_{n+1} = 2x T_n - T_{n-1}
      constexpr R two = R(2);
      return two * x * T_nm1 - T_nm2;
    }
  }();
};

template<typename R, std::size_t N>
constexpr auto chebyshev_t_n()
{
  return chebyshev_t_memo<R, N>::value;
}

// Chebyshev polynomials of the Second Kind: U_n(x)
// Recurrence:
//   U_0(x) = 1
//   U_1(x) = 2x
//   U_{n+1}(x) = 2x U_n(x) - U_{n-1}(x)
template<typename R, std::size_t N>
struct chebyshev_u_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0) // U_0(x) = 1
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1) // U_1(x) = 2x
      return lam::polynomial_nttp<R, 1>{{R(0), R(2)}};
    else
    {
      // Recursive case
      constexpr auto U_nm1 = chebyshev_u_memo<R, N - 1>::value;
      constexpr auto U_nm2 = chebyshev_u_memo<R, N - 2>::value;

      // x as a polynomial
      constexpr auto x = lam::make_monomial<R, 1>();

      // U_{n+1} = 2x U_n - U_{n-1}
      constexpr R two = R(2);
      return two * x * U_nm1 - U_nm2;
    }
  }();
};

template<typename R, std::size_t N>
constexpr auto chebyshev_u_n()
{
  return chebyshev_u_memo<R, N>::value;
}

} // namespace lam::orthogonal

// Runtime Reference Implementations
double reference_chebyshev_t(unsigned n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;

  double T_prev2 = 1.0;
  double T_prev1 = x;

  for (unsigned k = 1; k < n; ++k)
  {
    double T_curr = 2.0 * x * T_prev1 - T_prev2;
    T_prev2 = T_prev1;
    T_prev1 = T_curr;
  }
  return T_prev1;
}

double reference_chebyshev_u(unsigned n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return 2.0 * x;

  double U_prev2 = 1.0;
  double U_prev1 = 2.0 * x;

  for (unsigned k = 1; k < n; ++k)
  {
    double U_curr = 2.0 * x * U_prev1 - U_prev2;
    U_prev2 = U_prev1;
    U_prev1 = U_curr;
  }
  return U_prev1;
}

int main()
{
  using R = double;

  std::println("=== Chebyshev Polynomials ===\n");

  constexpr auto T3 = lam::orthogonal::chebyshev_t_n<R, 3>();
  constexpr auto T4 = lam::orthogonal::chebyshev_t_n<R, 4>();
  constexpr auto U3 = lam::orthogonal::chebyshev_u_n<R, 3>();
  constexpr auto U4 = lam::orthogonal::chebyshev_u_n<R, 4>();

  std::println("First Kind T_n(x):");
  std::println("T_3 = {}", T3); // 4x^3 - 3x
  std::println("T_4 = {}", T4); // 8x^4 - 8x^2 + 1

  std::println("\nSecond Kind U_n(x):");
  std::println("U_3 = {}", U3); // 8x^3 - 4x
  std::println("U_4 = {}", U4); // 16x^4 - 12x^2 + 1

  // Compile-time verification
  static_assert(lam::is_approx_equal(T3.coefficients[3], R(4)));
  static_assert(lam::is_approx_equal(T3.coefficients[1], R(-3)));
  static_assert(lam::is_approx_equal(T4.coefficients[0], R(1)));

  static_assert(lam::is_approx_equal(U3.coefficients[3], R(8)));
  static_assert(lam::is_approx_equal(U3.coefficients[1], R(-4)));
  static_assert(lam::is_approx_equal(U4.coefficients[0], R(1)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Comparison against reference
#ifdef HAS_BOOST_MATH
  std::println("Reference: Boost.Math (boost::math::chebyshev_t/u)");
#else
  std::println("Reference: Local fallback (recurrence relation)");
#endif

  double eps = std::numeric_limits<double>::epsilon();
  std::array<double, 7> test_points = {-1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0};

  auto run_check = [&](auto name, unsigned n, auto p, auto ref_func) {
    std::println("--- {} (n={}) ---", name, n);
    for (double xval : test_points)
    {
      double lam_val = p(xval);
      double ref_val = ref_func(n, xval);
      double diff = lam_val - ref_val;
      double ulps = 0.0;
      if (std::abs(ref_val) > std::numeric_limits<double>::min())
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);

      std::println("x={:>5.2f} lam={:>10.5f} ref={:>10.5f} ulps={:>5.1f}", xval, lam_val, ref_val, ulps);
    }
  };

  auto ref_t = [&](unsigned n, double x) {
#ifdef HAS_BOOST_MATH
    return boost::math::chebyshev_t(n, x);
#else
    return reference_chebyshev_t(n, x);
#endif
  };

  auto ref_u = [&](unsigned n, double x) {
#ifdef HAS_BOOST_MATH
    return boost::math::chebyshev_u(n, x);
#else
    return reference_chebyshev_u(n, x);
#endif
  };

  run_check("T", 3, T3, ref_t);
  run_check("U", 3, U3, ref_u);

  // Benchmarking
  constexpr int iterations = 10'000'000;
  using Clock = std::chrono::steady_clock;
  volatile double sink = 0.0;
  constexpr auto T20 = lam::orthogonal::chebyshev_t_n<R, 20>();

  std::println("\n=== Performance Benchmark (T_20, 10M evals) ===");

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

  benchmark("lam::chebyshev_t", [&](double val) { return T20(val); });
  benchmark("reference", [&](double val) { return ref_t(20, val); });

  // Static Table Verification
  std::println("\n=== Static Table Verification (Exact Values) ===");
  struct Check
  {
    double x;
    double expected;
  };

  // T_3(x) = 4x^3 - 3x
  constexpr auto T3_checks = std::to_array<Check>({
    {0.0, 0.0}, {1.0, 1.0}, {-1.0, -1.0}, {0.5, -1.0} // 4(0.125) - 3(0.5) = 0.5 - 1.5 = -1.0
  });

  // U_3(x) = 8x^3 - 4x
  constexpr auto U3_checks = std::to_array<Check>({
    {0.0, 0.0}, {1.0, 4.0}, {-1.0, -4.0}, {0.5, -1.0} // 8(0.125) - 4(0.5) = 1 - 2 = -1
  });

  auto verify_table = [&](auto name, auto poly, auto checks) {
    bool all_passed = true;
    for (auto [x, expected] : checks)
    {
      double val = poly(x);
      double err = std::abs(val - expected);
      constexpr auto tol = 1e-14;
      if (err > tol)
      {
        std::println("FAIL: {}({}) = {:.16f}, expected {:.16f}", name, x, val, expected);
        all_passed = false;
      }
    }
    if (all_passed)
      std::println("{} passed {} table checks.", name, checks.size());
  };

  verify_table("T_3", T3, T3_checks);
  verify_table("U_3", U3, U3_checks);

  // Differential Equation Verification (Analytic Proof)
  std::println("\n=== Differential Equation Verification (Analytic Proof) ===");

  using T = R;
  constexpr auto x = lam::make_monomial<T, 1>();
  constexpr auto one = lam::polynomial_nttp<T, 0>{T(1)};

  // Chebyshev T ODE: (1-x^2)T_n'' - xT_n' + n^2T_n = 0
  auto verify_ode_T = [&](auto name, unsigned n, auto poly) {
    auto dT = lam::derivative(poly);
    auto d2T = lam::derivative(dT);

    auto result = (one - x * x) * d2T - x * dT + T(n * n) * poly;

    bool all_zero = true;
    constexpr auto tol = 1e-12;
    for (auto c : result.coefficients)
      if (std::abs(c) > tol)
        all_zero = false;
    if (all_zero)
      std::println("{} satisfies Chebyshev T ODE analytically.", name);
    else
      std::println("{} FAILED Chebyshev T ODE check.", name);
  };

  // Chebyshev U ODE: (1-x^2)U_n'' - 3xU_n' + n(n+2)U_n = 0
  auto verify_ode_U = [&](auto name, unsigned n, auto poly) {
    auto dU = lam::derivative(poly);
    auto d2U = lam::derivative(dU);

    auto result = (one - x * x) * d2U - (T(3) * x) * dU + T(n * (n + 2)) * poly;

    bool all_zero = true;
    constexpr auto tol = 1e-12;
    for (auto c : result.coefficients)
      if (std::abs(c) > tol)
        all_zero = false;
    if (all_zero)
      std::println("{} satisfies Chebyshev U ODE analytically.", name);
    else
      std::println("{} FAILED Chebyshev U ODE check.", name);
  };

  verify_ode_T("T_3", 3, T3);
  verify_ode_U("U_3", 3, U3);
  verify_ode_T("T_20", 20, T20);
  // Re-generate U20 for check since it was only in benchmark
  constexpr auto U20 = lam::orthogonal::chebyshev_u_n<R, 20>();
  verify_ode_U("U_20", 20, U20);

  return 0;
}
