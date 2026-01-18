/*
 *  laguerre.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Laguerre polynomial generation example
 *    Demonstrates generating the first N Laguerre polynomials at compile time
 *    using recursive constexpr functions and std::index_sequence
 *
 *    Laguerre recurrence: L_{n+1}(x) = ((2n+1-x)L_n(x) - n*L_{n-1}(x)) / (n+1)
 */

#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <print>
#include <vector>

import lam.polynomial_nttp;

#ifdef HAS_BOOST_MATH
#include <boost/math/special_functions/laguerre.hpp>
#endif

// Error-free transformation: TwoProductFMA
// Computes x*y = p + e exactly
template<typename T>
constexpr std::pair<T, T> two_prod(T x, T y)
{
  T p = x * y;
  T e = std::fma(x, y, -p);
  return {p, e};
}

// Error-free transformation: TwoSum
// Computes a + b = s + e exactly
template<typename T>
constexpr std::pair<T, T> two_sum(T a, T b)
{
  T s = a + b;
  T z = s - a;
  T e = (a - (s - z)) + (b - z);
  return {s, e};
}

// Compensated Horner's Method
// Returns close to the precision of the type despite cancellation
template<typename R, std::size_t N>
constexpr R compensated_horner(const lam::polynomial_nttp<R, N>& poly, R x)
{
  R s = poly[N];
  R r = 0;

  for (int i = static_cast<int>(N) - 1; i >= 0; --i)
  {
    auto [p, pi] = two_prod(s, x);
    auto [sum, sigma] = two_sum(p, poly[i]);
    s = sum;
    r = r * x + (pi + sigma);
  }
  return s + r;
}

namespace lam::orthogonal
{

// Generate the nth Laguerre polynomial L_n(x)
// Uses the recurrence relation:
//   L_0(x) = 1
//   L_1(x) = 1 - x
//   L_{n+1}(x) = ((2n+1 - x) * L_n(x) - n * L_{n-1}(x)) / (n+1)
template<typename R, std::size_t N>
struct laguerre_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0) // L_0(x) = 1
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1) // L_1(x) = 1 - x
      return lam::polynomial_nttp<R, 1>{{R(1), R(-1)}};
    else
    {
      // Recursive case: need L_{n-1} and L_{n-2} (memoized)
      constexpr auto L_nm1 = laguerre_memo<R, N - 1>::value;
      constexpr auto L_nm2 = laguerre_memo<R, N - 2>::value;

      // x as a polynomial
      constexpr auto x = lam::make_monomial<R, 1>();

      // L_n = ((2n-1 - x) * L_{n-1} - (n-1) * L_{n-2}) / n
      // Note: when computing L_N, we use n = N-1 in the shifted recurrence
      constexpr R two_n_minus_1 = R(2 * (N - 1) + 1);
      constexpr R n_minus_1 = R(N - 1);
      constexpr R inv_n = R(1) / R(N);

      // (2n-1) * L_{n-1}
      constexpr auto term1 = two_n_minus_1 * L_nm1;
      // x * L_{n-1}
      constexpr auto term2 = x * L_nm1;
      // (n-1) * L_{n-2}
      constexpr auto term3 = n_minus_1 * L_nm2;

      // ((2n-1)*L_{n-1} - x*L_{n-1} - (n-1)*L_{n-2}) / n
      return (term1 - term2 - term3) * inv_n;
    }
  }();
};

template<typename R, std::size_t N>
constexpr auto laguerre_n()
{
  return laguerre_memo<R, N>::value;
}

// Helper to generate a tuple of Laguerre polynomials L_0 through L_{N-1}
template<typename R, std::size_t... Is>
constexpr auto make_laguerre_tuple_impl(std::index_sequence<Is...>)
{
  return std::make_tuple(laguerre_n<R, Is>()...);
}

// Generate tuple containing the first N Laguerre polynomials
template<typename R, std::size_t N>
constexpr auto first_n_laguerre()
{
  return make_laguerre_tuple_impl<R>(std::make_index_sequence<N>{});
}

// Helper to print tuple elements with their index using std::formatter
template<typename Tuple, std::size_t... Is>
void print_laguerre_tuple_impl(const Tuple& t, std::index_sequence<Is...>)
{
  (std::println("L_{} = {}", Is, std::get<Is>(t)), ...);
}

template<std::size_t N>
void print_laguerre_tuple(const auto& t)
{
  print_laguerre_tuple_impl(t, std::make_index_sequence<N>{});
}

} // namespace lam::orthogonal

// Runtime Laguerre implementation for comparison
// (Used as std::laguerre is unavailable in this environment)
double reference_laguerre(unsigned n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return 1.0 - x;

  double L_prev2 = 1.0;     // L_0
  double L_prev1 = 1.0 - x; // L_1

  for (unsigned k = 1; k < n; ++k)
  {
    // L_{k+1} = ((2k+1 - x)L_k - k*L_{k-1}) / (k+1)
    double L_curr = ((2 * k + 1 - x) * L_prev1 - k * L_prev2) / (k + 1);
    L_prev2 = L_prev1;
    L_prev1 = L_curr;
  }
  return L_prev1;
}

int main()
{
  using R = double;

  std::println("=== Laguerre Polynomials (Compile-Time Generation) ===\n");

  // Generate individual polynomials
  constexpr auto L0 = lam::orthogonal::laguerre_n<R, 0>();
  constexpr auto L1 = lam::orthogonal::laguerre_n<R, 1>();
  constexpr auto L2 = lam::orthogonal::laguerre_n<R, 2>();
  constexpr auto L3 = lam::orthogonal::laguerre_n<R, 3>();
  constexpr auto L4 = lam::orthogonal::laguerre_n<R, 4>();
  constexpr auto L5 = lam::orthogonal::laguerre_n<R, 5>();
  constexpr auto L10 = lam::orthogonal::laguerre_n<R, 10>();
  constexpr auto L20 = lam::orthogonal::laguerre_n<R, 20>();

  std::println("Individual Laguerre polynomials:");
  std::println("L_0 = {}", L0);
  std::println("L_1 = {}", L1);
  std::println("L_2 = {}", L2);
  std::println("L_3 = {}", L3);
  std::println("L_4 = {}", L4);
  std::println("L_5 = {}", L5);

  // Verify at compile time using static_assert
  // L_2(x) = 1 - 2x + (1/2)x^2
  // L_2(x) = 1 - 2x + (1/2)x^2
  static_assert(lam::is_approx_equal(L2.coefficients[0], R(1)));
  static_assert(lam::is_approx_equal(L2.coefficients[1], R(-2)));
  static_assert(lam::is_approx_equal(L2.coefficients[2], R(0.5)));

  // L_3(x) = 1 - 3x + (3/2)x^2 - (1/6)x^3
  static_assert(lam::is_approx_equal(L3.coefficients[0], R(1)));
  static_assert(lam::is_approx_equal(L3.coefficients[1], R(-3)));
  static_assert(lam::is_approx_equal(L3.coefficients[2], R(1.5)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Generate the first N polynomials as a tuple
  std::println("Generating first 8 Laguerre polynomials as a tuple:");
  constexpr auto laguerre_tuple = lam::orthogonal::first_n_laguerre<R, 8>();
  lam::orthogonal::print_laguerre_tuple<8>(laguerre_tuple);

  // Demonstrate formatting options
  std::println("\n--- Formatting examples ---");
  std::println("Default:      {}", L3);
  std::println("Variable 't': {:t}", L3);
  std::println("Precision 2:  {:.2}", L3);

  // Evaluate polynomials at a point
  std::println("\n--- Evaluation at x = 1.0 ---");
  constexpr R x = 1.0;
  std::println("L_0({}) = {}", x, L0(x));
  std::println("L_1({}) = {}", x, L1(x));
  std::println("L_2({}) = {}", x, L2(x));
  std::println("L_3({}) = {}", x, L3(x));
  std::println("L_4({}) = {}", x, L4(x));
  std::println("L_5({}) = {}", x, L5(x));

  // Note: Laguerre polynomials satisfy L_n(0) = 1 for all n
  std::println("\n--- Verifying L_n(0) = 1 ---");
  static_assert(lam::is_approx_equal(L0(0.0), 1.0));
  static_assert(lam::is_approx_equal(L1(0.0), 1.0));
  static_assert(lam::is_approx_equal(L2(0.0), 1.0));
  static_assert(lam::is_approx_equal(L3(0.0), 1.0));
  static_assert(lam::is_approx_equal(L4(0.0), 1.0));
  static_assert(lam::is_approx_equal(L5(0.0), 1.0));
  std::println("All L_n(0) = 1 verified at compile time!");

  // Compare against reference implementation
  std::println("\n=== Comparison with reference_laguerre (runtime) ===\n");

  double eps = std::numeric_limits<double>::epsilon();
  std::println("Machine epsilon: {:.4e}", eps);
  std::println("Note: 1 ULP is roughly epsilon * |value|\n");

  std::array<double, 7> test_points = {0.0, 0.5, 1.0, 2.0, 3.5, 5.0, 10.0};

  std::println("{:>4} {:>10} {:>20} {:>20} {:>15} {:>10}", "n", "x", "lam::", "ref", "diff", "ulps");
  std::println("{:-<85}", "");

  auto run_check = [&](unsigned n, auto p) {
    for (double xval : test_points)
    {
      double lam_val = p(xval);
#ifdef HAS_BOOST_MATH
      double ref_val = boost::math::laguerre(n, xval);
#else
      double ref_val = reference_laguerre(n, xval);
#endif
      double diff = lam_val - ref_val;

      double ulps = 0.0;
      if (ref_val != 0.0)
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);
      else if (diff != 0.0)
        ulps = std::abs(diff) / std::numeric_limits<double>::min();

      std::println("{:>4} {:>10.2f} {:>20.12f} {:>20.12f} {:>15.2e} {:>10.1f}", n, xval, lam_val, ref_val, diff, ulps);
    }
    std::println("");
  };

  run_check(0, L0);
  run_check(1, L1);
  run_check(2, L2);
  run_check(3, L3);
  run_check(4, L4);
  run_check(5, L5);
  constexpr auto L50 = lam::orthogonal::laguerre_n<R, 50>();
  run_check(50, L50);

  std::println("\n--- Summary ---");
  std::println("lam::laguerre: compile-time coefficients, Horner's evaluation");
#ifdef HAS_BOOST_MATH
  std::println("ref_laguerre:  boost::math::laguerre (verified external library)");
#else
  std::println("ref_laguerre:  runtime recurrence evaluation (local fallback)");
#endif

  std::println("\n=== Performance Benchmark (N=5, 10,000,000 evals) ===");
  using Clock = std::chrono::steady_clock;
  constexpr int iterations = 10'000'000;
  volatile double sink = 0.0; // prevent optimization

  auto benchmark = [&](auto&& name, auto&& func) {
    auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
    {
      double x = 1.0 + (i % 100) * 0.01;
      sink = func(x);
    }
    auto end = Clock::now();
    std::chrono::duration<double> diff = end - start;
    std::println("{:<15}: {:.3f} s  ({:.1f} M/s)", name, diff.count(), iterations / diff.count() / 1e6);
  };

  benchmark("lam::laguerre", [&](double x) { return L5(x); });
  benchmark("reference", [&](double x) {
#ifdef HAS_BOOST_MATH
    return boost::math::laguerre(5, x);
#else
      return reference_laguerre(5, x);
#endif
  });

  std::println("\n=== Performance Benchmark (N=10, 10,000,000 evals) ===");
  benchmark("lam::laguerre", [&](double x) { return L10(x); });
  benchmark("reference", [&](double x) {
#ifdef HAS_BOOST_MATH
    return boost::math::laguerre(10, x);
#else
      return reference_laguerre(10, x);
#endif
  });

  std::println("\n=== Performance Benchmark (N=20, 10,000,000 evals) ===");
  benchmark("lam::laguerre", [&](double x) { return L20(x); });
  benchmark("reference", [&](double x) {
#ifdef HAS_BOOST_MATH
    return boost::math::laguerre(20, x);
#else
      return reference_laguerre(20, x);
#endif
  });

  benchmark("compensated", [&](double x) { return compensated_horner(L20, x); });

  // Stability check for N=20
  std::println("\n--- Stability Check (N=20, x=10.0) ---");
  double val_20 = L20(10.0);
  double comp_20 = compensated_horner(L20, 10.0);
#ifdef HAS_BOOST_MATH
  double ref_20 = boost::math::laguerre(20, 10.0);
#else
  double ref_20 = reference_laguerre(20, 10.0);
#endif
  std::println("lam:  {:.16e}", val_20);
  std::println("comp: {:.16e}", comp_20);
  std::println("ref:  {:.16e}", ref_20);
  std::println("err (lam):  {:.2e}", std::abs(val_20 - ref_20));
  std::println("err (comp): {:.2e}", std::abs(comp_20 - ref_20));

  return 0;
}
