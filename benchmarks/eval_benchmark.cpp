/*
 *  eval_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Benchmarks single-point evaluation performance across different compiler
 *    flags and backends. Compares:
 *      - Generic Horner's Method (Scalar)
 *      - SIMD/BLAS Accelerated Evaluation (if enabled)
 *      - Boost.Math (if enabled)
 */

#include <chrono>
#include <cmath>
#include <print>
#include <vector>

import lam.polynomial_nttp;

#ifdef HAS_BOOST_MATH
#include <boost/math/tools/polynomial.hpp>
#endif

template<typename R, std::size_t N>
constexpr auto make_test_polynomial()
{
  lam::polynomial_nttp<R, N> p{};
  R factorial = R(1);
  for (std::size_t i = 0; i <= N; ++i)
  {
    if (i > 0)
      factorial *= R(i);
    p.coefficients[i] = R(1) / factorial;
  }
  return p;
}

template<std::size_t N>
void benchmark_degree(int iterations, int num_values)
{
  constexpr auto poly = make_test_polynomial<double, N>();

  std::vector<double> values(num_values);
  for (int i = 0; i < num_values; ++i)
    values[i] = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(num_values - 1);

  volatile double sink = 0.0;

  auto start = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& x : values)
      sink = poly(x);
  auto end = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  double evals_per_sec = static_cast<double>(iterations * num_values) / (static_cast<double>(duration.count()) / 1e6);

  std::println("  Degree {:3}: {:7} µs ({:.2f}M evals/sec)", N, duration.count(), evals_per_sec / 1e6);

#ifdef HAS_BOOST_MATH
  std::vector<double> coeffs_vec(poly.coefficients.begin(), poly.coefficients.end());
  boost::math::tools::polynomial<double> boost_poly(coeffs_vec.data(), coeffs_vec.size() - 1);

  volatile double sink_boost = 0.0;
  double max_diff = 0.0;

  auto start_boost = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& x : values)
    {
      double b_val = boost_poly.evaluate(x);
      sink_boost = b_val;

      if (iter == 0)
      {
        double my_val = poly(x);
        double diff = std::abs(my_val - b_val);
        if (diff > max_diff)
          max_diff = diff;
      }
    }
  auto end_boost = std::chrono::steady_clock::now();
  auto boost_us = std::chrono::duration_cast<std::chrono::microseconds>(end_boost - start_boost).count();
  double boost_mps = static_cast<double>(iterations * num_values) / (static_cast<double>(boost_us) / 1e6);

  std::println("        Boost: {:7} µs ({:.2f}M evals/sec) | Speedup: {:.2f}x | MaxDiff: {:e}", boost_us,
               boost_mps / 1e6, static_cast<double>(boost_us) / static_cast<double>(duration.count()), max_diff);
#endif
}

int main()
{
  constexpr int iterations = 100'000;
  constexpr int num_values = 100;

  std::println("=== Polynomial Evaluation Benchmark (Horner's Method) ===");
  std::println("    {} iterations × {} values\n", iterations, num_values);

#ifdef HAS_BOOST_MATH
  std::println("    [Boost Enabled]");
#else
  std::println("    [Boost DISABLED]");
#endif

  std::println("--- Correctness Check ---");
  constexpr auto poly = make_test_polynomial<double, 20>();
  constexpr double at_one = poly(1.0);
  std::println("  poly exp(1) = {}", at_one);
  std::println("  std::exp(1) = {}", std::exp(1.0));
  std::println("  diff        = {}\n", std::abs(at_one - std::exp(1.0)));

  static_assert(poly(0.0) == 1.0);
  std::println("  ✓ Compile-time evaluation verified");

  std::println("\n--- Backend Consistency Check ---");
  double val_parallel = poly.evaluate_parallel(1.0);
  double val_backend = poly(1.0);
  std::println("  Fallback : {}", val_parallel);
  std::println("  Backend  : {}", val_backend);
  double consistency_diff = std::abs(val_parallel - val_backend);
  if (consistency_diff < 1e-14)
  {
    std::println("  ✓ Backend matches Fallback (diff = {})", consistency_diff);
  }
  else
  {
    std::println("  ✗ Backend MISMATCH (diff = {})", consistency_diff);
    std::exit(1);
  }
  std::println("");

  std::println("--- Benchmark by Degree ---");
  benchmark_degree<5>(iterations, num_values);
  benchmark_degree<10>(iterations, num_values);
  benchmark_degree<15>(iterations, num_values);
  benchmark_degree<20>(iterations, num_values);
  benchmark_degree<25>(iterations, num_values);
  benchmark_degree<30>(iterations, num_values);
  benchmark_degree<50>(iterations, num_values);
  benchmark_degree<100>(iterations, num_values);

  std::println("\n--- vs std::exp ---");
  std::vector<double> values(num_values);
  for (int i = 0; i < num_values; ++i)
    values[i] = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(num_values - 1);

  volatile double sink = 0.0;

  auto start_std = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& x : values)
      sink = std::exp(x);
  auto end_std = std::chrono::steady_clock::now();
  auto std_us = std::chrono::duration_cast<std::chrono::microseconds>(end_std - start_std).count();

  std::println("  std::exp: {} µs", std_us);

  return 0;
}
