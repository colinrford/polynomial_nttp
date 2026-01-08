/*
 *  eval_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Benchmarks Horner's method for various polynomial degrees
 */

import std;
import lam.polynomial_nttp;

// Generate a polynomial with coefficients 1/i! for testing (exp series)
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
}

int main()
{
  constexpr int iterations = 100'000;
  constexpr int num_values = 100;

  std::println("=== Polynomial Evaluation Benchmark (Horner's Method) ===");
  std::println("    {} iterations × {} values\n", iterations, num_values);

  // Verify correctness
  std::println("--- Correctness Check ---");
  constexpr auto poly = make_test_polynomial<double, 20>();
  constexpr double at_one = poly(1.0);
  std::println("  poly exp(1) = {}", at_one);
  std::println("  std::exp(1) = {}", std::exp(1.0));
  std::println("  diff        = {}\n", std::abs(at_one - std::exp(1.0)));

  static_assert(poly(0.0) == 1.0);
  std::println("  ✓ Compile-time evaluation verified\n");

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
