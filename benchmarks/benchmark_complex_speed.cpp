/*
 *  benchmark_complex_speed.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Micro-benchmark for complex polynomial evaluation, measuring potential speedups
 *    from Accelerate/BLAS optimized paths versus generic Horner's method.
 */
#include <chrono>
#include <complex>
#include <print>
#include <random>
#include <vector>

import lam.polynomial_nttp;

using namespace lam::polynomial;

template<std::size_t N>
void run_benchmark()
{
  univariate::polynomial_nttp<std::complex<double>, N> poly;

  // Deterministic pseudo-random fill
  for (std::size_t i = 0; i <= N; ++i)
  {
    poly.coefficients[i] = {static_cast<double>(i % 10), static_cast<double>((i * 2) % 7)};
  }

  std::complex<double> x{0.99, 0.01};
  constexpr int iterations = 50000;

  // Warmup
  std::complex<double> sink = poly(x);

  // 1. Generic Horner
  auto start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    sink += poly.evaluate_horner(x);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto horner_dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  // 2. Accelerate (Optimized)
  start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    sink += poly.evaluate_accelerate(x);
  }
  end = std::chrono::high_resolution_clock::now();
  auto accel_dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  // Metrics Calculation
  // Complex FMA approx: 4 muls + 4 adds = 8 FLOPs per degree
  double total_flops = 8.0 * N * iterations;
  double gflops = (total_flops / (accel_dur / 1e6)) / 1e9;

  // Bandwidth: We read N+1 complex coefficients (16 bytes) per iteration + x (16 bytes)
  // All in L1 cache likely, but still measures "throughput" of data
  double total_bytes = (double)(N + 1) * 16.0 * iterations;
  double gb_s = (total_bytes / (accel_dur / 1e6)) / 1e9;

  std::print("Degree N = {:<4} | Horner: {:>7} us | Accelerate: {:>7} us | Speedup: {:.2f}x | Performance: {:.2f} "
             "GFLOPS, {:.2f} GB/s (sink={})\n",
             N, horner_dur, accel_dur, (double)horner_dur / accel_dur, gflops, gb_s, sink.real());
}

int main()
{
  std::print("=== Complex Polynomial Evaluation Benchmark ===\n");
  std::print("Iterations: 50,000\n\n");

  run_benchmark<100>();
  run_benchmark<500>();
  run_benchmark<1000>();
  run_benchmark<2000>();
  run_benchmark<5000>();

  return 0;
}
