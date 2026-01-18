/*
 *  benchmark_operator_speed.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Benchmarks the latency/throughput of basic arithmetic operators (+, -, *)
 *    for small-degree polynomials to ensure minimal abstraction overhead.
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial;

template<std::size_t N>
void bench_operator(const char* label)
{
  using T = std::complex<double>;
  std::array<T, N + 1> coeffs;
  for (size_t i = 0; i <= N; ++i)
    coeffs[i] = {1.0, 1.0};

  polynomial_nttp<T, N> p(coeffs);
  polynomial_nttp<T, N> q(coeffs);

  // Warmup
  auto res_warmup = p * q;
  volatile auto sink1 = res_warmup[0];

  int iterations = 100;
  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    auto res = p * q;
    volatile auto sink = res[0];
  }
  auto end = std::chrono::steady_clock::now();

  double avg_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / (double)iterations;
  std::println("Degree N={:<4} | {:<20} | Time: {:>8.1f} us", N, label, avg_us);
}

int main()
{
  std::println("=== Operator* Performance Verification ===");
  std::println("(Should match FFT speeds for N >= 32)");

  bench_operator<16>("Small (Naive)");
  bench_operator<32>("Threshold (FFT?)");
  bench_operator<64>("Medium (FFT)");
  bench_operator<1024>("Large (FFT)");
  bench_operator<16384>("Sunspots (16k)");

  return 0;
}
