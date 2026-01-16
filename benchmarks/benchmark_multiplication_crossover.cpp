/*
 *  benchmark_multiplication_crossover.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Determines the crossover point (N) where FFT-based multiplication O(N log N)
 *    becomes faster than naive convolution O(N^2).
 */
#include <chrono>
#include <complex>
#include <print>
#include <random>
#include <vector>

import lam.polynomial_nttp;

// Import FFT module directly for the benchmark comparisons
// import :univariate.fft; // Removed invalid partition import

using namespace lam::polynomial;

// Naive O(N^2) Multiplication (simulating current operator*)
template<typename T>
std::vector<T> naive_multiply(const std::vector<T>& a, const std::vector<T>& b)
{
  std::size_t n = a.size();
  std::size_t m = b.size();
  std::vector<T> res(n + m - 1, 0);
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < m; ++j)
    {
      res[i + j] += a[i] * b[j];
    }
  }
  return res;
}

// FFT O(N log N) Multiplication
// 1. Pad to power of 2 >= sizeA + sizeB - 1
// 2. FFT(A), FFT(B)
// 3. Pointwise Multiply
// 4. IFFT(Result)
template<typename T>
std::vector<std::complex<double>> fft_multiply(const std::vector<T>& a, const std::vector<T>& b)
{
  std::size_t res_size = a.size() + b.size() - 1;
  std::size_t n = 1;
  while (n < res_size)
    n <<= 1;

  // Helper to convert/pad
  auto prepare = [&](const std::vector<T>& in) {
    std::vector<std::complex<double>> out(n, {0, 0});
    for (size_t i = 0; i < in.size(); ++i)
    {
      if constexpr (std::is_same_v<T, std::complex<double>>)
        out[i] = in[i];
      else
        out[i] = {in[i], 0};
    }
    return out;
  };

  auto fa = prepare(a);
  auto fb = prepare(b);

  // Call library FFT
  fa = univariate::fft::fft(std::move(fa), false);
  fb = univariate::fft::fft(std::move(fb), false);

  // Pointwise
  for (size_t i = 0; i < n; ++i)
    fa[i] *= fb[i];

  // Inverse
  fa = univariate::fft::fft(std::move(fa), true);

  // Trim
  fa.resize(res_size);
  return fa;
}

void run_benchmark(std::size_t N)
{
  using T = std::complex<double>;
  std::vector<T> a(N), b(N);
  for (size_t i = 0; i < N; ++i)
  {
    a[i] = {1.0, 1.0};
    b[i] = {1.0, -1.0};
  }

  // Warmup & Iteration Count Heuristic
  int iterations = 1000;
  if (N > 1000)
    iterations = 100;
  if (N > 5000)
    iterations = 10;

  auto start_naive = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    auto res = naive_multiply(a, b);
    volatile auto sink = res[0];
  }
  auto end_naive = std::chrono::high_resolution_clock::now();

  auto start_fft = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    auto res = fft_multiply(a, b);
    volatile auto sink = res[0];
  }
  auto end_fft = std::chrono::high_resolution_clock::now();

  auto dur_naive =
    std::chrono::duration_cast<std::chrono::microseconds>(end_naive - start_naive).count() / (double)iterations;
  auto dur_fft =
    std::chrono::duration_cast<std::chrono::microseconds>(end_fft - start_fft).count() / (double)iterations;

  std::print("N={:<4} | Naive: {:>7.1f} us | FFT: {:>7.1f} us | Winner: {}\n", N, dur_naive, dur_fft,
             (dur_fft < dur_naive ? "FFT" : "Naive"));
}

int main()
{
  std::print("=== Multiplication Crossover Benchmark ===\n");
  std::print("(Comparing Naive O(N^2) vs Accelerate FFT O(log N))\n\n");

  std::vector<size_t> sizes = {16, 32, 64, 80, 96, 128, 256, 512, 1024, 2048, 4096};

  for (auto n : sizes)
  {
    run_benchmark(n);
  }
  return 0;
}
