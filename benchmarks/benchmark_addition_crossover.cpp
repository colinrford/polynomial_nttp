/*
 *  benchmark_addition_crossover.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Determines the crossover point (N) where vDSP accelerated addition
 *    becomes faster than naive loop addition.
 */
#include <algorithm>
#include <chrono>
#include <complex>
#include <print>
#include <vector>

// Standalone vDSP declarations for benchmarking without module internal access
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

template<typename T>
void naive_add(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& res)
{
  std::size_t n = a.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    res[i] = a[i] + b[i];
  }
}

template<typename T>
void accelerated_add(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& res)
{
#ifdef __APPLE__
  if constexpr (std::is_same_v<T, double>)
  {
    vDSP_vaddD(a.data(), 1, b.data(), 1, res.data(), 1, a.size());
  }
  else if constexpr (std::is_same_v<T, float>)
  {
    vDSP_vadd(a.data(), 1, b.data(), 1, res.data(), 1, a.size());
  }
  else if constexpr (std::is_same_v<T, std::complex<double>>)
  {
    vDSP_vaddD(reinterpret_cast<const double*>(a.data()), 1, reinterpret_cast<const double*>(b.data()), 1,
               reinterpret_cast<double*>(res.data()), 1, a.size() * 2);
  }
#else
  // Fallback for non-Apple (just to compile)
  naive_add(a, b, res);
#endif
}

void run_benchmark(std::size_t N)
{
  using T = double;
  std::vector<T> a(N, 1.0), b(N, 2.0), res(N);

  int iterations = 100000;
  if (N > 1000)
    iterations = 10000;
  if (N > 100000)
    iterations = 1000;

  auto start_naive = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    naive_add(a, b, res);
    volatile auto sink = res[0];
  }
  auto end_naive = std::chrono::steady_clock::now();

  auto start_acc = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    accelerated_add(a, b, res);
    volatile auto sink = res[0];
  }
  auto end_acc = std::chrono::steady_clock::now();

  auto dur_naive =
    std::chrono::duration_cast<std::chrono::nanoseconds>(end_naive - start_naive).count() / (double)iterations;
  auto dur_acc = std::chrono::duration_cast<std::chrono::nanoseconds>(end_acc - start_acc).count() / (double)iterations;

  std::print("N={:<4} | Naive: {:>7.1f} ns | vDSP: {:>7.1f} ns | Winner: {}\n", N, dur_naive, dur_acc,
             (dur_acc < dur_naive ? "vDSP" : "Naive"));
}

int main()
{
  std::print("=== Addition Crossover Benchmark ===\n");
  std::vector<size_t> sizes = {16, 32, 64, 80, 96, 128, 256, 512, 1024, 4096, 16384};

  for (auto n : sizes)
  {
    run_benchmark(n);
  }
  return 0;
}
