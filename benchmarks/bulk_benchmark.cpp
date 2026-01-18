/*
 *  bulk_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Benchmarks high-throughput evaluation of polynomials on large datasets using:
 *      - TBB (std::execution::par_unseq equivalent)
 *      - SIMD Backends (Accelerate/BLAS)
 *    Comparison against Boost.Math parallel execution.
 */
#include <algorithm>
#include <chrono>
#include <print>
#include <random>
#include <vector>

#ifdef HAS_BOOST_MATH
#include <boost/math/tools/polynomial.hpp>
#endif

// Guard TBB include if not already done by module?
// No, module does internal things. We need it here for the MANUAL boost loop.
#ifdef LAM_USE_TBB
#include <tbb/parallel_for.h>
#endif

import lam.polynomial_nttp;

// Generate random points
std::vector<double> generate_random_points(std::size_t n)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-2.0, 2.0);
  std::vector<double> v(n);
  std::generate(v.begin(), v.end(), [&]() { return dis(gen); });
  return v;
}

// Benchmark runner for a specific degree
template<std::size_t N>
void run_benchmark(const std::vector<double>& inputs, std::vector<double>& outputs)
{
  using R = double;
  lam::polynomial_nttp<R, N> poly;
  for (int i = 0; i <= N; ++i)
    poly.coefficients[i] = 1.0;

  // Warmup
  poly.evaluate_bulk(std::span{inputs}.subspan(0, 100), std::span{outputs}.subspan(0, 100));

  auto start = std::chrono::steady_clock::now();
  poly.evaluate_bulk(inputs, outputs);
  auto end = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  double mps = static_cast<double>(inputs.size()) / (static_cast<double>(duration.count()) / 1000.0);

  std::println("  lam::poly<{:3}>: {:4} ms | Throughput: {:6.2f} M/sec", N, duration.count(), mps / 1e6);
}

#ifdef HAS_BOOST_MATH
template<std::size_t N>
void run_boost_benchmark(const std::vector<double>& inputs, std::vector<double>& outputs)
{
  using R = double;
  // Create Boost polynomial (all 1.0)
  std::vector<R> coeffs(N + 1, 1.0);
  boost::math::tools::polynomial<R> poly(coeffs.data(), N);

  auto start = std::chrono::steady_clock::now();

  if constexpr (lam::polynomial::config::use_tbb)
  {
#ifdef LAM_USE_TBB
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, inputs.size()), [&](const tbb::blocked_range<std::size_t>& r) {
      for (std::size_t i = r.begin(); i != r.end(); ++i)
      {
        outputs[i] = poly.evaluate(inputs[i]);
      }
    });
#endif
  }
  else
  {
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
      outputs[i] = poly.evaluate(inputs[i]);
    }
  }

  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  double mps = static_cast<double>(inputs.size()) / (static_cast<double>(duration.count()) / 1000.0);

  std::println("  Boost<{:3}>    : {:4} ms | Throughput: {:6.2f} M/sec", N, duration.count(), mps / 1e6);
}
#endif

template<std::size_t N>
void run_complex_benchmark(std::size_t M)
{
  using R = std::complex<double>;
  std::vector<R> inputs(M);
  std::vector<R> outputs(M);

  // Fill inputs
  for (size_t i = 0; i < M; ++i)
    inputs[i] = {static_cast<double>(i % 100) * 0.1, static_cast<double>(i % 100) * 0.1};

  // Create poly with random complex coeffs
  std::array<R, N + 1> coeffs;
  for (auto& c : coeffs)
    c = {1.0, 1.0};
  lam::polynomial::polynomial_nttp<R, N> poly(coeffs);

  auto start = std::chrono::steady_clock::now();
  poly.evaluate_bulk(inputs, outputs);
  auto end = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  double mps = static_cast<double>(M) / (static_cast<double>(duration.count()) / 1000.0);

  std::println("  Cplx<{:3}>     : {:4} ms | Throughput: {:6.2f} M/sec", N, duration.count(), mps / 1e6);
}

template<std::size_t N>
void run_complex_scalar_benchmark(std::size_t M)
{
  using R = std::complex<double>;
  std::vector<R> inputs(M);
  std::vector<R> outputs(M);
  for (size_t i = 0; i < M; ++i)
    inputs[i] = {static_cast<double>(i % 100) * 0.1, static_cast<double>(i % 100) * 0.1};
  std::array<R, N + 1> coeffs;
  for (auto& c : coeffs)
    c = {1.0, 1.0};
  lam::polynomial::polynomial_nttp<R, N> poly(coeffs);

  // Force scalar path by calling evaluate_horner manually inside a loop
  auto start = std::chrono::steady_clock::now();
  for (size_t i = 0; i < M; ++i)
  {
    outputs[i] = poly.evaluate_horner(inputs[i]);
  }
  auto end = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  double mps = static_cast<double>(M) / (static_cast<double>(duration.count()) / 1000.0);

  std::println("  CplxScalar<{:3}> : {:4} ms | Throughput: {:6.2f} M/sec", N, duration.count(), mps / 1e6);
}

#ifdef HAS_BOOST_MATH
template<std::size_t N>
void run_boost_complex_benchmark(std::size_t M)
{
  using R = std::complex<double>;
  std::vector<R> inputs(M);
  std::vector<R> outputs(M);
  for (size_t i = 0; i < M; ++i)
    inputs[i] = {static_cast<double>(i % 100) * 0.1, static_cast<double>(i % 100) * 0.1};

  std::vector<R> coeffs(N + 1);
  for (auto& c : coeffs)
    c = {1.0, 1.0};
  boost::math::tools::polynomial<R> poly(coeffs.data(), N);

  auto start = std::chrono::steady_clock::now();
  // TBB Parallel wrapper for Boost
#ifdef LAM_USE_TBB
  if constexpr (lam::polynomial::config::use_tbb)
  {
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, M), [&](const tbb::blocked_range<std::size_t>& r) {
      for (std::size_t i = r.begin(); i != r.end(); ++i)
      {
        outputs[i] = poly.evaluate(inputs[i]);
      }
    });
  }
  else
  {
    for (size_t i = 0; i < M; ++i)
      outputs[i] = poly.evaluate(inputs[i]);
  }
#else
  for (size_t i = 0; i < M; ++i)
    outputs[i] = poly.evaluate(inputs[i]);
#endif
  auto end = std::chrono::steady_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  double mps = static_cast<double>(M) / (static_cast<double>(duration.count()) / 1000.0);

  std::println("  BoostCplx<{:3}>  : {:4} ms | Throughput: {:6.2f} M/sec", N, duration.count(), mps / 1e6);
}
#endif

int main()
{
  constexpr std::size_t M = lam::polynomial::config::is_tsan_build ? 100'000 : 10'000'000;

  // Warmup
  std::println("Warming up...");
  // Check Configuration
  std::print("=== Bulk Evaluation Benchmark (M={}) [", M);
  if constexpr (lam::polynomial::config::use_tbb)
    std::print("TBB:ON ");
  else
    std::print("TBB:OFF ");
  if constexpr (lam::polynomial::config::use_accelerate)
    std::print("ACC:ON ");
  else
    std::print("ACC:OFF ");
  if constexpr (lam::polynomial::config::use_blas)
    std::print("BLAS:ON ");
  else
    std::print("BLAS:OFF ");
  std::println("] ===");

  std::vector<double> inputs = generate_random_points(M);
  std::vector<double> outputs(M);

  run_benchmark<10>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<10>(inputs, outputs);
#endif
  run_benchmark<15>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<15>(inputs, outputs);
#endif
  run_benchmark<20>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<20>(inputs, outputs);
#endif
  run_benchmark<30>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<30>(inputs, outputs);
#endif
  run_benchmark<40>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<40>(inputs, outputs);
#endif
  run_benchmark<50>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<50>(inputs, outputs);
#endif
  run_benchmark<100>(inputs, outputs);
#ifdef HAS_BOOST_MATH
  run_boost_benchmark<100>(inputs, outputs);
#endif
  run_benchmark<500>(inputs, outputs);

  std::println("\n=== Complex Evaluation Benchmark ===");
  run_complex_benchmark<10>(M);
#ifdef HAS_BOOST_MATH
  run_boost_complex_benchmark<10>(M);
#endif

  // N=80: Check crossover
  run_complex_benchmark<80>(M);
  run_complex_scalar_benchmark<80>(M);

  // N=100: Check acceleration
  run_complex_benchmark<100>(M);
  run_complex_scalar_benchmark<100>(M);
#ifdef HAS_BOOST_MATH
  run_boost_complex_benchmark<100>(M);
#endif

  return 0;
}
