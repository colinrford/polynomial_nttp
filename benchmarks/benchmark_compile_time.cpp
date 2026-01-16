/*
 *  benchmark_compile_time.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Benchmarks compilation time and compile-time evaluation limits for:
 *      - FFT Multiplication (vs Naive)
 *      - Large Compile-Time Datasets (Sunspots)
 *      - Polynomial Division and GCD
 */
#include <array>
#include <complex>
#include <vector>

import lam.polynomial_nttp;

// We need to access the internal implementations if we want to isolate them
// But the module exports them mostly via the main types.
// For this benchmark, we'll reimplement the Naive/FFT constexpr logic LOCALLY
// to ensure we are purely measuring the algorithm cost at compile time
// without fighting module visibility of implementation details.

using Complex = std::complex<double>;

// --- 1. Naive Constexpr Multiplication ---
template<std::size_t N>
consteval auto naive_mult_bench()
{
  std::array<Complex, N> a, b;
  for (size_t i = 0; i < N; ++i)
  {
    a[i] = {1, 0};
    b[i] = {1, 0};
  }

  std::array<Complex, 2 * N> res{};
  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < N; ++j)
    {
      res[i + j] += a[i] * b[j];
    }
  }
  return res[0]; // Return something to prevent DCE
}

// --- 2. FFT Constexpr Multiplication ---
// We use the library's implementation because it has valid constexpr sin/cos support
using namespace lam::polynomial::univariate::fft;

template<std::size_t N>
consteval auto fft_mult_bench()
{
  // Pad to power of 2
  size_t sz = 1;
  while (sz < 2 * N)
    sz <<= 1;

  std::vector<Complex> fa(sz), fb(sz);
  for (size_t i = 0; i < N; ++i)
  {
    fa[i] = {1, 0};
    fb[i] = {1, 0};
  }

  // Call Library FFT (constexpr)
  fa = fft(std::move(fa), false);
  fb = fft(std::move(fb), false);

  for (size_t i = 0; i < sz; ++i)
    fa[i] *= fb[i];

  fa = fft(std::move(fa), true);

  return fa[0];
}

#if BENCH_MODE >= 2
// Sunspots Benchmark Data
#include "../examples/data/sunspot_data.hpp"

consteval auto sunspot_bench()
{
  constexpr std::size_t N = lam::examples::SunspotCount;
  std::vector<Complex> fa(N);
  for (size_t i = 0; i < N; ++i)
  {
    fa[i] = {lam::examples::sunspot_data[i], 0};
  }

  fa = fft(std::move(fa), false);
  for (size_t i = 0; i < N; ++i)
    fa[i] *= fa[i];
  fa = fft(std::move(fa), true);

  return fa[0];
}

consteval auto sunspot_bench_32k()
{
  constexpr std::size_t N = 32768; // 32k
  std::vector<Complex> fa(N);
  // Fill with data repeating
  for (size_t i = 0; i < N; ++i)
  {
    fa[i] = {lam::examples::sunspot_data[i % lam::examples::SunspotCount], 0};
  }

  fa = fft(std::move(fa), false);
  for (size_t i = 0; i < N; ++i)
    fa[i] *= fa[i];
  fa = fft(std::move(fa), true);

  return fa[0];
}
#endif

#if BENCH_MODE >= 4
// --- 4. Division (NTTP) & 5. GCD (Function) ---

// Helper to generate a deterministic polynomial at compile time
template<std::size_t N>
consteval auto make_poly_pattern() {
    lam::polynomial::univariate::polynomial_nttp<double, N> p{};
    for(size_t i=0; i<=N; ++i) {
        p.coefficients[i] = (double)((i % 17) + 1); // Avoid zeros
    }
    return p;
}

#if BENCH_MODE == 4
// Division uses NTTPs, so operands must be static constexpr or similar validity scope
constexpr auto div_a = make_poly_pattern<BENCH_N>();
// Create a divisor of smaller degree approx N/2
constexpr auto div_b = make_poly_pattern<BENCH_N/2>();

consteval auto division_bench() {
    // division_prototype returns pair<quotient, remainder>
    constexpr auto res = lam::division_prototype<double, BENCH_N, div_a, BENCH_N/2, div_b>();
    // Access result to ensure evaluation
    return res.second[0]; // First coeff of remainder
}
#endif

#if BENCH_MODE == 5
consteval auto gcd_bench() {
    constexpr size_t N = BENCH_N;
    constexpr auto a = make_poly_pattern<N>();
    // Make b slightly different to avoid trivial GCD
    auto b = make_poly_pattern<N>(); 
    b.coefficients[0] += 1.0; 
    
    // poly_gcd returns polynomial_nttp
    auto res = lam::polynomial::poly_gcd(a, b);
    return res[0];
}
#endif

#endif

#include <iostream>

int main()
{
  // BENCH_MODE defined via -D compiler flag
  // 0 = Naive, 1 = FFT, 2 = Sunspots (FFT)

#if BENCH_MODE == 0
  constexpr auto res = naive_mult_bench<BENCH_N>();
  std::cout << "Naive Result (N=" << BENCH_N << "): " << res.real() << std::endl;
#elif BENCH_MODE == 1
  constexpr auto res = fft_mult_bench<BENCH_N>();
  std::cout << "FFT Result (N=" << BENCH_N << "): " << res.real() << std::endl;
#elif BENCH_MODE == 2
  constexpr auto res = sunspot_bench();
  std::cout << "FFT Result (Compile-Time): " << res.real() << std::endl;
  // ... verification verification ...
  double exact = 0;
  constexpr auto& data = lam::examples::sunspot_data;
  constexpr size_t N = 16384;
  exact += data[0] * data[0];
  for(size_t k=1; k<N; ++k) exact += data[k] * data[N-k];
  std::cout << "Exact Result (Runtime):    " << exact << std::endl;
  std::cout << "Difference:                " << (res.real() - exact) << std::endl;

#elif BENCH_MODE == 3
  constexpr auto res = sunspot_bench_32k();
  std::cout << "FFT Result 32k (Compile-Time): " << res.real() << std::endl;
#endif

#if BENCH_MODE == 4
// --- 4. Division (NTTP) ---
    constexpr auto res = division_bench();
    std::cout << "Division Result (N=" << BENCH_N << "): " << res << std::endl;
#elif BENCH_MODE == 5
// --- 5. GCD (Function) ---
    constexpr auto res = gcd_bench();
    std::cout << "GCD Result (N=" << BENCH_N << "): " << res << std::endl;
#endif
  return 0;
}
// force rebuild naive
