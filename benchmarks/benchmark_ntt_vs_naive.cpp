/*
 *  benchmark_ntt_vs_naive.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Use a prime that supports large NTTs
// P = c * 2^k + 1
// 0xFFFFFFFF00000001 = 2^64 - 2^32 + 1 (Solinas prime) - NO, that's not necessarily NTT friendly for all sizes
// Let's use a specific NTT prime.
// 4179340454199820289 = 29 * 2^57 + 1
constexpr auto ntt_prime = 4179340454199820289_Z;
using field = decltype(lam::cbn::Zq(ntt_prime));

template<std::size_t N>
void run_benchmark()
{
  using poly = lam::polynomial_nttp<field, N>;

  // Generate deterministic random-like polynomials
  poly p, q;
  for (std::size_t i = 0; i <= N; ++i)
  {
    // ZqElement constructible from long
    p.coefficients[i] = field(static_cast<long>((i * 12345) ^ 0xDEADBEEF));
    q.coefficients[i] = field(static_cast<long>((i * 67890) ^ 0xCAFEBABE));
  }

  auto start = std::chrono::steady_clock::now();

  constexpr int iterations = 1000;
  // Volatile to prevent optimization
  volatile int dummy = 0;

  for (int i = 0; i < iterations; ++i)
  {
    auto res = p * q;
    // Access via .data (big_int) -> [0] (limb)
    dummy += (int)(res[0].data[0]);
  }

  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  double avg_us = static_cast<double>(duration) / iterations;

  std::println("Degree N={}: {:.3f} us", N, avg_us);
}

int main()
{
  std::println("Benchmarking Polynomial Multiplication (NTT vs Naive)");
  std::println("Threshold is likely around N=64");

  // Below Threshold (Naive)
  run_benchmark<32>();

  // Around Threshold
  run_benchmark<64>();
  run_benchmark<128>();

  // Above Threshold (NTT)
  run_benchmark<256>();
  run_benchmark<512>();
  run_benchmark<1024>();

  return 0;
}
