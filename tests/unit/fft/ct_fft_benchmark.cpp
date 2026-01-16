/*
 *  ct_fft_benchmark.cpp
 *  Benchmarks default compile-time limits for FFT.
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Helper to ensure forcing compile-time evaluation
template<std::size_t N>
consteval bool test_fft_size()
{
  std::vector<std::complex<double>> data(N, {1.0, 0.0});
  auto res = fft::fft(data, false);
  return res.size() == N;
}

int main()
{
  // 256 is safe
  constexpr bool res256 = test_fft_size<256>();
  static_assert(res256);

  // 1024 is safe
  constexpr bool res1024 = test_fft_size<1024>();
  static_assert(res1024);

  // 4096 usually hits default steps limit, so we don't test it here to avoid CI failure
  // unless user boosted -fconstexpr-steps.

  if constexpr (res256 && res1024)
  {
    return 0; // pass
  }
  else
  {
    return 1; // fail
  }
}
