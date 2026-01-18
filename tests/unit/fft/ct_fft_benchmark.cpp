/*
 *  ct_fft_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  an FFT benchmark for compile-time evaluation
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
  constexpr bool res256 = test_fft_size<256>();
  static_assert(res256);

  constexpr bool res1024 = test_fft_size<1024>();
  static_assert(res1024);
  // 4096 usually hits default steps limit
  if constexpr (res256 && res1024)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}
