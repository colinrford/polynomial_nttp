/*
 *  benchmark_ct_sunspots.cpp
 *    Benchmarks compilation time for Large Dataset FFT (Sunspots N=16384)
 */
import std;
import lam.polynomial_nttp;

#include "../examples/data/sunspot_data.hpp"

using namespace lam::polynomial::univariate::fft;
using complex = std::complex<double>;

consteval auto sunspot_bench()
{
  constexpr std::size_t N = lam::examples::SunspotCount;
  std::vector<complex> fa(N);
  for (std::size_t i = 0; i < N; ++i)
  {
    fa[i] = {lam::examples::sunspot_data[i], 0};
  }

  fa = fft(std::move(fa), false);
  for (std::size_t i = 0; i < N; ++i)
    fa[i] *= fa[i];
  fa = fft(std::move(fa), true);

  return fa[0];
}

int main()
{
  constexpr auto res = sunspot_bench();
  std::println("FFT Result (Compile-Time): {}", res.real());

  double exact = 0;
  constexpr auto& data = lam::examples::sunspot_data;
  constexpr std::size_t N = 16384;
  exact += data[0] * data[0];
  for (std::size_t k = 1; k < N; ++k)
    exact += data[k] * data[N - k];
    
  std::println("Exact Result (Runtime):    {}", exact);
  std::println("Difference:                {}", (res.real() - exact));
  return 0;
}
