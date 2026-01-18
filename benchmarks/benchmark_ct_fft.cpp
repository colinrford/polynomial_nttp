/*
 *  benchmark_ct_fft.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate::fft;
using complex = std::complex<double>;
constexpr std::size_t N = 256;

template<std::size_t Size>
consteval auto fft_mult_bench()
{
  std::size_t sz = 1;
  while (sz < 2 * Size)
    sz <<= 1;

  std::vector<complex> fa(sz), fb(sz);
  for (std::size_t i = 0; i < Size; ++i)
  {
    fa[i] = {1.0, 0.0};
    fb[i] = {1.0, 0.0};
  }

  fa = fft(std::move(fa), false);
  fb = fft(std::move(fb), false);

  for (std::size_t i = 0; i < sz; ++i)
    fa[i] *= fb[i];

  fa = fft(std::move(fa), true);

  return fa[0];
}

int main()
{
  constexpr auto res = fft_mult_bench<N>();
  std::println("FFT Result (N={}): {}", N, res.real());
  return 0;
}
