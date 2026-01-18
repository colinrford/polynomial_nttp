/*
 *  benchmark_ct_naive.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */
import std;
import lam.polynomial_nttp;

using complex = std::complex<double>;
constexpr std::size_t N = 256;

template<std::size_t Size>
consteval auto naive_mult_bench()
{
  std::array<complex, Size> a, b;
  for (std::size_t i = 0; i < Size; ++i)
  {
    a[i] = {1.0, 0.0};
    b[i] = {1.0, 0.0};
  }

  std::array<complex, 2 * Size> res{};
  for (std::size_t i = 0; i < Size; ++i)
  {
    for (std::size_t j = 0; j < Size; ++j)
    {
      res[i + j] += a[i] * b[j];
    }
  }
  return res[0];
}

int main()
{
  constexpr auto res = naive_mult_bench<N>();
  std::println("Naive Result (N={}): {}", N, res.real());
  return 0;
}
