/*
 *  benchmark_ct_div.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */
import std;
import lam.polynomial_nttp;

constexpr std::size_t N = 1024;

template<std::size_t Size>
consteval auto make_poly_pattern()
{
  lam::polynomial::univariate::polynomial_nttp<double, Size> p{};
  for (std::size_t i = 0; i <= Size; ++i)
  {
    p.coefficients[i] = (double)((i % 17) + 1);
  }
  return p;
}

constexpr auto div_a = make_poly_pattern<N>();
constexpr auto div_b = make_poly_pattern<N / 2>();

consteval auto division_bench()
{
  constexpr auto res = lam::division_prototype<double, N, div_a, N / 2, div_b>();
  return res.second[0];
}

int main()
{
  constexpr auto res = division_bench();
  std::println("Division Result (N={}): {}", N, res);
  return 0;
}
