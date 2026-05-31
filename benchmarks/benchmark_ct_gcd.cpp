// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: 2025-2026 Colin Ford

/*
 *  benchmark_ct_gcd.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */

import std;
import lam.polynomial.nttp;

constexpr std::size_t N = 256;

template<std::size_t Size>
consteval auto make_poly_pattern()
{
  lam::polynomial::nttp::univariate::polynomial_nttp<double, Size> p{};
  for (std::size_t i = 0; i <= Size; ++i)
  {
    p.coefficients[i] = static_cast<double>((i % 17) + 1);
  }
  return p;
}

consteval auto gcd_bench()
{
  constexpr auto a = make_poly_pattern<N>();
  auto b = make_poly_pattern<N>();
  b.coefficients[0] += 1.0;

  auto res = lam::polynomial::nttp::poly_gcd(a, b);
  return res;
}

int main()
{
  constexpr auto res = gcd_bench();
  std::println("GCD Result (N={}): {}", N, res);
  return 0;
}
