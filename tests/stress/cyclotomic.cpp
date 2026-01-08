/*
 *  cyclotomic.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Stress test for polynomial_nttp cyclotomic polynomials.
 */
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Cyclotomic Polynomial Phi_30(x)
// Degree phi(30) = 8.
// Roots are primitive 30th roots of unity. All complex.
// P(x) = x^8 + x^7 - x^5 - x^4 - x^3 + x + 1
int main()
{
  std::println("Stress Test: Cyclotomic Polynomials");
  std::println("===================================");

  using Real = double;

  // x^8 + x^7 - x^5 - x^4 - x^3 + x + 1
  constexpr polynomial_nttp<Real, 8> phi30{{1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 1.0}};

  std::println("Phi_30(x) = {}", phi30);

  // Solve with real arithmetic
  auto roots30 = lam::polynomial::univariate::roots::roots(phi30);

  std::println("Found {} roots:", roots30.size());
  for (auto r : roots30)
  {
    std::println("  {}", r.value);
  }

  // Expect 0 real roots
  if (roots30.empty())
  {
    std::println("SUCCESS: Found 0 real roots as expected.");
    return 0;
  }
  else
  {
    std::println("FAILED: Found phantom real roots.");
    return 1;
  }
}
