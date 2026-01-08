
/*
 *  mignotte.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Stress test for polynomial_nttp Mignotte polynomials.
 */
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Mignotte Polynomial: P(x) = x^n - 2(ax - 1)^2
// With n=5, a=100.
// P(x) = x^5 - 2(10000x^2 - 200x + 1)
//      = x^5 - 20000x^2 + 400x - 2
int main()
{
  std::println("Stress Test: Mignotte Polynomials");
  std::println("=================================");

  // n=5, a=100
  using Real = double;
  constexpr polynomial_nttp<Real, 5> m5{{-2.0, 400.0, -20000.0, 0.0, 0.0, 1.0}};

  std::println("Mignotte(5, 100): P(x) = {}", m5);

  auto roots5 = lam::polynomial::univariate::roots::roots(m5);
  std::println("Found {} roots:", roots5.size());
  std::vector<Real> close_roots;

  for (auto r : roots5)
  {
    std::println("  {:.16f}", r.value);
    if (std::abs(r.value - 0.01) < 1e-4)
    {
      close_roots.push_back(r.value);
    }
  }

  // Check for cluster near 1/a = 0.01
  // Expected: 2 real roots very close to 0.01
  if (close_roots.size() >= 2)
  {
    std::println("SUCCESS: Found cluster of {} roots near 0.01", close_roots.size());
    double sep = std::abs(close_roots[0] - close_roots[1]);
    std::println("Separation: {:.6e}", sep);
    // Theoretical separation approx sqrt(2) * a^(-n/2 - 1)
    // n=5, a=100 -> 1.41 * 100^(-3.5) = 1.41 * 10^-7
    if (sep < 1e-5 && sep > 1e-9)
    {
      std::println("Separation consistent with theory.");
      return 0;
    }
  }
  else
  {
    std::println("FAILED: Did not find root cluster near 0.01");
    return 1;
  }

  return 0;
}
