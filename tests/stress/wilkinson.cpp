/*
 *  wilkinson.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Stress test for polynomial_nttp Wilkinson polynomials.
 */
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Helper to expand (x - 1)(x - 2)...(x - N)
// This is done at runtime for simplicity, or we can use small N constants
template<typename T, std::size_t N>
constexpr auto make_wilkinson()
{
  // Start with P(x) = 1
  polynomial_nttp<T, 0> p{{T(1)}};

  auto expand_recursive = [&]<std::size_t I>(auto&& self, std::integral_constant<std::size_t, I>) {
    if constexpr (I > N)
    {
      return p;
    }
    else
    {
      // Multiply by (x - I)
      // Current p has degree I-1. New p has degree I.
      // But types are static. This recursion approach is tricky with static types.
      // Better approach: hardcode or use a fixed max size container?
      // Actually, simply constructing coefficients from roots is easier if we know them.
      // But computing elementary symmetric polynomials for N=20 is huge.
      return self(self, std::integral_constant<std::size_t, I + 1>{});
    }
  };
  // ... complicated.

  // Easier approach: Just define the coefficients for small N, or build it dynamically
  // but the library requires static extent.
  // For N=10, we could precompute.
}

// Let's manually define Wilkinson 5, 10 for now using roots expansion logic if possible.
// Or just test specific known cases.

int main()
{
  std::println("Stress Test: Wilkinson Polynomials");
  std::println("==================================");

  // Wilkinson 5: (x-1)(x-2)(x-3)(x-4)(x-5)
  // = x^5 - 15x^4 + 85x^3 - 225x^2 + 274x - 120
  using Real = double;

  constexpr polynomial_nttp<Real, 5> w5{{-120.0, 274.0, -225.0, 85.0, -15.0, 1.0}};

  std::println("\nWilkinson N=5:");
  std::println("P(x) = {}", w5);

  auto roots5 = lam::polynomial::univariate::roots::roots(w5);
  std::println("Found {} roots:", roots5.size());
  for (auto r : roots5)
  {
    std::println("  {}", r.value);
  }

  // Check if close to integer roots 1..5
  bool pass5 = true;
  for (int i = 1; i <= 5; ++i)
  {
    bool found = false;
    for (auto r : roots5)
    {
      if (std::abs(r.value - i) < 1e-4)
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      std::println("FAILED to find root {}", i);
      pass5 = false;
    }
  }

  if (pass5)
    std::println("Wilkinson N=5 PASSED");
  else
    std::println("Wilkinson N=5 FAILED");

  return pass5 ? 0 : 1;
}

// Helper to check approximate roots
template<typename RootsContainer>
constexpr bool check_roots_contain(const RootsContainer& roots, std::initializer_list<double> expected)
{
    if (roots.size() != expected.size()) return false;
    for (double exp : expected) {
        bool found = false;
        for (auto r : roots) {
            if (r.value > exp - 0.1 && r.value < exp + 0.1) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

// Compile-time verification
// N=3: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
constexpr bool test_w3_constexpr() {
    using Real = double;
    constexpr polynomial_nttp<Real, 3> w3{{-6.0, 11.0, -6.0, 1.0}};
    auto r = lam::polynomial::univariate::roots::roots(w3);
    return check_roots_contain(r, {1.0, 2.0, 3.0});
}
static_assert(test_w3_constexpr(), "Wilkinson(3) [Cardano] should be solvable at compile-time");

// N=4: (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
constexpr bool test_w4_constexpr() {
    using Real = double;
    constexpr polynomial_nttp<Real, 4> w4{{24.0, -50.0, 35.0, -10.0, 1.0}};
    auto r = lam::polynomial::univariate::roots::roots(w4);
    return check_roots_contain(r, {1.0, 2.0, 3.0, 4.0});
}
static_assert(test_w4_constexpr(), "Wilkinson(4) [Newton] should be solvable at compile-time");

// N=5: (x-1)...(x-5)
constexpr bool test_w5_constexpr() {
    using Real = double;
    constexpr polynomial_nttp<Real, 5> w5{{-120.0, 274.0, -225.0, 85.0, -15.0, 1.0}};
    auto r = lam::polynomial::univariate::roots::roots(w5);
    return check_roots_contain(r, {1.0, 2.0, 3.0, 4.0, 5.0});
}
static_assert(test_w5_constexpr(), "Wilkinson(5) [Newton Rec.] should be solvable at compile-time");

