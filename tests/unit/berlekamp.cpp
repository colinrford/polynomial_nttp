/*
 *  berlekamp.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp berlekamp_factor function
 */
 
import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::polynomial::univariate;
using namespace lam::polynomial::univariate::berlekamp;
using namespace lam::cbn;
using namespace lam::cbn::literals;

// ============================================================
// Test: Polynomial GCD over GF(5)
// ============================================================

bool test_poly_gcd()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), four(4_Z);

  // gcd(x^2 - 1, x - 1) = x - 1
  // x^2 - 1 = x^2 + 4 in GF(5)
  polynomial_nttp<GF, 2> a{{four, zero, one}}; // x^2 + 4
  polynomial_nttp<GF, 2> b{{four, one, zero}}; // x + 4

  auto g = poly_gcd(a, b);

  // Result should be monic (x - 1), i.e., x + 4 in GF(5)
  // Check: g[1] should be 1, g[0] should be 4
  if (g[1] != one)
  {
    std::println("GCD: expected g[1]=1, got something else");
    return false;
  }
  if (g[0] != four)
  {
    std::println("GCD: expected g[0]=4, got something else");
    return false;
  }

  std::println("poly_gcd: gcd(x^2-1, x-1) = x-1 ✓");
  return true;
}

// ============================================================
// Test: Polynomial GCD - coprime case
// ============================================================

bool test_poly_gcd_coprime()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), two(2_Z);

  // gcd(x + 1, x + 2) should be 1 (coprime)
  polynomial_nttp<GF, 2> a{{one, one, zero}}; // x + 1
  polynomial_nttp<GF, 2> b{{two, one, zero}}; // x + 2

  auto g = poly_gcd(a, b);

  // Result should be constant 1
  // Since both are linear and distinct, gcd is 1
  // Monic constant = 1
  if (g[0] != one)
  {
    std::println("Coprime GCD: expected g[0]=1");
    return false;
  }
  // Higher coefficients should be zero
  if (g[1] != zero || g[2] != zero)
  {
    std::println("Coprime GCD: expected higher coeffs to be 0");
    return false;
  }

  std::println("poly_gcd: gcd(x+1, x+2) = 1 ✓");
  return true;
}

// ============================================================
// Test: GCD with common factor
// ============================================================

bool test_poly_gcd_common_factor()
{
  constexpr auto mod = 7_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), two(2_Z), three(3_Z), six(6_Z);

  // a = (x - 1)(x - 2) = x^2 - 3x + 2 = x^2 + 4x + 2 in GF(7)
  // b = (x - 1)(x - 3) = x^2 - 4x + 3 = x^2 + 3x + 3 in GF(7)
  // gcd = (x - 1) = x + 6 in GF(7)

  polynomial_nttp<GF, 2> a{{two, GF(4_Z), one}}; // x^2 + 4x + 2
  polynomial_nttp<GF, 2> b{{three, three, one}}; // x^2 + 3x + 3

  auto g = poly_gcd(a, b);

  // Should be monic linear: x + 6
  if (g[1] != one)
  {
    std::println("Common factor GCD: expected g[1]=1");
    return false;
  }
  if (g[0] != six)
  {
    std::println("Common factor GCD: expected g[0]=6");
    return false;
  }

  std::println("poly_gcd: gcd((x-1)(x-2), (x-1)(x-3)) = x-1 ✓");
  return true;
}

// ============================================================
// Test: Berlekamp Factorization
// ============================================================

bool test_berlekamp_factor()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), four(4_Z);

  // Factor x^2 - 1 = (x-1)(x+1) = (x+4)(x+1) over GF(5)
  polynomial_nttp<GF, 2> f{{four, zero, one}}; // x^2 + 4 = x^2 - 1

  // P = 5 (field size), pass zero and one for element construction
  auto [factors, count] = berlekamp_factor<GF, 5, 2>(f, zero, one);

  if (count != 2)
  {
    std::println("berlekamp_factor: expected 2 factors, got {}", count);
    return false;
  }

  std::println("berlekamp_factor: x^2-1 factors into {} irreducibles ✓", count);
  return true;
}

// ============================================================
// Compile-Time Tests
// ============================================================

consteval bool test_poly_gcd_constexpr()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), four(4_Z);

  // gcd(x^2 - 1, x - 1) = x - 1
  polynomial_nttp<GF, 2> a{{four, zero, one}};
  polynomial_nttp<GF, 2> b{{four, one, zero}};

  auto g = poly_gcd(a, b);

  return g[1] == one && g[0] == four;
}

consteval bool test_poly_gcd_coprime_constexpr()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), two(2_Z);

  // gcd(x + 1, x + 2) should be 1 (coprime)
  polynomial_nttp<GF, 2> a{{one, one, zero}}; // x + 1
  polynomial_nttp<GF, 2> b{{two, one, zero}}; // x + 2

  auto g = poly_gcd(a, b);

  // Result should be constant 1, higher coeffs 0
  return g[0] == one && g[1] == zero && g[2] == zero;
}

consteval bool test_poly_gcd_common_factor_constexpr()
{
  constexpr auto mod = 7_Z;
  using GF = decltype(Zq(mod));

  GF one(1_Z), two(2_Z), three(3_Z), six(6_Z);

  // a = (x - 1)(x - 2) = x^2 + 4x + 2 in GF(7)
  // b = (x - 1)(x - 3) = x^2 + 3x + 3 in GF(7)
  // gcd = (x - 1) = x + 6 in GF(7)

  polynomial_nttp<GF, 2> a{{two, GF(4_Z), one}};
  polynomial_nttp<GF, 2> b{{three, three, one}};

  auto g = poly_gcd(a, b);

  // Should be monic linear: x + 6
  return g[1] == one && g[0] == six;
}

consteval bool test_berlekamp_factor_constexpr()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  GF zero(0_Z), one(1_Z), four(4_Z);

  // Factor x^2 - 1 = (x-1)(x+1) over GF(5)
  polynomial_nttp<GF, 2> f{{four, zero, one}};

  auto [factors, count] = berlekamp_factor<GF, 5, 2>(f, zero, one);

  return count == 2;
}

// ============================================================
// Main
// ============================================================

int main()
{
  // Compile-time verification
  static_assert(test_poly_gcd_constexpr(), "poly_gcd should work at compile time");
  static_assert(test_poly_gcd_coprime_constexpr(), "poly_gcd coprime should work at compile time");
  static_assert(test_poly_gcd_common_factor_constexpr(), "poly_gcd common factor should work at compile time");
  static_assert(test_berlekamp_factor_constexpr(), "berlekamp_factor should work at compile time");

  std::println("Running Berlekamp algorithm tests with ctbignum...\n");
  std::println("Compile-time static_asserts: PASSED ✓");

  if (!test_poly_gcd())
  {
    std::println("test_poly_gcd: FAILED");
    return 1;
  }

  if (!test_poly_gcd_coprime())
  {
    std::println("test_poly_gcd_coprime: FAILED");
    return 1;
  }

  if (!test_poly_gcd_common_factor())
  {
    std::println("test_poly_gcd_common_factor: FAILED");
    return 1;
  }

  if (!test_berlekamp_factor())
  {
    std::println("test_berlekamp_factor: FAILED");
    return 1;
  }

  std::println("\nAll Berlekamp tests passed.");
  return 0;
}
