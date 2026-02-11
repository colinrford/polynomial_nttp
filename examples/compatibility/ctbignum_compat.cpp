/*
 *  ctbignum_compat.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Compile-time verification of polynomial_nttp compatibility with ctbignum.
 */
import std;
import lam.polynomial_nttp;
import lam.concepts;
import lam.ctbignum;

namespace compat = lam::polynomial::univariate::compat;

int main()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  try
  {
    constexpr auto modulus = 17_Z;
    using GF = decltype(Zq(modulus));

    constexpr GF zero = GF(0_Z);
    constexpr GF one = GF(1_Z);
    constexpr GF two = GF(2_Z);
    constexpr GF three = GF(3_Z);

    std::println("Verifying structural and algebraic compatibility...");

    static_assert(lam::concepts::experimental::ring_element_c_weak<GF>, "Zq must satisfy ring_element_c_weak");
    static_assert(lam::concepts::experimental::field_element_c_weak<GF>, "Zq must satisfy field_element_c_weak");

    compat::verify_exact_structure(one, two, three);
    compat::verify_exact_algebra(one);
    compat::verify_exact_field_algebra(one, two);

    // finite field tests
    compat::verify_exact_poly_evaluation(zero, one, two, three);
    compat::verify_exact_poly_mult_coeffs(zero, one);
    compat::verify_exact_algebraic_identity(zero, one);
    compat::verify_exact_higher_degree(zero, one, two, three);

    std::println("Runtime check passed: ctbignum Zq elements work with polynomial_nttp!");
  }
  catch (const char* msg)
  {
    std::print(std::cerr, "Test FAILED: {}\n", msg);
    return 1;
  }
  catch (...)
  {
    std::print(std::cerr, "Test FAILED: Unknown exception\n");
    return 1;
  }
  return 0;
}

// Compile-time check
constexpr bool test_ctbignum_constexpr()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;
  constexpr auto mod = 17_Z;
  using GF = decltype(Zq(mod));

  static_assert(lam::concepts::experimental::ring_element_c_weak<GF>);
  static_assert(lam::concepts::experimental::field_element_c_weak<GF>);

  compat::verify_exact_structure(GF(1_Z), GF(2_Z), GF(3_Z));
  compat::verify_exact_algebra(GF(1_Z));
  compat::verify_exact_field_algebra(GF(1_Z), GF(2_Z));

  // GF(17) specific: verify 1/2 = 9
  compat::verify_exact_field_inverse_value(GF(1_Z), GF(2_Z), GF(9_Z));

  // finite field tests at compile time
  compat::verify_exact_poly_evaluation(GF(0_Z), GF(1_Z), GF(2_Z), GF(3_Z));
  compat::verify_exact_poly_mult_coeffs(GF(0_Z), GF(1_Z));
  compat::verify_exact_algebraic_identity(GF(0_Z), GF(1_Z));
  compat::verify_exact_higher_degree(GF(0_Z), GF(1_Z), GF(2_Z), GF(3_Z));

  // Test division_prototype with NTTP polynomials
  // (x^2 - 1) / (x - 1) = (x + 1)
  constexpr GF neg_one = GF(0_Z) - GF(1_Z); // -1 mod 17 = 16
  constexpr GF zero = GF(0_Z);
  constexpr GF one = GF(1_Z);
  constexpr GF two = GF(2_Z);

  constexpr lam::polynomial_nttp<GF, 2> dividend{{neg_one, zero, one}}; // x^2 - 1
  constexpr lam::polynomial_nttp<GF, 1> divisor{{neg_one, one}};        // x - 1

  auto quotient = compat::verify_division_nttp<GF, 2, dividend, 1, divisor>();
  // quotient should be (x + 1) i.e. {1, 1}
  if (quotient[0] != one)
    throw "Division quotient [0] failed";
  if (quotient[1] != one)
    throw "Division quotient [1] failed";

  // Division with remainder: x^2 / (x - 1) = x + 1 remainder 1
  constexpr lam::polynomial_nttp<GF, 2> x_squared{{zero, zero, one}};
  auto [q2, r2] = compat::verify_division_with_remainder<GF, 2, x_squared, 1, divisor>();
  // remainder should be 1
  if (r2[0] != one)
    throw "Division remainder failed";

  // ============================================================
  // NON-TRIVIAL DIVISION TESTS
  // ============================================================

  constexpr GF three = GF(3_Z);
  constexpr GF four = GF(4_Z);
  constexpr GF five = GF(5_Z);

  // Test 1: Higher degree division
  // (x^4 + x^3 + x^2 + x + 1) / (x^2 + 1)
  // In GF(17): x^4 + x^3 + x^2 + x + 1 = (x^2 + 1)(x^2 + x) + (1)
  // So quotient = x^2 + x, remainder = 1
  constexpr lam::polynomial_nttp<GF, 4> high_deg{{one, one, one, one, one}};
  constexpr lam::polynomial_nttp<GF, 2> quadratic_divisor{{one, zero, one}}; // x^2 + 1

  auto [q3, r3] = compat::verify_division_with_remainder<GF, 4, high_deg, 2, quadratic_divisor>();
  // Quotient should be x^2 + x = {0, 1, 1}
  if (q3[0] != zero)
    throw "High degree division quotient [0] failed";
  if (q3[1] != one)
    throw "High degree division quotient [1] failed";
  if (q3[2] != one)
    throw "High degree division quotient [2] failed";

  // Test 2: Division requiring modular inverse
  // (2x^2 + 4x + 2) / (2x + 2) in GF(17)
  // = (2(x^2 + 2x + 1)) / (2(x + 1)) = (x^2 + 2x + 1) / (x + 1) = x + 1
  // But with leading coeff 2, we need modular inverse of 2 (which is 9 in GF(17))
  constexpr lam::polynomial_nttp<GF, 2> scaled_poly{{two, four, two}}; // 2x^2 + 4x + 2
  constexpr lam::polynomial_nttp<GF, 1> scaled_divisor{{two, two}};    // 2x + 2

  auto [q4, r4] = compat::verify_division_with_remainder<GF, 2, scaled_poly, 1, scaled_divisor>();
  // Quotient should be x + 1 = {1, 1}
  if (q4[0] != one)
    throw "Modular inverse division quotient [0] failed";
  if (q4[1] != one)
    throw "Modular inverse division quotient [1] failed";
  // Remainder should be 0
  if (r4[0] != zero)
    throw "Modular inverse division remainder failed";

  // Test 3: Polynomial division with remainder
  // x^3 + 1 divided by x^2 + x + 1 in GF(17)
  // Long division:
  //   x^3 + 1 = (x^2 + x + 1) * q + r
  //   Leading term: x^3 / x^2 = x, so first term of q is x
  //   x * (x^2 + x + 1) = x^3 + x^2 + x
  //   (x^3 + 1) - (x^3 + x^2 + x) = -x^2 - x + 1 = 16x^2 + 16x + 1 in GF(17)
  //   16x^2 / x^2 = 16, so second term of q is 16
  //   16 * (x^2 + x + 1) = 16x^2 + 16x + 16
  //   (16x^2 + 16x + 1) - (16x^2 + 16x + 16) = 1 - 16 = -15 = 2 in GF(17)
  // So: q = x + 16 = {16, 1}, r = 2 = {2}
  constexpr lam::polynomial_nttp<GF, 3> cubic{{one, zero, zero, one}};  // x^3 + 1
  constexpr lam::polynomial_nttp<GF, 2> irred_divisor{{one, one, one}}; // x^2 + x + 1

  auto [q5, r5] = compat::verify_division_with_remainder<GF, 3, cubic, 2, irred_divisor>();
  // Quotient = x + 16 = {16, 1}
  if (q5[0] != neg_one)
    throw "Irreducible division quotient [0] failed"; // 16 = -1 mod 17
  if (q5[1] != one)
    throw "Irreducible division quotient [1] failed";
  // Remainder = 2
  if (r5[0] != two)
    throw "Irreducible division remainder failed";

  return true;
}

static_assert(test_ctbignum_constexpr(), "ctbignum should work at compile time");

// ============================================================
// LINEAR ROOT TESTS IN MULTIPLE FINITE FIELDS
// Each test demonstrates non-trivial modular arithmetic
// ============================================================

// GF(2) - Characteristic 2
// x + 1 = 0 → x = -1 = 1 (since -1 ≡ 1 mod 2)
constexpr bool test_linear_roots_gf2()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;
  using GF2 = decltype(Zq(2_Z));

  constexpr GF2 one = GF2(1_Z);

  // In GF(2): 1 + 1 = 0, so -1 = 1
  // Polynomial: x + 1 = 0 → x = 1
  return compat::verify_exact_linear_roots(one, one, one);
}

// GF(7) - Characteristic 7
// 5x + 3 = 0 → x = -3/5
// -3 mod 7 = 4
// 5⁻¹ mod 7 = 3 (since 5 * 3 = 15 = 1 mod 7)
// x = 4 * 3 = 12 mod 7 = 5
// Verification: 5*5 + 3 = 25 + 3 = 28 = 0 mod 7 ✓
constexpr bool test_linear_roots_gf7()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;
  using GF7 = decltype(Zq(7_Z));

  constexpr GF7 three = GF7(3_Z);
  constexpr GF7 five = GF7(5_Z);

  // 5x + 3 = 0 → x = 5 (interesting: 5*5 = 25 mod 7 = 4, and 4 + 3 = 7 = 0)
  return compat::verify_exact_linear_roots(five, three, five);
}

// GF(11) - Characteristic 11
// 4x + 9 = 0 → x = -9/4
// -9 mod 11 = 2
// 4⁻¹ mod 11 = 3 (since 4 * 3 = 12 = 1 mod 11)
// x = 2 * 3 = 6
// Verification: 4*6 + 9 = 24 + 9 = 33 = 0 mod 11 ✓
constexpr bool test_linear_roots_gf11()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;
  using GF11 = decltype(Zq(11_Z));

  constexpr GF11 four = GF11(4_Z);
  constexpr GF11 nine = GF11(9_Z);
  constexpr GF11 six = GF11(6_Z);

  return compat::verify_exact_linear_roots(four, nine, six);
}

// GF(13) - Characteristic 13
// 8x + 5 = 0 → x = -5/8
// -5 mod 13 = 8
// 8⁻¹ mod 13 = 5 (since 8 * 5 = 40 = 1 mod 13)
// x = 8 * 5 = 40 mod 13 = 1
// Verification: 8*1 + 5 = 13 = 0 mod 13 ✓
constexpr bool test_linear_roots_gf13()
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;
  using GF13 = decltype(Zq(13_Z));

  constexpr GF13 one = GF13(1_Z);
  constexpr GF13 five = GF13(5_Z);
  constexpr GF13 eight = GF13(8_Z);

  return compat::verify_exact_linear_roots(eight, five, one);
}

// Compile-time verification of all finite field root tests
static_assert(test_linear_roots_gf2(), "Linear root finding should work in GF(2)");
static_assert(test_linear_roots_gf7(), "Linear root finding should work in GF(7)");
static_assert(test_linear_roots_gf11(), "Linear root finding should work in GF(11)");
static_assert(test_linear_roots_gf13(), "Linear root finding should work in GF(13)");

// ============================================================
// Concept Verification: has_sqrt / has_cbrt
// ============================================================

// GF(17) type alias for concept checking
using namespace lam::cbn::literals;
constexpr auto mod_17 = 17_Z;
using GF17 = decltype(lam::cbn::Zq(mod_17));

// Verify: double HAS sqrt/cbrt (via std::sqrt/std::cbrt)
static_assert(compat::has_sqrt<double>, "double should satisfy has_sqrt");
static_assert(compat::has_cbrt<double>, "double should satisfy has_cbrt");

// Verify: finite fields do NOT have sqrt/cbrt (no Tonelli-Shanks/AMM yet)
static_assert(!compat::has_sqrt<GF17>, "GF(17) should NOT satisfy has_sqrt (no Tonelli-Shanks)");
static_assert(!compat::has_cbrt<GF17>, "GF(17) should NOT satisfy has_cbrt (no AMM)");

// Verify: ctbignum supports optional sqrt (found unexpectedly!)
static_assert(compat::has_optional_sqrt<GF17>, "GF(17) SHOULD satisfy has_optional_sqrt (ctbignum provides it!)");

// Verify: check if ctbignum has optional cbrt
static_assert(compat::has_optional_cbrt<GF17>, "GF(17) should satisfy has_optional_cbrt?");

// Test direct usage of sqrt in GF(17)
constexpr bool test_zq_sqrt()
{
  using namespace lam::cbn::literals;
  auto z13 = GF17(13_Z);
  auto z8 = GF17(8_Z);

  auto res = sqrt(z13);
  if (!res.has_value())
    return false;

  // sqrt(13) is 8 or 9 (since 9*9=81=13, 8*8=64=13)
  return (*res == z8) || (*res == -z8);
}
static_assert(test_zq_sqrt(), "Zq sqrt implementation should work correctly");

// Helper for quadratic roots
// x^2 + 2x + 2 = 0 in GF(17)
// Roots should be 3 and 12
constexpr bool test_quadratic_roots_gf17()
{
  using namespace lam::cbn::literals;

  auto one = GF17(1_Z);
  auto two = GF17(2_Z);

  // Polynomial: 1x^2 + 2x + 2 = 0
  lam::polynomial::univariate::polynomial_nttp<GF17, 2> p{{two, two, one}};

  // Solve using roots() which now should use ctbignum's sqrt!
  // Note: must call roots_degree_2 directly or fully qualified roots to avoid ambiguity
  auto r = lam::polynomial::univariate::roots::roots_degree_2(p);

  if (r.size() != 2)
    return false;

  // Roots are likely unordered, so check existence
  // Expected: 3 and 12
  auto r1 = GF17(3_Z);
  auto r2 = GF17(12_Z);

  bool found_r1 = (r[0].value == r1 || r[1].value == r1);
  bool found_r2 = (r[0].value == r2 || r[1].value == r2);

  if (!found_r1 || !found_r2)
    return false;

  // Verify p(root) == 0
  auto val1 = p(r[0].value);
  auto val2 = p(r[1].value);
  auto zero = GF17(0_Z);

  return (val1 == zero) && (val2 == zero);
}
// ============================================================
// LARGE PRIME NTT TEST (N >= 64)
// Verify that the O(N log N) path works at compile time
// ============================================================

// Must specialize finite_field_traits so NTT engine recognizes ZqElement
namespace lam::polynomial::univariate
{
template<typename T, T... Modulus>
struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
{
  static constexpr bool is_finite_field = true;
  // Extract first limb of modulus (assuming single-limb for this test)
  static constexpr T modulus = []() {
    constexpr T mods[] = {Modulus...};
    return mods[0];
  }();
};
} // namespace lam::polynomial::univariate

constexpr bool test_large_ntt_compile_time()
{
  using namespace lam::cbn::literals;
  // Solinas Prime: 29 * 2^57 + 1
  constexpr auto large_prime = 4179340454199820289_Z;
  using large_field = decltype(lam::cbn::Zq(large_prime));

  // Create polynomials of degree 32 (result degree 64 -> triggers NTT threshold 64)
  // Actually, to trigger N >= 64 in operator* (M+N >= 64), we can use M=32, N=32.
  constexpr std::size_t Deg = 32;
  lam::polynomial::univariate::polynomial_nttp<large_field, Deg> p;
  lam::polynomial::univariate::polynomial_nttp<large_field, Deg> q;

  // Simple inputs: p = x + 1, q = x - 1
  // p * q = x^2 - 1
  using wrapper = lam::polynomial::univariate::finite_field_traits<large_field>;

  // Initialize to zero
  for (std::size_t i = 0; i <= Deg; ++i)
  {
    p.coefficients[i] = large_field(0_Z);
    q.coefficients[i] = large_field(0_Z);
  }

  p.coefficients[0] = large_field(1_Z);
  p.coefficients[1] = large_field(1_Z); // p = 1 + x

  // q = x - 1 = -1 + x
  // -1 mod P = P - 1
  auto neg_one = large_field(0_Z) - large_field(1_Z);
  q.coefficients[0] = neg_one;
  q.coefficients[1] = large_field(1_Z);

  // Multiply
  // This should trigger the NTT path because Deg+Deg = 64 >= 64
  auto result = p * q;

  // Expected: -1 + 0x + x^2 + 0x^3 ...
  if (result[0] != neg_one)
    return false;
  if (result[1] != large_field(0_Z))
    return false;
  if (result[2] != large_field(1_Z))
    return false;

  // Verify higher terms are zero
  for (std::size_t i = 3; i <= 2 * Deg; ++i)
  {
    if (result[i] != large_field(0_Z))
      return false;
  }

  return true;
}

static_assert(test_large_ntt_compile_time(), "Large Prime NTT (N=64) must work at compile time!");
