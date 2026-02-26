/*
 *  poly_rem.cpp – Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit tests for poly_rem over double
 *
 *  Tests are structured as consteval functions verified by static_assert
 *  (compile-time) and also exercised with volatile-sourced dividend
 *  coefficients to force runtime execution (runtime coverage).
 *
 *  Tests resistant to compile-time evaluation: none for these degree ranges.
 *  Tests resistant to runtime: none — Mod is always NTTP (compile-time),
 *    but dividend coefficients can be sourced at runtime.
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)

constexpr bool approx_eq(double a, double b)
{
  constexpr auto absdiff = [](double x, double y) { return x > y ? x - y : y - x; };
  return absdiff(a, b) < 1e-9;
}

// ============================================================
// Degenerate modulus: N == 0
// poly_rem<Mod>(a) returns polynomial_nttp<R, 0>{} for any a
// ============================================================

consteval bool test_poly_rem_zero_degree_mod_constant_dividend()
{
  constexpr polynomial_nttp<double, 0> mod{1.};
  constexpr polynomial_nttp<double, 0> a{5.};
  constexpr auto r = poly_rem<double, 0, mod>(a);
  // N == 0 branch: always returns the zero degree-0 polynomial
  return approx_eq(r[0], 0.);
}

consteval bool test_poly_rem_zero_degree_mod_cubic_dividend()
{
  constexpr polynomial_nttp<double, 0> mod{3.};
  constexpr polynomial_nttp<double, 3> a{{1., 2., 3., 4.}};
  constexpr auto r = poly_rem<double, 0, mod>(a);
  return approx_eq(r[0], 0.);
}

// ============================================================
// M < N: dividend already lies in the reduced space
// Result is the dividend zero-padded to N - 1 coefficients
// ============================================================

consteval bool test_poly_rem_constant_dividend_linear_mod()
{
  // deg(a) = 0 < N = 1: remainder is the constant itself
  constexpr polynomial_nttp<double, 1> mod{{-1., 1.}}; // x - 1
  constexpr polynomial_nttp<double, 0> a{7.};
  constexpr auto r = poly_rem<double, 1, mod>(a);
  // return type is polynomial_nttp<double, 0>
  return approx_eq(r[0], 7.);
}

consteval bool test_poly_rem_constant_dividend_quadratic_mod()
{
  // deg(a) = 0 < N = 2: result is degree-1 type, zero-padded
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}}; // x^2 + 1
  constexpr polynomial_nttp<double, 0> a{3.};
  constexpr auto r = poly_rem<double, 2, mod>(a);
  // return type is polynomial_nttp<double, 1>
  return approx_eq(r[0], 3.) && approx_eq(r[1], 0.);
}

consteval bool test_poly_rem_linear_dividend_quadratic_mod()
{
  // deg(a) = 1 < N = 2: remainder is a itself in the degree-1 result type
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}}; // x^2 + 1
  constexpr polynomial_nttp<double, 1> a{{3., -2.}};      // -2x + 3
  constexpr auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 3.) && approx_eq(r[1], -2.);
}

consteval bool test_poly_rem_linear_dividend_cubic_mod()
{
  // deg(a) = 1 < N = 3: zero-padded into degree-2 result type
  constexpr polynomial_nttp<double, 3> mod{{1., 0., 0., 1.}}; // x^3 + 1
  constexpr polynomial_nttp<double, 1> a{{5., 4.}};           // 4x + 5
  constexpr auto r = poly_rem<double, 3, mod>(a);
  // return type is polynomial_nttp<double, 2>
  return approx_eq(r[0], 5.) && approx_eq(r[1], 4.) && approx_eq(r[2], 0.);
}

// ============================================================
// M == N: single elimination step
// ============================================================

consteval bool test_poly_rem_self_mod_self()
{
  // a mod a = 0
  constexpr polynomial_nttp<double, 1> mod{{-1., 1.}}; // x - 1
  constexpr polynomial_nttp<double, 1> a{{-1., 1.}};
  constexpr auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], 0.);
}

consteval bool test_poly_rem_degree_n_dividend_no_remainder()
{
  // (x - 1)(x + 1) = x^2 - 1; mod (x + 1) = 0
  constexpr polynomial_nttp<double, 1> mod{{1., 1.}}; // x + 1
  constexpr polynomial_nttp<double, 1> a{{-2., 1.}};  // x - 2  (not a multiple of x+1)
  // f(-1) = -1 - 2 = -3; remainder should be -3
  constexpr auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], -3.);
}

// ============================================================
// M > N: full synthetic division
// ============================================================

consteval bool test_poly_rem_cubic_linear_exact()
{
  // x^3 - 7x + 6 has root x = 2: (x - 2) divides it exactly
  // coefficients in ascending degree: 6, -7, 0, 1
  constexpr polynomial_nttp<double, 1> mod{{-2., 1.}};       // x - 2
  constexpr polynomial_nttp<double, 3> a{{6., -7., 0., 1.}}; // x^3 - 7x + 6
  constexpr auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], 0.);
}

consteval bool test_poly_rem_cubic_linear_nonzero()
{
  // x^3 - 7x + 6 evaluated at x = -1 is -1 + 7 + 6 = 12
  // so remainder of division by (x + 1) is 12
  constexpr polynomial_nttp<double, 1> mod{{1., 1.}};        // x + 1
  constexpr polynomial_nttp<double, 3> a{{6., -7., 0., 1.}}; // x^3 - 7x + 6
  constexpr auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], 12.);
}

consteval bool test_poly_rem_quartic_quadratic_exact()
{
  // (x^2 - 1)(x^2 + 1) = x^4 - 1; mod (x^2 + 1): (x^4 - 1) = (x^2 - 1)(x^2 + 1) + 0
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};        // x^2 + 1
  constexpr polynomial_nttp<double, 4> a{{-1., 0., 0., 0., 1.}}; // x^4 - 1
  constexpr auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 0.) && approx_eq(r[1], 0.);
}

consteval bool test_poly_rem_quartic_quadratic_nonzero()
{
  // x^4 + x^2 + 1 mod (x^2 + 1):
  // x^4 + x^2 + 1 = (x^2 + 1)(x^2) + 1, so remainder = 1
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};       // x^2 + 1
  constexpr polynomial_nttp<double, 4> a{{1., 0., 1., 0., 1.}}; // x^4 + x^2 + 1
  constexpr auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 1.) && approx_eq(r[1], 0.);
}

consteval bool test_poly_rem_quartic_quadratic_linear_remainder()
{
  // x^4 + x mod (x^2 + 1):
  // x^4 = (x^2+1)(x^2) - x^2; so x^4 + x = (x^2)(x^2+1) - x^2 + x
  // then -x^2 + x mod (x^2+1): -x^2 = -(x^2+1) + 1, so -x^2 + x = -(x^2+1) + 1 + x
  // remainder = x + 1
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};       // x^2 + 1
  constexpr polynomial_nttp<double, 4> a{{0., 1., 0., 0., 1.}}; // x^4 + x
  constexpr auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 1.) && approx_eq(r[1], 1.);
}

consteval bool test_poly_rem_zero_dividend()
{
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};
  constexpr polynomial_nttp<double, 3> a{}; // zero polynomial
  constexpr auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 0.) && approx_eq(r[1], 0.);
}

// ============================================================
// Reconstruction: a == quotient * Mod + poly_rem<Mod>(a)
// verified by comparing poly_rem with division_prototype's remainder
// ============================================================

consteval bool test_poly_rem_matches_division_prototype_linear()
{
  constexpr polynomial_nttp<double, 1> mod{{-2., 1.}}; // x - 2
  constexpr polynomial_nttp<double, 3> a{{6., -7., 0., 1.}};
  constexpr auto qr = division_prototype<double, 3, a, 1, mod>();
  constexpr auto r_div = qr.second;
  constexpr auto r_rem = poly_rem<double, 1, mod>(a);
  return approx_eq(r_div[0], r_rem[0]);
}

consteval bool test_poly_rem_matches_division_prototype_quadratic()
{
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}}; // x^2 + 1
  constexpr polynomial_nttp<double, 4> a{{1., 0., 1., 0., 1.}};
  constexpr auto qr = division_prototype<double, 4, a, 2, mod>();
  constexpr auto r_div = qr.second;
  constexpr auto r_rem = poly_rem<double, 2, mod>(a);
  return approx_eq(r_div[0], r_rem[0]) && approx_eq(r_div[1], r_rem[1]);
}

consteval bool test_poly_rem_reconstruction_cubic_linear()
{
  // verify a == quotient * mod + remainder
  constexpr polynomial_nttp<double, 1> mod{{-2., 1.}};
  constexpr polynomial_nttp<double, 3> a{{6., -7., 0., 1.}};
  constexpr auto qr = division_prototype<double, 3, a, 1, mod>();
  constexpr auto q = qr.first;
  constexpr auto r_from_div = qr.second;
  constexpr auto r_from_rem = poly_rem<double, 1, mod>(a);
  // consistency between the two remainder computations
  if (!approx_eq(r_from_div[0], r_from_rem[0]))
    return false;
  // reconstruction: q * mod + r should equal a coefficient-wise
  constexpr auto reconstructed = q * mod + r_from_div;
  return approx_eq(reconstructed[0], a[0]) && approx_eq(reconstructed[1], a[1]) && approx_eq(reconstructed[2], a[2]) &&
         approx_eq(reconstructed[3], a[3]);
}

// ============================================================
// Runtime coverage: volatile-sourced dividend coefficients force
// the same code paths to execute at runtime
// ============================================================

bool runtime_test_poly_rem_m_less_than_n()
{
  volatile double v0 = 3., v1 = -2.;
  double c0 = v0, c1 = v1;
  polynomial_nttp<double, 1> a{{c0, c1}};
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};
  auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 3.) && approx_eq(r[1], -2.);
}

bool runtime_test_poly_rem_cubic_linear_exact()
{
  volatile double c0 = 6., c1 = -7., c2 = 0., c3 = 1.;
  polynomial_nttp<double, 3> a{{c0, c1, c2, c3}};
  constexpr polynomial_nttp<double, 1> mod{{-2., 1.}};
  auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], 0.);
}

bool runtime_test_poly_rem_cubic_linear_nonzero()
{
  volatile double c0 = 6., c1 = -7., c2 = 0., c3 = 1.;
  polynomial_nttp<double, 3> a{{c0, c1, c2, c3}};
  constexpr polynomial_nttp<double, 1> mod{{1., 1.}};
  auto r = poly_rem<double, 1, mod>(a);
  return approx_eq(r[0], 12.);
}

bool runtime_test_poly_rem_quartic_quadratic_nonzero()
{
  volatile double c0 = 1., c1 = 0., c2 = 1., c3 = 0., c4 = 1.;
  polynomial_nttp<double, 4> a{{c0, c1, c2, c3, c4}};
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};
  auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 1.) && approx_eq(r[1], 0.);
}

bool runtime_test_poly_rem_zero_dividend()
{
  volatile double z = 0.;
  polynomial_nttp<double, 3> a{{z, z, z, z}};
  constexpr polynomial_nttp<double, 2> mod{{1., 0., 1.}};
  auto r = poly_rem<double, 2, mod>(a);
  return approx_eq(r[0], 0.) && approx_eq(r[1], 0.);
}

// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

static_assert(test_poly_rem_zero_degree_mod_constant_dividend(), "poly_rem: N=0 constant dividend");
static_assert(test_poly_rem_zero_degree_mod_cubic_dividend(), "poly_rem: N=0 cubic dividend");
static_assert(test_poly_rem_constant_dividend_linear_mod(), "poly_rem: M=0 < N=1");
static_assert(test_poly_rem_constant_dividend_quadratic_mod(), "poly_rem: M=0 < N=2, zero-pad");
static_assert(test_poly_rem_linear_dividend_quadratic_mod(), "poly_rem: M=1 < N=2");
static_assert(test_poly_rem_linear_dividend_cubic_mod(), "poly_rem: M=1 < N=3, zero-pad");
static_assert(test_poly_rem_self_mod_self(), "poly_rem: a mod a = 0");
static_assert(test_poly_rem_degree_n_dividend_no_remainder(), "poly_rem: M=N, nonzero remainder");
static_assert(test_poly_rem_cubic_linear_exact(), "poly_rem: cubic mod linear, exact");
static_assert(test_poly_rem_cubic_linear_nonzero(), "poly_rem: cubic mod linear, f(-1)=12");
static_assert(test_poly_rem_quartic_quadratic_exact(), "poly_rem: quartic mod quadratic, exact");
static_assert(test_poly_rem_quartic_quadratic_nonzero(), "poly_rem: quartic mod quadratic, remainder=1");
static_assert(test_poly_rem_quartic_quadratic_linear_remainder(), "poly_rem: quartic mod quadratic, linear remainder");
static_assert(test_poly_rem_zero_dividend(), "poly_rem: zero dividend");
static_assert(test_poly_rem_matches_division_prototype_linear(),
              "poly_rem: consistent with division_prototype (linear)");
static_assert(test_poly_rem_matches_division_prototype_quadratic(),
              "poly_rem: consistent with division_prototype (quadratic)");
static_assert(test_poly_rem_reconstruction_cubic_linear(), "poly_rem: reconstruction a = q*mod + r");

int main()
{
  constexpr bool ct_zero_mod_constant = test_poly_rem_zero_degree_mod_constant_dividend();
  constexpr bool ct_zero_mod_cubic = test_poly_rem_zero_degree_mod_cubic_dividend();
  constexpr bool ct_const_lin = test_poly_rem_constant_dividend_linear_mod();
  constexpr bool ct_const_quad = test_poly_rem_constant_dividend_quadratic_mod();
  constexpr bool ct_linear_quad = test_poly_rem_linear_dividend_quadratic_mod();
  constexpr bool ct_linear_cubic = test_poly_rem_linear_dividend_cubic_mod();
  constexpr bool ct_self = test_poly_rem_self_mod_self();
  constexpr bool ct_degree_n = test_poly_rem_degree_n_dividend_no_remainder();
  constexpr bool ct_cubic_exact = test_poly_rem_cubic_linear_exact();
  constexpr bool ct_cubic_nonzero = test_poly_rem_cubic_linear_nonzero();
  constexpr bool ct_quartic_exact = test_poly_rem_quartic_quadratic_exact();
  constexpr bool ct_quartic_nonzero = test_poly_rem_quartic_quadratic_nonzero();
  constexpr bool ct_quartic_linear_rem = test_poly_rem_quartic_quadratic_linear_remainder();
  constexpr bool ct_zero_dividend = test_poly_rem_zero_dividend();
  constexpr bool ct_matches_div_lin = test_poly_rem_matches_division_prototype_linear();
  constexpr bool ct_matches_div_quad = test_poly_rem_matches_division_prototype_quadratic();
  constexpr bool ct_reconstruction = test_poly_rem_reconstruction_cubic_linear();

  bool rt_m_lt_n = runtime_test_poly_rem_m_less_than_n();
  bool rt_cubic_ex = runtime_test_poly_rem_cubic_linear_exact();
  bool rt_cubic_nz = runtime_test_poly_rem_cubic_linear_nonzero();
  bool rt_quartic = runtime_test_poly_rem_quartic_quadratic_nonzero();
  bool rt_zero = runtime_test_poly_rem_zero_dividend();

  if constexpr (ct_zero_mod_constant && ct_zero_mod_cubic && ct_const_lin && ct_const_quad && ct_linear_quad &&
                ct_linear_cubic && ct_self && ct_degree_n && ct_cubic_exact && ct_cubic_nonzero && ct_quartic_exact &&
                ct_quartic_nonzero && ct_quartic_linear_rem && ct_zero_dividend && ct_matches_div_lin &&
                ct_matches_div_quad && ct_reconstruction)
  {
    if (rt_m_lt_n && rt_cubic_ex && rt_cubic_nz && rt_quartic && rt_zero)
      return 0;
    return 1;
  }
  else
  {
    return 1;
  }
}
