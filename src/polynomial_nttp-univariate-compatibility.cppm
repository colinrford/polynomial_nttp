/*
 *  polynomial_nttp-univariate-compatibility.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module lam.polynomial_nttp:univariate.compat;

import std;
import :univariate.structure;
import :univariate.algebra;
import :univariate.roots;

namespace lam::polynomial::univariate::compat
{

// ============================================================
// Type Capability Concepts for Root Finding
// ============================================================

export template<typename T>
concept has_sqrt = requires(T x) 
{
  { std::sqrt(x) } -> std::convertible_to<T>;
};

export template<typename T>
concept has_optional_sqrt = requires(T x) 
{
  { sqrt(x) } -> std::same_as<std::optional<T>>;
};

export template<typename T>
concept has_sqrt_like = has_sqrt<T> || has_optional_sqrt<T>;

export template<typename T>
concept has_cbrt = requires(T x) 
{
  { std::cbrt(x) } -> std::convertible_to<T>;
};

export template<typename T>
concept has_optional_cbrt = requires(T x) 
{
  { cbrt(x) } -> std::same_as<std::optional<T>>;
};

export template<typename T>
concept has_cbrt_like = has_cbrt<T> || has_optional_cbrt<T>;
export template<typename T>
concept has_nth_root = requires(T x, std::size_t n) 
{
  { nth_root(x, n) } -> std::same_as<std::optional<T>>;
};

// Generic compatibility verification for type T
// T must verify:
// 1. Structure (construction, access)
// 2. Arithmetic (+, -, *) with polynomials
// 3. Calculus (derivative, antiderivative)
// 4. Division (division_prototype)
// 5. Roots (roots())

export 
template<ring_element_c_weak T>
constexpr void verify_structure(T one, T two, T three)
{
  constexpr std::size_t N = 2;
  polynomial_nttp<T, N> p{{one, two, three}};
  static_assert(N == 2);

  if (p[0] != one)
    throw "Access check failed at index 0";
  if (p[1] != two)
    throw "Access check failed at index 1";
}

export 
template<ring_element_c_weak T>
constexpr void verify_algebra(T one)
{
  polynomial_nttp<T, 1> p1{{one, one}};
  polynomial_nttp<T, 1> p2{{one, -one}};
  auto sum = p1 + p2;
  auto diff = p1 - p2;
  auto prod = p1 * p2;
}

// Requires field_element_c_weak
export 
template<field_element_c_weak T>
constexpr void verify_field_algebra(T one, T two)
{
  auto q = two / two;
  if (q != one)
    throw "Division failed: two / two != one";

  auto four = two + two;
  if (four / two != two)
    throw "Division failed: four / two != two";

  auto inv_two = one / two;
  if (two * inv_two != one)
    throw "Multiplicative inverse failed";
}

// GF(17) specific check: verify 1/2 = 9
export 
template<field_element_c_weak T>
constexpr void verify_field_inverse_value(T one, T two, T nine)
{
  auto inv_two = one / two;
  if (inv_two != nine)
    throw "GF(17): 1/2 should equal 9";
  if (two * nine != one)
    throw "GF(17): 2 * 9 should equal 1";
}

// Polynomial evaluation test
export 
template<ring_element_c_weak T>
constexpr void verify_poly_evaluation(T zero, T one, T two, T three)
{
  polynomial_nttp<T, 2> p{{one, two, three}};

  auto val_at_0 = p(zero);
  if (val_at_0 != one)
    throw "Evaluation p(0) failed";

  auto val_at_1 = p(one);
  auto six = one + two + three;
  if (val_at_1 != six)
    throw "Evaluation p(1) failed";

  auto val_at_2 = p(two);
}

// Polynomial multiplication coefficient verification
export 
template<ring_element_c_weak T>
constexpr void verify_poly_mult_coeffs(T zero, T one)
{
  polynomial_nttp<T, 1> p{{one, one}};
  auto p_squared = p * p;

  auto two = one + one;
  if (p_squared[0] != one)
    throw "Mult coeff [0] failed";
  if (p_squared[1] != two)
    throw "Mult coeff [1] failed";
  if (p_squared[2] != one)
    throw "Mult coeff [2] failed";
}

// Algebraic identity: (a + b)(a - b) = a² - b²
export 
template<ring_element_c_weak T>
constexpr void verify_algebraic_identity(T zero, T one)
{
  polynomial_nttp<T, 1> a_plus_b{{one, one}};
  auto neg_one = zero - one;
  polynomial_nttp<T, 1> a_minus_b{{neg_one, one}};

  auto product = a_plus_b * a_minus_b;

  if (product[0] != neg_one)
    throw "Identity coeff [0] failed";
  if (product[1] != zero)
    throw "Identity coeff [1] failed";
  if (product[2] != one)
    throw "Identity coeff [2] failed";
}

// Division with non-trivial remainder
export 
template<ring_element_c_weak R, std::size_t M, polynomial_nttp<R, M> dividend, std::size_t N, polynomial_nttp<R, N> divisor>
constexpr auto verify_division_with_remainder()
{
  auto [q, r] = division_prototype<R, M, dividend, N, divisor>();
  return std::make_pair(q, r);
}

// Higher degree polynomial test
export 
template<ring_element_c_weak T>
constexpr void verify_higher_degree(T zero, T one, T two, T three)
{
  polynomial_nttp<T, 4> p{{one, two, three, two, one}};

  auto dp = derivative(p);
  auto four = two + two;
  auto six = three + three;

  if (dp[0] != two)
    throw "Higher degree derivative [0] failed";
  if (dp[1] != six)
    throw "Higher degree derivative [1] failed";
  if (dp[2] != six)
    throw "Higher degree derivative [2] failed";
  if (dp[3] != four)
    throw "Higher degree derivative [3] failed";
}

export 
template<ring_element_c_weak T, ring_element_c_weak Zero>
constexpr void verify_calculus(Zero zero, T one)
{
  polynomial_nttp<T, 2> p{{zero, zero, one}};
  auto d = derivative(p);
}

// Division_prototype requires NTTPs
export 
template<ring_element_c_weak R, std::size_t M, polynomial_nttp<R, M> dividend, std::size_t N, polynomial_nttp<R, N> divisor>
constexpr auto verify_division_nttp()
{
  auto [q, r] = division_prototype<R, M, dividend, N, divisor>();
  return q;
}

// Linear root test: ax + b = 0 → x = -b/a
// This tests the roots_degree_1() function for degree-1 polynomials
export 
template<ring_element_c_weak T>
constexpr bool verify_linear_roots(T a, T b, T expected_root)
{
  // Create linear polynomial: b + a*x
  polynomial_nttp<T, 1> p{{b, a}};

  // Find roots using degree-1 solver
  auto r = roots_degree_1(p);

  // Should find exactly 1 root
  if (r.size() != 1)
    return false;

  // Check root value matches expected
  if (r[0].value != expected_root)
    return false;

  // Verify: p(root) should equal zero
  auto zero = b - b; // Get zero element
  if (p(expected_root) != zero)
    return false;

  return true;
}

} // namespace lam::polynomial::univariate::compat
