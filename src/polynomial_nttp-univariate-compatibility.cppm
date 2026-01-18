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
// Type Capability Concepts
// ============================================================

export template<typename T>
concept is_exact_arithmetic = std::integral<T> ||                         // Built-in integers
                              std::numeric_limits<T>::is_exact ||         // Explicitly specialized types
                              (!std::numeric_limits<T>::is_specialized && // If no specialization exists...
                               !std::floating_point<T>);                  // ...and not a built-in float, assume exact

export template<typename T>
concept is_approximate_arithmetic = std::floating_point<T> || !is_exact_arithmetic<T>;

export template<typename T>
concept has_sqrt = requires(T x) {
  { std::sqrt(x) } -> std::convertible_to<T>;
};

export template<typename T>
concept has_optional_sqrt = requires(T x) {
  { sqrt(x) } -> std::same_as<std::optional<T>>;
};

export template<typename T>
concept has_sqrt_like = has_sqrt<T> || has_optional_sqrt<T>;

export template<typename T>
concept has_cbrt = requires(T x) {
  { std::cbrt(x) } -> std::convertible_to<T>;
};

export template<typename T>
concept has_optional_cbrt = requires(T x) {
  { cbrt(x) } -> std::same_as<std::optional<T>>;
};

export template<typename T>
concept has_cbrt_like = has_cbrt<T> || has_optional_cbrt<T>;

export template<typename T>
concept has_nth_root = requires(T x, std::size_t n) {
  { nth_root(x, n) } -> std::same_as<std::optional<T>>;
};

// ============================================================
// Capability Report Struct
// ============================================================

export template<typename T>
struct type_capabilities
{
  static constexpr bool is_ring = ring_element_c_weak<T>;
  static constexpr bool is_field = field_element_c_weak<T>;
  static constexpr bool is_exact = is_exact_arithmetic<T>;
  static constexpr bool is_approximate = is_approximate_arithmetic<T>;
  static constexpr bool can_sqrt = has_sqrt_like<T>;
  static constexpr bool can_cbrt = has_cbrt_like<T>;
  static constexpr bool can_nth_root = has_nth_root<T>;

  // Derived capabilities
  static constexpr bool can_find_linear_roots = is_field;
  static constexpr bool can_find_quadratic_roots = is_field && can_sqrt;
  static constexpr bool can_find_cubic_roots = is_field && can_cbrt;
};

// Generic compatibility verification for type T
// T must verify:
// 1. Structure (construction, access)
// 2. Arithmetic (+, -, *) with polynomials
// 3. Calculus (derivative, antiderivative)
// 4. Division (division_prototype)
// 5. Roots (roots())

// ============================================================
// Exact Arithmetic Verification
// ============================================================

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_structure(T one, T two, T three)
{
  constexpr std::size_t N = 2;
  polynomial_nttp<T, N> p{{one, two, three}};
  static_assert(N == 2);

  if (p[0] != one)
    throw "Access check failed at index 0";
  if (p[1] != two)
    throw "Access check failed at index 1";
}

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_algebra(T one)
{
  polynomial_nttp<T, 1> p1{{one, one}};
  polynomial_nttp<T, 1> p2{{one, -one}};
  auto sum = p1 + p2;
  auto diff = p1 - p2;
  auto prod = p1 * p2;
}

export template<field_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_field_algebra(T one, T two)
{
  auto q = two / two;
  if (q != one)
    throw "Division failed";

  auto four = two + two;
  if (four / two != two)
    throw "Division of sum failed";

  auto inv_two = one / two;
  if (two * inv_two != one)
    throw "Multiplicative inverse failed";
}

// GF(17) specific check
export template<field_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_field_inverse_value(T one, T two, T nine)
{
  auto inv_two = one / two;
  if (inv_two != nine)
    throw "GF(17): 1/2 != 9";
  if (two * nine != one)
    throw "GF(17): 2 * 9 != 1";
}

// ============================================================
// Approximate Arithmetic Verification
// ============================================================

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_structure(T one, T two, T three, std::optional<double> tol = std::nullopt)
{
  constexpr std::size_t N = 2;
  polynomial_nttp<T, N> p{{one, two, three}};

  if (!is_approx_equal(p[0], one, tol))
    throw "Access check failed at index 0";
  if (!is_approx_equal(p[1], two, tol))
    throw "Access check failed at index 1";
}

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_algebra(T one, std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 1> p1{{one, one}};
  polynomial_nttp<T, 1> p2{{one, -one}};
  auto sum = p1 + p2;
  auto diff = p1 - p2;
  auto prod = p1 * p2;
}

export template<field_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_field_algebra(T one, T two, std::optional<double> tol = std::nullopt)
{
  auto q = two / two;
  if (!is_approx_equal(q, one, tol))
    throw "Division failed";

  auto four = two + two;
  if (!is_approx_equal(four / two, two, tol))
    throw "Division of sum failed";

  auto inv_two = one / two;
  if (!is_approx_equal(two * inv_two, one, tol))
    throw "Multiplicative inverse failed";
}

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_poly_evaluation(T zero, T one, T two, T three)
{
  polynomial_nttp<T, 2> p{{one, two, three}};

  auto val_at_0 = p(zero);
  if (val_at_0 != one)
    throw "Evaluation p(0) failed";

  auto val_at_1 = p(one);
  auto six = one + two + three;
  if (val_at_1 != six)
    throw "Evaluation p(1) failed";
}

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_poly_evaluation(T zero, T one, T two, T three,
                                                  std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 2> p{{one, two, three}};

  auto val_at_0 = p(zero);
  if (!is_approx_equal(val_at_0, one, tol))
    throw "Evaluation p(0) failed";

  auto val_at_1 = p(one);
  auto six = one + two + three;
  if (!is_approx_equal(val_at_1, six, tol))
    throw "Evaluation p(1) failed";
}

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_poly_mult_coeffs(T zero, T one)
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

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_poly_mult_coeffs(T zero, T one, std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 1> p{{one, one}};
  auto p_squared = p * p;

  auto two = one + one;
  if (!is_approx_equal(p_squared[0], one, tol))
    throw "Mult coeff [0] failed";
  if (!is_approx_equal(p_squared[1], two, tol))
    throw "Mult coeff [1] failed";
  if (!is_approx_equal(p_squared[2], one, tol))
    throw "Mult coeff [2] failed";
}

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_algebraic_identity(T zero, T one)
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

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_algebraic_identity(T zero, T one, std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 1> a_plus_b{{one, one}};
  auto neg_one = zero - one;
  polynomial_nttp<T, 1> a_minus_b{{neg_one, one}};

  auto product = a_plus_b * a_minus_b;

  if (!is_approx_equal(product[0], neg_one, tol))
    throw "Identity coeff [0] failed";
  if (!is_approx_equal(product[1], zero, tol))
    throw "Identity coeff [1] failed";
  if (!is_approx_equal(product[2], one, tol))
    throw "Identity coeff [2] failed";
}

// Division utilities (kept generic as they don't perform assertions)
export template<ring_element_c_weak R, std::size_t M, polynomial_nttp<R, M> dividend, std::size_t N,
                polynomial_nttp<R, N> divisor>
constexpr auto verify_division_with_remainder()
{
  auto [q, r] = division_prototype<R, M, dividend, N, divisor>();
  return std::make_pair(q, r);
}

export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr void verify_exact_higher_degree(T zero, T one, T two, T three)
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

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr void verify_approximate_higher_degree(T zero, T one, T two, T three, std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 4> p{{one, two, three, two, one}};

  auto dp = derivative(p);
  auto four = two + two;
  auto six = three + three;

  if (!is_approx_equal(dp[0], two, tol))
    throw "Higher degree derivative [0] failed";
  if (!is_approx_equal(dp[1], six, tol))
    throw "Higher degree derivative [1] failed";
  if (!is_approx_equal(dp[2], six, tol))
    throw "Higher degree derivative [2] failed";
  if (!is_approx_equal(dp[3], four, tol))
    throw "Higher degree derivative [3] failed";
}

export template<ring_element_c_weak T, ring_element_c_weak Zero>
constexpr void verify_calculus(Zero zero, T one)
{
  polynomial_nttp<T, 2> p{{zero, zero, one}};
  auto d = derivative(p);
  (void)d; // Suppress unused var
}

export template<ring_element_c_weak R, std::size_t M, polynomial_nttp<R, M> dividend, std::size_t N,
                polynomial_nttp<R, N> divisor>
constexpr auto verify_division_nttp()
{
  auto [q, r] = division_prototype<R, M, dividend, N, divisor>();
  return q;
}

// Linear root tests
export template<ring_element_c_weak T>
  requires is_exact_arithmetic<T>
constexpr bool verify_exact_linear_roots(T a, T b, T expected_root)
{
  polynomial_nttp<T, 1> p{{b, a}}; // ax + b

  auto r = roots_degree_1(p);
  if (r.size() != 1)
    return false;
  if (r[0].value != expected_root)
    return false;

  auto zero = b - b;
  if (p(expected_root) != zero)
    return false;

  return true;
}

export template<ring_element_c_weak T>
  requires is_approximate_arithmetic<T>
constexpr bool verify_approximate_linear_roots(T a, T b, T expected_root, std::optional<double> tol = std::nullopt)
{
  polynomial_nttp<T, 1> p{{b, a}}; // ax + b

  auto r = roots_degree_1(p);
  if (r.size() != 1)
    return false;
  if (!is_approx_equal(r[0].value, expected_root, tol))
    return false;

  auto zero = b - b;
  if (!is_approx_equal(p(expected_root), zero, tol))
    return false;

  return true;
}

} // namespace lam::polynomial::univariate::compat
