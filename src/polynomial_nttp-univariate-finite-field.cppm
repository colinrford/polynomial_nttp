/*
 *  polynomial_nttp-univariate-finite-field.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module lam.polynomial_nttp:univariate.finite_field;

import std;
import lam.concepts;
import :univariate.structure;
import :univariate.algebra;

namespace lam::polynomial::univariate::finite_field
{

template<typename K>
concept field_element_c_weak = lam::concepts::experimental::field_element_c_weak<K>;

using lam::polynomial::univariate::polynomial_nttp;
using lam::polynomial::univariate::algebra::poly_inv;
using lam::polynomial::univariate::algebra::poly_rem;

// Pick the multiplicative identity of K, preferring K::one() over K(1).
template<field_element_c_weak K>
constexpr K k_one()
{
  if constexpr (requires {
                  { K::one() } -> std::same_as<K>;
                })
    return K::one();
  else
    return K(1);
}

/**
 * finite_field_extension<K, N, Irreducible> represents the quotient ring K[x]/(Irreducible),
 * where Irreducible is an irreducible polynomial of degree N over the field K.
 * Elements are stored as polynomials of degree < N (i.e., polynomial_nttp<K, N-1>).
 *
 * Satisfies field_element_c_weak via:
 *   additive identity:       ::zero()
 *   multiplicative identity: ::one()
 *   inverse:                 inv(a)  (free function, found via ADL)
 *   division:                a / b   (defined as a * inv(b))
 */
export 
template<field_element_c_weak K, std::size_t N, polynomial_nttp<K, N> Irreducible>
  requires(N >= 1)
struct finite_field_extension
{
  polynomial_nttp<K, N - 1> data{};

  constexpr finite_field_extension() = default;

  explicit constexpr finite_field_extension(polynomial_nttp<K, N - 1> p) noexcept : data{p} {}

  static constexpr finite_field_extension zero() noexcept { return {}; }

  static constexpr finite_field_extension one() noexcept
  { return finite_field_extension{polynomial_nttp<K, N - 1>{k_one<K>()}}; }

  constexpr finite_field_extension operator+(const finite_field_extension& rhs) const noexcept
  { return finite_field_extension{data + rhs.data}; }

  constexpr finite_field_extension operator-(const finite_field_extension& rhs) const noexcept
  { return finite_field_extension{data - rhs.data}; }

  constexpr finite_field_extension operator-() const noexcept { return finite_field_extension{-data}; }

  // Multiplication: multiply as polynomials, then reduce mod Irreducible.
  constexpr finite_field_extension operator*(const finite_field_extension& rhs) const noexcept
  { return finite_field_extension{poly_rem<K, N, Irreducible>(data * rhs.data)}; }

  // Equality: compare underlying coefficient arrays directly.
  constexpr bool operator==(const finite_field_extension& rhs) const noexcept
  { return data.coefficients == rhs.data.coefficients; }

  // Coefficient access (index 0 = constant term).
  constexpr K operator[](std::size_t i) const noexcept { return data[i]; }
};

// inv(a): multiplicative inverse of a in K[x]/(Irreducible), via extended GCD.
// Found by ADL from the finite_field_extension namespace.
export 
template<field_element_c_weak K, std::size_t N, polynomial_nttp<K, N> Irreducible>
constexpr auto inv(finite_field_extension<K, N, Irreducible> a)
{ return finite_field_extension<K, N, Irreducible>{poly_inv<K, N, Irreducible>(a.data)}; }

// Division: a / b = a * inv(b).
export 
template<field_element_c_weak K, std::size_t N, polynomial_nttp<K, N> Irreducible>
constexpr auto operator/(finite_field_extension<K, N, Irreducible> a, finite_field_extension<K, N, Irreducible> b)
{ return a * inv(b); }

} // end namespace lam::polynomial::univariate::finite_field

namespace lam::polynomial
{
export using univariate::finite_field::finite_field_extension;
} // end namespace lam::polynomial

namespace lam
{
export using polynomial::univariate::finite_field::finite_field_extension;
} // end namespace lam
