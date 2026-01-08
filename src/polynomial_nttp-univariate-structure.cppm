/*
 *  polynomial_nttp-univariate-structure.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module lam.polynomial_nttp:univariate.structure;

import std;
import lam.concepts;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace lam::polynomial::univariate
{

template<typename R>
concept ring_element_c_weak = lam::concepts::experimental::ring_element_c_weak<R>;
template<typename R>
concept field_element_c_weak = lam::concepts::experimental::field_element_c_weak<R>;

export
template<typename T>
constexpr auto get_epsilon() 
{
  if constexpr (std::is_floating_point_v<T>) 
    return std::numeric_limits<T>::epsilon();
  else if constexpr (requires { typename T::value_type; }) 
    if constexpr (std::is_floating_point_v<typename T::value_type>)
      return std::numeric_limits<typename T::value_type>::epsilon();
  return 0.; // dummy
}

export
template<typename T>
constexpr bool is_approx_equal(const T& a, const T& b, std::optional<double> abs_tol = std::nullopt) 
{
  if constexpr (std::is_floating_point_v<T>) 
  {
    auto diff = a - b;
    if (diff < static_cast<T>(0)) 
      diff = -diff;
    if (abs_tol) 
      return diff < static_cast<T>(*abs_tol);
    return diff < (get_epsilon<T>() * static_cast<T>(2));
  } else 
  {
    if constexpr (requires { typename T::value_type; }) 
    {
      using value_t = typename T::value_type;
      if constexpr (std::is_floating_point_v<value_t> && 
                    requires(T r) { 
                      { r.real() } -> std::convertible_to<value_t>;
                      { r.imag() } -> std::convertible_to<value_t>; 
                    }) 
      {
        auto diff = a - b;
        auto real_diff = diff.real();
        auto imag_diff = diff.imag();
        if (real_diff < static_cast<value_t>(0)) 
          real_diff = -real_diff;
        if (imag_diff < static_cast<value_t>(0)) 
          imag_diff = -imag_diff;
        
        if (abs_tol) 
        {
          auto tol = static_cast<value_t>(*abs_tol);
          return (real_diff < tol) && (imag_diff < tol);
        }

        auto eps = get_epsilon<T>() * static_cast<value_t>(2);
        return (real_diff < eps) && (imag_diff < eps);
      }
    }
    return a == b;
  }
}

export
template<typename T>
constexpr bool is_negligible(const T& val, std::optional<double> abs_tol = std::nullopt) 
{ return is_approx_equal(val, T(0), abs_tol); }

export
template<ring_element_c_weak R = double,
         std::size_t N = 0>
struct polynomial_nttp
{
  using coefficient_t = R; // std::type_identity_t<R>?
  using index_t = decltype(N);

  index_t degree = N;
  std::array<coefficient_t, N + 1> coefficients{};

  // possibly allegedly std::fill may not perform this sometimes?
  // could be it's only a thing with std::vector
  constexpr polynomial_nttp() noexcept
  {
    stdr::fill(stdr::begin(coefficients),
               stdr::end(coefficients),
               coefficient_t{0});
  }

  constexpr polynomial_nttp(std::initializer_list<coefficient_t> li) noexcept
  { std::copy(li.begin(), li.end(), coefficients.begin()); }

  explicit constexpr
  polynomial_nttp(std::array<coefficient_t, N + 1> coeffs) noexcept
  : coefficients(coeffs)
  { }

  explicit constexpr polynomial_nttp(coefficient_t constant) noexcept
  {
    if constexpr (N > 0)
    {
      stdr::fill(stdr::begin(coefficients),
                 stdr::end(coefficients),
                 coefficient_t(0));
    }
    coefficients[0] = constant;
  }

  using begin_type = decltype(stdr::begin(coefficients));
  using end_type = decltype(stdr::end(coefficients));
  using cbegin_type = decltype(stdr::cbegin(coefficients));
  using cend_type = decltype(stdr::cend(coefficients));
  using rbegin_type = decltype(stdr::rbegin(coefficients));
  using rend_type = decltype(stdr::rend(coefficients));
  using crbegin_type = decltype(stdr::crbegin(coefficients));
  using crend_type = decltype(stdr::crend(coefficients));
  constexpr auto begin()
  { return std::forward<begin_type>(stdr::begin(coefficients)); }
  constexpr auto end()
  { return std::forward<end_type>(stdr::end(coefficients)); }
  constexpr auto begin() const
  { return std::forward<cbegin_type>(stdr::begin(coefficients)); }
  constexpr auto end() const
  { return std::forward<cend_type>(stdr::end(coefficients)); }
  constexpr auto cbegin() const
  { return std::forward<cbegin_type>(stdr::cbegin(coefficients)); }
  constexpr auto cend() const
  { return std::forward<cend_type>(stdr::cend(coefficients)); }
  constexpr auto rbegin()
  { return std::forward<rbegin_type>(stdr::rbegin(coefficients)); }
  constexpr auto rend()
  { return std::forward<rend_type>(stdr::rend(coefficients)); }
  constexpr auto crbegin() const
  { return std::forward<crbegin_type>(stdr::crbegin(coefficients)); }
  constexpr auto crend() const
  { return std::forward<crend_type>(stdr::crend(coefficients)); }

  constexpr coefficient_t operator[](index_t index) noexcept
  { return index <= N ? coefficients.at(index) : coefficient_t{0}; }

  constexpr coefficient_t operator[](index_t index) const noexcept
  { return index <= N ? coefficients.at(index) : coefficient_t{0}; }

  // Horner's method: p(x) = a₀ + x(a₁ + x(a₂ + x(...)))
  // O(n) multiplications, sequential dependency chain
  constexpr coefficient_t operator()(const coefficient_t x) const noexcept
  {
    coefficient_t result = coefficients[N];
    for (index_t i = N; i > 0; --i)
      result = result * x + coefficients[i - 1];
    return result;
  }

  //raw loops may be better in constexpr context
  constexpr polynomial_nttp operator-() const noexcept
  {
    polynomial_nttp negated_polynomial{};
    stdr::transform(stdr::cbegin(coefficients),
                    stdr::cend(coefficients),
                    stdr::begin(negated_polynomial),
                    std::negate());
    return negated_polynomial;
  };
};

// ============================================================
// Roots Types (Moved from roots module to break circular dependency)
// ============================================================

// Root with algebraic multiplicity
export template<typename K>
struct root_with_multiplicity
{
  K value;
  std::size_t multiplicity{1};

  constexpr bool operator==(const root_with_multiplicity&) const = default;
};

// Fixed-capacity roots container for compile-time and runtime use
export template<typename K, std::size_t max_roots>
struct roots_result
{
  std::array<root_with_multiplicity<K>, max_roots> data{};
  std::size_t count{0};

  // Iteration
  constexpr auto begin() const { return data.begin(); }
  constexpr auto end() const { return data.begin() + count; }
  constexpr auto begin() { return data.begin(); }
  constexpr auto end() { return data.begin() + count; }
  constexpr std::size_t size() const { return count; }
  constexpr bool empty() const { return count == 0; }
  constexpr const auto& operator[](std::size_t i) const { return data[i]; }
  constexpr auto& operator[](std::size_t i) { return data[i]; }

  // Mutation
  constexpr void push(root_with_multiplicity<K> r)
  {
    if (count < max_roots)
      data[count++] = r;
  }
  constexpr void push(K value, std::size_t mult = 1) { push({value, mult}); }

  // Append all roots from another result
  template<std::size_t M>
  constexpr void append(const roots_result<K, M>& other)
  {
    for (std::size_t i = 0; i < other.count && count < max_roots; ++i)
      data[count++] = other.data[i];
  }

  // Conversion to vector (for runtime convenience)
  operator std::vector<root_with_multiplicity<K>>() const
  {
    return std::vector<root_with_multiplicity<K>>{begin(), end()};
  }

  constexpr std::vector<root_with_multiplicity<K>> to_vector() const
  {
    return std::vector<root_with_multiplicity<K>>{begin(), end()};
  }

  // Extract just the root values
  constexpr std::array<K, max_roots> values() const
  {
    std::array<K, max_roots> vals{};
    for (std::size_t i = 0; i < count; ++i)
      vals[i] = data[i].value;
    return vals;
  }
};

} // end namespace lam::polynomial::univariate

// Export to lam::polynomial for convenient access
namespace lam::polynomial {
export using univariate::polynomial_nttp;
export using univariate::is_approx_equal;
export using univariate::is_negligible;
export using univariate::get_epsilon;
} // end namespace lam::polynomial

// Export to lam for the simplest access
namespace lam {
export using polynomial::univariate::polynomial_nttp;
export using polynomial::univariate::is_approx_equal;
export using polynomial::univariate::is_negligible;
export using polynomial::univariate::get_epsilon;
} // end namespace lam
