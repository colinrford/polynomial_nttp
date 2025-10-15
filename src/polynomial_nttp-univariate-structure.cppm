/*
 *  polynomial_nttp-univariate-structure.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module polynomial_nttp:univariate.structure;

import std;
import experimental.concepts;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace math_nttp
{

namespace polynomial
{

namespace univariate
{

template<typename R>
concept ring_element_c_weak = experimental::concepts::ring_element_c_weak<R>;

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

  constexpr polynomial_nttp(std::initializer_list<R> li) noexcept
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

  constexpr coefficient_t operator[](index_t index) noexcept
  { return index <= N ? coefficients[index] : coefficient_t{0}; }

  constexpr coefficient_t operator[](index_t index) const noexcept
  { return index <= N ? coefficients[index] : coefficient_t{0}; }

  constexpr coefficient_t operator()(const coefficient_t x) const noexcept
  { // obviously insufficient :-)
    coefficient_t value = coefficient_t{0};
    for (index_t i = 0; i != N + 1; ++i)
    {
      coefficient_t x_power = i == 0 ? 1 : x;
      for (index_t j = 1; j < i; ++j)
        x_power *= x;
      value += coefficients[i] * x_power;
    }
    return value;
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

} // end namespace univariate
} // end namespace polynomial
  export using namespace polynomial::univariate;
  //using polynomial::univariate::polynomial_nttp;
} // end namespace math_nttp
