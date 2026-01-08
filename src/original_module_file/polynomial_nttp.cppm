/*
 *  polynomial_nttp.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

module;

import std;
import lam.experimental.concepts;
using namespace lam::experimental::concepts;

export module polynomial_nttp;

namespace stdr = std::ranges;
namespace stdv = std::views;

export namespace math_nttp
{

namespace polynomial
{

auto indexing_set = [](auto n) {
  return stdv::iota(decltype(n)(0), n);
};

auto indexing_set_from_to = [](auto m, auto n) {
  return stdv::iota(decltype(n)(m), n);
};

template<experimental::concepts::ring_element_c_weak R = double,
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
  }

};

// enable syntax for adding polynomials in R[X]
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_plus_q{};
    for (auto&& i : indexing_set(max_degree + 1))
      p_plus_q.coefficients[i] = p[i] + q[i];
    return p_plus_q;
  }();
}

// add a constant on the left
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator+(R&& r,
                         const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N> r_plus_p = p;
    r_plus_p.coefficients[0] = r + r_plus_p.coefficients[0];
    return r_plus_p;
  }();
}

// add a constant on the right
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, N>& p,
                         R&& r) noexcept
{
  return [&]() {
    polynomial_nttp<R, N> p_plus_r = p;
    p_plus_r.coefficients[0] = p_plus_r.coefficients[0] + r;
    return p_plus_r;
  }();
}
// enable syntax for subtracting polynomials in R[X]
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_minus_q{};
    for (auto&& i : indexing_set(max_degree + 1))
      p_minus_q.coefficients[i] = p[i] - q[i];
    return p_minus_q;
  }();
}

// subtract a constant on the left
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator-(R&& r,
                         const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N> r_minus_p = -p;
    r_minus_p.coefficients[0] += r; // R needs +=
    return r_minus_p;
  }();
}

// subtract a constant on the right
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, N>& p,
                         R&& r) noexcept
{
  return [&]() {
    polynomial_nttp<R, N> p_minus_r = p;
    p_minus_r.coefficients[0] = p_minus_r.coefficients[0] - r;
    return p_minus_r;
  }();
}
// enable syntax for multiplying polynomials in R[X]
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  return [&]() {
    polynomial_nttp<R, M + N> p_times_q{};
    for (auto&& k : indexing_set(M + N + 1))
      for (auto&& i : indexing_set(k + 1))
        p_times_q.coefficients[k] += p[i] * q[k - i];
    return p_times_q;
  }();
}

// multiply a constant on the left
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator*(R&& r,
                         const polynomial_nttp<R, N>& p) noexcept
{
  polynomial_nttp<R, N> r_times_p{};
  stdr::transform(stdr::cbegin(p),
                  stdr::cend(p),
                  stdr::begin(r_times_p),
                  [&](auto&& p_i) {
                    return r * p_i;
                  });
  return r_times_p;
}

// multiply a constant on the right
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, N>& p,
                         R&& r) noexcept
{
  polynomial_nttp<R, N> p_times_r{};
  stdr::transform(stdr::cbegin(p),
                  stdr::cend(p),
                  stdr::begin(p_times_r),
                  [&](auto&& p_i) {
                    return p_i * r;
                  });
  return p_times_r;
}

/*
 *  the norm of a polynomial is its degree
 *  a `polynomial_nttp` can be of degree `N > 0` with all zero coefficients
 */
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr auto norm([[maybe_unused]] const polynomial_nttp<R, N>& p) noexcept
{ return N; }
/* returns the leading coefficient, even if it is zero! */
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr R leading(const polynomial_nttp<R, N>& p) noexcept
{ return p.coefficients[N]; }
/* returns a new monomial of degree `N` */
template<experimental::concepts::ring_element_c_weak R = double,
         std::size_t N>
constexpr polynomial_nttp<R, N> make_monomial() noexcept
{
  polynomial_nttp<R, N> new_monomial{};
  new_monomial.coefficients[N] = static_cast<R>(1);
  return new_monomial;
}

/*
 *  division_prototype()
 *    takes several NTTPs:
 *      `std::size_t M`                 – the degree of the first polynomial
 *      `polynomial_nttp<R, M> a_of_x`  – polynomial a(x) of degree M in R[X]
 *      `std::size_t N`                 – the degree of the second polynomial
 *      `polynomial_nttp<R, N> b_of_x`  – polynomial b(x) of degree N in R[X]
 *
 *    returns a `std::pair` composed of the quotient `q` and the remainder `r`
 *      the exact return type is determined inside the body of the function
 *      in particular, it is determined within the immediately invoked lambda
 *      `sizes_and_oversized_arrays_q_and_r`, a `std::pair` of `std::pair`s
 *
 *    the quotient `q` is identically `0`, and
 *    the remainder `r` is identically `a_of_x` if
 *      - `b_of_x` is identically `0`, or
 *      - `N > M`, i.e., `norm(b_of_x) > norm(a_of_x)`
 *
 *    division_prototype() will trigger a compile-time failure if one chooses
 *      to create (directly or inadvertently) an `N` degree polynomial with
 *      0 as leading coefficient; the author sees this as a drawback of the
 *      current version and will introduce better behavior in the future
\*                                                                          */
template<experimental::concepts::field_element_c_weak R = double,
         std::size_t M,
         polynomial_nttp<R, M> a_of_x,
         std::size_t N,
         polynomial_nttp<R, N> b_of_x>
constexpr auto division_prototype()
noexcept
{
  constexpr auto comparatore = [&](auto&& a, auto&& b_i) {
    auto abs_val = a - b_i;
    abs_val = abs_val > 0 ? abs_val : -abs_val;
    return abs_val < 1e-7 ? true : false;
  };
  constexpr auto b_is_not_zero = [&]() {
    bool is_zero = true;
    for (auto&& b_i : b_of_x)
      if (!comparatore(static_cast<R>(0), b_i))
      {
        is_zero = false;
        break;
      }
    return !is_zero;
  }();

  constexpr auto sizes_and_oversized_arrays_q_and_r = [&]() {
    std::array<R, M + 1> oversized_quotient{};
    std::array<R, M + 1> oversized_remainder{};
    std::vector<R> quotient{};
    std::vector<R> remainder{};
    for (std::size_t i = 0; i < a_of_x.coefficients.size(); ++i)
      remainder.push_back(a_of_x.coefficients[i]);
    std::size_t quotient_size = 0;
    if (N <= M && b_is_not_zero)
    {
      auto deg_b = norm(b_of_x);
      auto lead_coeff_b = leading(b_of_x); // will fail at compile time if = 0
      while (remainder.size() - 1 >= deg_b)
      {
        auto s = remainder.back() / lead_coeff_b; // see 3 lines above
        quotient.push_back(s);
        std::size_t i_offset = remainder.size() - deg_b - 1;
        for (std::size_t i = 0; i <= N; ++i)
          remainder[i + i_offset] = remainder[i + i_offset]
                                    - (s * b_of_x.coefficients[i]);
        remainder.pop_back();
      }
      if (comparatore(static_cast<R>(0), remainder.back()))
        remainder.pop_back();
      std::reverse(quotient.begin(), quotient.end());
      quotient_size = quotient.size() - 1;
    } else // will dereference a null without something like this
    { quotient.push_back(static_cast<R>(0)); }
    for (std::size_t i = 0; i <= quotient_size; ++i)
      oversized_quotient[i] = quotient[i];
    std::size_t remainder_size = remainder.size() - 1;
    for (std::size_t i = 0; i <= remainder_size; ++i)
      oversized_remainder[i] = remainder[i];
    return std::make_pair(std::make_pair(quotient_size,
                                         remainder_size),
                          std::make_pair(oversized_quotient,
                                         oversized_remainder));
  }();

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.first> q{};
  for (std::size_t i = 0; i <= norm(q); ++i)
    q.coefficients[i] = sizes_and_oversized_arrays_q_and_r.second.first[i];

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.second> r{};
  for (std::size_t i = 0; i <= norm(r); ++i)
    r.coefficients[i] = sizes_and_oversized_arrays_q_and_r.second.second[i];

  return std::make_pair(q, r);
}

/*  R needs to be a field and the cast R(i) needs to make sense
 *
 *  return type is `auto` because:
 *    the derivative of a constant is 0, and as a polynomial is of degree 0
 *    the derivative of a degree `N > 0` polynomial is of degree `N - 1`
 *    `auto` and `if constexpr` are how we choose to avoid `0 - 1` underflow
 */
template<experimental::concepts::field_element_c_weak R = double,
         std::size_t N>
constexpr auto derivative(const polynomial_nttp<R, N>& p) noexcept
{
  if constexpr (N == 0)
    return polynomial_nttp<R, 0>{};
  else
  {
    polynomial_nttp<R, N - 1> ddxp{};
    auto ids = indexing_set(N);
    stdr::transform(stdr::cbegin(ids),
                    stdr::cend(ids),
                    stdr::begin(ddxp),
                    [&](auto&& index) {
                      return static_cast<R>(index + 1) * p[index + 1];
                    });
    return ddxp;
  }
}

/*
 *  this function returns the antiderivative of p for which a_0 = 0
 *  R needs to be a field and the cast R(i) needs to make sense
 */
template<experimental::concepts::field_element_c_weak R = double,
         std::size_t N>
constexpr polynomial_nttp<R, N + 1>
antiderivative(const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N + 1> antiderivative_of_p{};
    auto ids = indexing_set_from_to(1, N + 1 + 1);
    for (auto i : ids)
      antiderivative_of_p.coefficients[i] = 1 / static_cast<R>(i) * p[i - 1];
    return antiderivative_of_p;
  }();
}
} // end namespace polynomial
  using namespace polynomial; // lol rip gcc c++20 modules implementation
} // end namespace math_nttp
