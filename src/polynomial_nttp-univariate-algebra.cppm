/*
 *  polynomial_nttp-univariate-algebra.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module polynomial_nttp:univariate.algebra;

import std;
import experimental.concepts;
import :univariate.structure;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace math_nttp
{

namespace polynomial
{

namespace univariate
{

const auto indexing_set = [](auto n) {
  return stdv::iota(static_cast<decltype(n)>(0), n);
};

const auto indexing_set_from_to = [](auto m, auto n) {
  return stdv::iota(static_cast<decltype(n)>(m), n);
};

// ring_element_c_weak is defined in :univariate.structure
template<typename K>
concept field_element_c_weak = experimental::concepts::field_element_c_weak<K>;

// enable syntax for adding polynomials in R[X]
export
template<ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q)
noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_plus_q{};
    for (auto&& i : indexing_set(max_degree + 1))
      p_plus_q.coefficients.at(i) = p[i] + q[i];
    return p_plus_q;
  }();
}

// add a constant on the left
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator+(const R& r, const polynomial_nttp<R, N>& p)
noexcept
{
  return [&]() {
    polynomial_nttp<R, N> r_plus_p = p;
    r_plus_p.coefficients[0] = r + r_plus_p.coefficients[0];
    return r_plus_p;
  }();
}

// add a constant on the right
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, N>& p, const R& r)
noexcept
{
  return [&]() {
    polynomial_nttp<R, N> p_plus_r = p;
    p_plus_r.coefficients[0] = p_plus_r.coefficients[0] + r;
    return p_plus_r;
  }();
}

// enable syntax for subtracting polynomials in R[X]
export
template<ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q)
noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_minus_q{};
    for (auto&& i : indexing_set(max_degree + 1))
      p_minus_q.coefficients.at(i) = p[i] - q[i];
    return p_minus_q;
  }();
}

// subtract a constant on the left
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator-(const R& r, const polynomial_nttp<R, N>& p)
noexcept
{
  return [&]() {
    polynomial_nttp<R, N> r_minus_p = -p;
    r_minus_p.coefficients[0] += r; // R needs +=
    return r_minus_p;
  }();
}

// subtract a constant on the right
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, N>& p, const R& r)
noexcept
{
  return [&]() {
    polynomial_nttp<R, N> p_minus_r = p;
    p_minus_r.coefficients[0] = p_minus_r.coefficients[0] - r;
    return p_minus_r;
  }();
}

// enable syntax for multiplying polynomials in R[X]
export
template<ring_element_c_weak R = double,
         std::size_t M,
         std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q)
noexcept
{
  return [&]() {
    polynomial_nttp<R, M + N> p_times_q{};
    for (auto&& k : indexing_set(M + N + 1))
      for (auto&& i : indexing_set(k + 1))
        p_times_q.coefficients.at(k) += p[i] * q[k - i];
    return p_times_q;
  }();
}

// multiply a constant on the left
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator*(const R& r, const polynomial_nttp<R, N>& p)
noexcept
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
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, N>& p, const R& r)
noexcept
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
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr auto norm([[maybe_unused]] const polynomial_nttp<R, N>& p)
noexcept
{ return N; }
/* returns the leading coefficient, even if it is zero! */
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr R leading(const polynomial_nttp<R, N>& p)
noexcept
{ return p.coefficients[N]; }
/* returns the leading coefficient, even if it is zero! */
//export
template<ring_element_c_weak R = double,
         std::size_t N,
         polynomial_nttp<R, N> p_of_x>
constexpr R leading()
noexcept
{
  if constexpr (p_of_x.coefficients[N] != static_cast<R>(0)) // placeholder
  {
    return leading(p_of_x);
  } else
  {
    return p_of_x.coefficients[N];
  }
}
/* returns a new monomial of degree `N` */
export
template<ring_element_c_weak R = double,
         std::size_t N>
constexpr polynomial_nttp<R, N> make_monomial()
noexcept
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
export
template<field_element_c_weak R = double,
         std::size_t M,
         polynomial_nttp<R, M> a_of_x,
         std::size_t N,
         polynomial_nttp<R, N> b_of_x>
constexpr auto division_prototype()
// noexcept // don't think this can be marked noexcept!
{
  constexpr auto comparatore = [&](auto&& a, auto&& b_i) {
    auto abs_val = a - b_i;
    abs_val = abs_val > 0 ? abs_val : -abs_val;
    return abs_val < 1e-7 ? true : false; // NOLINT
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
    //remainder.reserve(a_of_x.coefficients.size()); // have not compared
    stdr::copy(a_of_x.coefficients.begin(),
               a_of_x.coefficients.end(),
               std::back_inserter(remainder));
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
                                    - (s * b_of_x.coefficients.at(i));
        remainder.pop_back();
      }
      if (comparatore(static_cast<R>(0), remainder.back()))
        remainder.pop_back();
      std::reverse(quotient.begin(), quotient.end());
      quotient_size = quotient.size() - 1;
    } else // will dereference a null without something like this
    { quotient.push_back(static_cast<R>(0)); }
    for (std::size_t i = 0; i <= quotient_size; ++i)
      oversized_quotient.at(i) = quotient[i];
    std::size_t remainder_size = remainder.size() - 1;
    for (std::size_t i = 0; i <= remainder_size; ++i)
      oversized_remainder.at(i) = remainder[i];
    return std::make_pair(std::make_pair(quotient_size,
                                         remainder_size),
                          std::make_pair(oversized_quotient,
                                         oversized_remainder));
  }();

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.first> q{};
  for (std::size_t i = 0; i <= norm(q); ++i)
    q.coefficients.at(i) = sizes_and_oversized_arrays_q_and_r.second.first.at(i);

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.second> r{};
  for (std::size_t i = 0; i <= norm(r); ++i)
    r.coefficients.at(i) = sizes_and_oversized_arrays_q_and_r.second.second.at(i);

  return std::make_pair(q, r);
}

/*  R needs to model a field and static_cast<R>(i) needs to make sense
 *
 *  return type is `auto` because:
 *    the derivative of a constant is 0, and as a polynomial is of degree 0
 *    the derivative of a degree `N > 0` polynomial is of degree `N - 1`
 *    `auto` and `if constexpr` are how we choose to avoid `0 - 1` underflow
 */
export
template<field_element_c_weak R = double,
         std::size_t N>
constexpr auto derivative(const polynomial_nttp<R, N>& p)
noexcept
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
 *  R needs to model a field and static_cast<R>(i) needs to make sense
 */
export
template<field_element_c_weak R = double,
         std::size_t N>
constexpr polynomial_nttp<R, N + 1>
antiderivative(const polynomial_nttp<R, N>& p)
noexcept
{
  return [&]() {
    polynomial_nttp<R, N + 1> antiderivative_of_p{};
    auto ids = indexing_set_from_to(1, N + 1 + 1);
    for (auto&& i : ids)
      antiderivative_of_p.coefficients.at(i) = 1 / static_cast<R>(i) * p[i - 1];
    return antiderivative_of_p;
  }();
}

} // end namespace univariate
} // end namespace polynomial
  export using namespace polynomial::univariate;
} // end namespace math_nttp
