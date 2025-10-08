
#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace math
{

auto indexing_set = [](auto n) {
  return stdv::iota(decltype(n)(0), n);
};

auto indexing_set_from_to = [](auto m, auto n) {
  return stdv::iota(decltype(n)(m), n);
};

template<typename R = double, std::size_t N = 0>
struct polynomial_nttp
{
  using coefficient_t = R;
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
  polynomial_nttp(const std::array<coefficient_t, N + 1> coeffs) noexcept
  : coefficients(coeffs)
  { }

  constexpr auto begin()
  { return stdr::begin(coefficients); }
  constexpr auto end()
  { return stdr::end(coefficients); }
  constexpr auto begin() const
  { return stdr::begin(coefficients); }
  constexpr auto end() const
  { return stdr::end(coefficients); }
  constexpr auto cbegin() const
  { return stdr::cbegin(coefficients); }
  constexpr auto cend() const
  { return stdr::cend(coefficients); }

  constexpr coefficient_t operator[](const index_t index) noexcept
  { return index <= N ? coefficients[index] : coefficient_t{0}; }

  constexpr coefficient_t operator[](const index_t index) const noexcept
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

template<typename R = double, std::size_t M, std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_plus_q{};
    for (auto i : indexing_set(max_degree + 1))
      p_plus_q.coefficients[i] = p[i] + q[i];
    return p_plus_q;
  }();
}

template<typename R = double, std::size_t M, std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);

  return [&]() {
    polynomial_nttp<R, max_degree> p_minus_q{};
    for (auto i : indexing_set(max_degree + 1))
      p_minus_q.coefficients[i] = p[i] - q[i];
    return p_minus_q;
  }();
}

template<typename R = double, std::size_t M, std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, M>& p,
                         const polynomial_nttp<R, N>& q) noexcept
{
  return [&]() {
    polynomial_nttp<R, M + N> p_times_q{};
    for (auto k : indexing_set(M + N + 1))
      for (auto i : indexing_set(k))
        p_times_q.coefficients[k] += p[i] * q[k - i];
    return p_times_q;
  }();
}

template<typename R = double, std::size_t N>
constexpr auto operator*(const R& r,
                         const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N> r_times_p{};
    for (auto k : indexing_set(N + 1))
      r_times_p.coefficients[k] = r * p[k];
    return r_times_p;
  }();
}

template<typename R = double, std::size_t N>
constexpr auto norm([[maybe_unused]] const polynomial_nttp<R, N>& p) noexcept
{ return N; } //p.degree; }

template<typename R = double, std::size_t N>
constexpr auto leading(const polynomial_nttp<R, N>& p) noexcept
{ return p.coefficients[N]; }

template<typename R = double, std::size_t N>
constexpr auto make_monomial() noexcept
{
  polynomial_nttp<R, N> new_monomial{};
  new_monomial[N] = static_cast<R>(1);
  return new_monomial;
}

// R needs to be a field
template<typename R = double,
         std::size_t M,
         polynomial_nttp<R, M> a_of_x,
         std::size_t N,
         polynomial_nttp<R, N> b_of_x>
constexpr auto division_prototype()
noexcept
{
  auto comparatore = [](const auto& a, const auto& b_i) {
    auto abs_val = a - b_i;
    abs_val = abs_val > 0 ? abs_val : -abs_val;
    return abs_val < 1e-7 ? true : false;
  };
  constexpr auto b_is_not_zero = [&]() {
    bool is_zero = true;
    for (const auto& b_i : b_of_x)
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
    if (N <= M && b_is_not_zero)
    {
      auto deg_b = norm(b_of_x);
      auto lead_coeff_b = leading(b_of_x);
      while (remainder.size() - 1 >= deg_b)
      {
        auto s = remainder.back() / lead_coeff_b;
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
    }
    std::size_t quotient_size = quotient.size() - 1;
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

// R needs to be a field and the cast R(i) needs to make sense
template<typename R = double, std::size_t N>
constexpr auto
derivative(const polynomial_nttp<R, N>& p) noexcept
{
  if constexpr (N == 0)
    return polynomial_nttp<R, 0>{};
  else {
    return [&]() {
      polynomial_nttp<R, N - 1> ddxp{};
      for (auto i : indexing_set(N))
        ddxp.coefficients[i] = (static_cast<R>(i) + 1) * p[i + 1];
      return ddxp;
    }();
  }
}

/*
 * this function returns the antiderivative of p for which a_0 = 0
 * R needs to be a field and the cast R(i) needs to make sense
 */
template<typename R = double, std::size_t N>
constexpr polynomial_nttp<R, N + 1>
antiderivative(const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N + 1> antiderivative_of_p{};
    for (auto i : indexing_set_from_to(1, N + 1 + 1))
      antiderivative_of_p.coefficients[i] = 1 / static_cast<R>(i) * p[i - 1];
    return antiderivative_of_p;
  }();
}

} // end namespace math
