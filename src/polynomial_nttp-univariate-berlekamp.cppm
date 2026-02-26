/*
 *  polynomial_nttp-univariate-berlekamp.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 *
 *  Berlekamp's algorithm for polynomial factorization over finite fields
 */

export module lam.polynomial_nttp:univariate.berlekamp;

import std;
import lam.concepts;
import :univariate.structure;
import :univariate.algebra;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace lam::polynomial::univariate::berlekamp
{

// ============================================================
// Polynomial Degree Utilities
// ============================================================

// Compute effective degree (ignoring trailing zeros)
template<field_element_c_weak K, std::size_t N>
constexpr std::size_t effective_degree(const polynomial_nttp<K, N>& p)
{
  for (std::size_t i = N; i > 0; --i)
    if (!is_negligible(p[i]))
      return i;
  return 0;
}

template<field_element_c_weak K, std::size_t N>
constexpr bool is_zero_poly(const polynomial_nttp<K, N>& p)
{
  for (std::size_t i = 0; i <= N; ++i)
    if (!is_negligible(p[i]))
      return false;
  return true;
}

// ============================================================
// Polynomial GCD (Euclidean Algorithm)
// ============================================================

// GCD for polynomials of the same static degree
// Returns monic GCD (leading coefficient = 1)
export template<field_element_c_weak K, std::size_t N>
constexpr auto poly_gcd(polynomial_nttp<K, N> a, polynomial_nttp<K, N> b) -> polynomial_nttp<K, N>
{
  // Euclidean algorithm
  while (!is_zero_poly(b))
  {
    auto deg_a = effective_degree(a);
    auto deg_b = effective_degree(b);

    if (deg_a < deg_b)
      std::swap(a, b);

    polynomial_nttp<K, N> remainder = a;
    K lead_b = b[effective_degree(b)];

    while (effective_degree(remainder) >= effective_degree(b) && !is_zero_poly(remainder))
    {
      auto deg_r = effective_degree(remainder);
      auto deg_divisor = effective_degree(b);

      if (deg_r < deg_divisor)
        break;

      K coeff = remainder[deg_r] / lead_b;
      std::size_t shift = deg_r - deg_divisor;

      for (std::size_t i = 0; i <= deg_divisor; ++i)
        if (i + shift <= N)
          remainder.coefficients[i + shift] = remainder[i + shift] - coeff * b[i];
    }

    a = b;
    b = remainder;
  }

  auto deg = effective_degree(a);
  if (deg > 0 || !is_negligible(a[0]))
  {
    K lead = a[deg];
    if (!is_negligible(lead))
      for (std::size_t i = 0; i <= deg; ++i)
        a.coefficients[i] = a[i] / lead;
  }

  return a;
}

// ============================================================
// Power of x modulo f(x)
// ============================================================

template<field_element_c_weak K, std::size_t N>
constexpr auto power_mod(const polynomial_nttp<K, N>& base, std::size_t exp, const polynomial_nttp<K, N>& modulus)
  -> polynomial_nttp<K, N>
{
  auto reduce_mod = [&modulus](const auto& poly) -> polynomial_nttp<K, N> {
    // Copy coefficients to dynamic vector to handle high degrees properly
    std::vector<K> rem(poly.coefficients.begin(), poly.coefficients.end());

    K lead_mod = modulus[effective_degree(modulus)];
    std::size_t deg_mod = effective_degree(modulus);

    while (true)
    {
      // Find effective degree of remainder
      std::size_t deg_r = 0;
      bool is_zero = true;
      for (std::size_t i = rem.size(); i-- > 0;)
      {
        if (!is_negligible(rem[i]))
        {
          deg_r = i;
          is_zero = false;
          break;
        }
      }

      if (is_zero || deg_r < deg_mod)
        break;

      K coeff = rem[deg_r] / lead_mod;
      std::size_t shift = deg_r - deg_mod;

      for (std::size_t i = 0; i <= deg_mod; ++i)
      {
        rem[i + shift] = rem[i + shift] - coeff * modulus[i];
      }
    }

    polynomial_nttp<K, N> result{};
    for (std::size_t i = 0; i <= N && i < rem.size(); ++i)
      result.coefficients[i] = rem[i];
    return result;
  };

  polynomial_nttp<K, N> result{};
  result.coefficients[0] = K(1);
  polynomial_nttp<K, N> current = base;

  while (exp > 0)
  {
    if (exp & 1)
    {
      auto product = result * current;
      result = reduce_mod(product);
    }

    {
      auto product = current * current;
      current = reduce_mod(product);
    }

    exp >>= 1;
  }

  return result;
}

// ============================================================
// Berlekamp Matrix Construction
// ============================================================

// Build Berlekamp matrix B where B[i] = coefficients of x^(i*p) mod f
// P = field characteristic
export template<field_element_c_weak K, std::size_t P, std::size_t N>
constexpr auto build_berlekamp_matrix(const polynomial_nttp<K, N>& f) -> std::array<std::array<K, N>, N>
{
  std::array<std::array<K, N>, N> B{};

  polynomial_nttp<K, N> x_poly{};
  x_poly.coefficients[1] = K(1);

  for (std::size_t i = 0; i < N; ++i)
  { // Compute x^(i*P) mod f
    auto x_ip = power_mod(x_poly, i * P, f);

    for (std::size_t j = 0; j < N; ++j)
      B[i][j] = x_ip[j];
  }

  return B;
}

// ============================================================
// Gaussian Elimination for Null Space
// ============================================================

export template<field_element_c_weak K, std::size_t N>
constexpr auto berlekamp_null_space(std::array<std::array<K, N>, N> B)
  -> std::pair<std::array<std::array<K, N>, N>, std::size_t>
{
  for (std::size_t i = 0; i < N; ++i)
    B[i][i] = B[i][i] - K(1);

  std::array<std::size_t, N> pivot_col{};
  std::size_t rank = 0;

  for (std::size_t col = 0; col < N && rank < N; ++col)
  {
    std::size_t pivot_row = rank;
    while (pivot_row < N && is_negligible(B[pivot_row][col]))
      ++pivot_row;

    if (pivot_row == N)
      continue;

    if (pivot_row != rank)
      std::swap(B[rank], B[pivot_row]);

    pivot_col[rank] = col;

    K inv = K(1) / B[rank][col];
    for (std::size_t j = col; j < N; ++j)
      B[rank][j] = B[rank][j] * inv;

    for (std::size_t i = 0; i < N; ++i)
    {
      if (i != rank && !is_negligible(B[i][col]))
      {
        K factor = B[i][col];
        for (std::size_t j = col; j < N; ++j)
          B[i][j] = B[i][j] - factor * B[rank][j];
      }
    }

    ++rank;
  }

  std::size_t null_dim = N - rank;

  std::array<std::array<K, N>, N> null_basis{};
  std::size_t basis_idx = 0;

  null_basis[0][0] = K(1);
  if (null_dim > 0)
    basis_idx = 1;

  std::array<bool, N> is_pivot{};
  for (std::size_t i = 0; i < rank; ++i)
    is_pivot[pivot_col[i]] = true;

  for (std::size_t free_col = 1; free_col < N && basis_idx < null_dim; ++free_col)
  {
    if (is_pivot[free_col])
      continue;

    std::array<K, N> basis_vec{};
    basis_vec[free_col] = K(1);

    for (std::size_t i = rank; i-- > 0;)
    {
      K sum = K(0);
      for (std::size_t j = pivot_col[i] + 1; j < N; ++j)
        sum = sum + B[i][j] * basis_vec[j];
      basis_vec[pivot_col[i]] = K(0) - sum;
    }

    null_basis[basis_idx++] = basis_vec;
  }

  return {null_basis, null_dim};
}

// ============================================================
// Berlekamp Factorization
// ============================================================

export template<field_element_c_weak K, std::size_t P, std::size_t N>
constexpr auto berlekamp_factor(const polynomial_nttp<K, N>& f, K zero, K one)
  -> std::pair<std::array<polynomial_nttp<K, N>, N>, std::size_t>
{
  std::array<polynomial_nttp<K, N>, N> factors{};
  std::size_t factor_count = 0;

  auto deg = effective_degree(f);
  if (deg <= 1)
  {
    factors[0] = f;
    return {factors, 1};
  }

  auto B = build_berlekamp_matrix<K, P, N>(f);

  auto [null_basis, null_dim] = berlekamp_null_space(B);

  if (null_dim == 1)
  { // f is irreducible
    factors[0] = f;
    return {factors, 1};
  }

  std::array<polynomial_nttp<K, N>, N> to_factor{};
  std::size_t queue_size = 1;
  to_factor[0] = f;

  for (std::size_t k = 1; k < null_dim && queue_size > 0; ++k)
  { // Convert basis vector to polynomial h(x)
    polynomial_nttp<K, N> h{};
    for (std::size_t i = 0; i < N; ++i)
      h.coefficients[i] = null_basis[k][i];

    std::array<polynomial_nttp<K, N>, N> new_queue{};
    std::size_t new_queue_size = 0;

    for (std::size_t q = 0; q < queue_size; ++q)
    {
      auto& g = to_factor[q];
      auto g_deg = effective_degree(g);

      if (g_deg <= 1)
      {
        factors[factor_count++] = g;
        continue;
      }

      bool split = false;

      K c = zero;
      for (std::size_t c_val = 0; c_val < P; ++c_val)
      {
        polynomial_nttp<K, N> h_minus_c = h;
        h_minus_c.coefficients[0] = h[0] - c;

        auto d = poly_gcd(g, h_minus_c);
        auto d_deg = effective_degree(d);

        if (d_deg > 0 && d_deg < g_deg)
        { // Found a non-trivial factor!
          new_queue[new_queue_size++] = d;

          polynomial_nttp<K, N> quotient{};
          polynomial_nttp<K, N> remainder = g;
          K lead_d = d[d_deg];

          while (effective_degree(remainder) >= d_deg && !is_zero_poly(remainder))
          {
            auto deg_r = effective_degree(remainder);
            if (deg_r < d_deg)
              break;

            K coeff = remainder[deg_r] / lead_d;
            quotient.coefficients[deg_r - d_deg] = coeff;

            for (std::size_t i = 0; i <= d_deg; ++i)
            {
              if (i + deg_r - d_deg <= N)
                remainder.coefficients[i + deg_r - d_deg] = remainder[i + deg_r - d_deg] - coeff * d[i];
            }
          }

          new_queue[new_queue_size++] = quotient;
          split = true;
          break;
        }
        c = c + one;
      }

      if (!split)
        new_queue[new_queue_size++] = g;
    }

    to_factor = new_queue;
    queue_size = new_queue_size;
  }

  for (std::size_t i = 0; i < queue_size; ++i)
    factors[factor_count++] = to_factor[i];

  return {factors, factor_count};
}

// ============================================================
// Find Roots via Berlekamp Factorization
// ============================================================

export template<field_element_c_weak K, std::size_t P, std::size_t N>
constexpr auto roots_berlekamp(const polynomial_nttp<K, N>& f, K zero, K one) -> roots_result<K, N>
{
  roots_result<K, N> result;

  // Check for inseparable polynomial: f'(x) == 0
  // This implies f(x) = g(x^P)
  auto deg_f = effective_degree(f);
  if (deg_f >= P)
  {
    auto df = derivative(f);
    if (is_zero_poly(df))
    {
      polynomial_nttp<K, N> g{};
      for (std::size_t i = 0; i * P <= deg_f; ++i)
      {
        g.coefficients[i] = f[i * P];
      }

      auto sub_roots = roots_berlekamp<K, P, N>(g, zero, one);

      for (std::size_t i = 0; i < sub_roots.count; ++i)
      {
        result.push(sub_roots[i].value, sub_roots[i].multiplicity * P);
      }
      return result;
    }
  }

  auto [factors, count] = berlekamp_factor<K, P, N>(f, zero, one);

  for (std::size_t i = 0; i < count; ++i)
  {
    auto& factor = factors[i];
    if (effective_degree(factor) == 1)
    { // Linear factor: a*x + b = 0 => x = -b/a
      K a = factor[1];
      K b = factor[0];
      if (!is_negligible(a))
      {
        K root = (zero - b) / a;
        result.push(root, 1);
      }
    }
  }

  return result;
}

} // namespace lam::polynomial::univariate::berlekamp

// Export symbols to parent namespace
namespace lam::polynomial
{
export using lam::polynomial::univariate::berlekamp::poly_gcd;
export using lam::polynomial::univariate::berlekamp::berlekamp_factor;
export using lam::polynomial::univariate::berlekamp::roots_berlekamp;
} // namespace lam::polynomial
