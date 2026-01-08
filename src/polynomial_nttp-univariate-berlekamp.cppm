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
// Finite Field Concept
// ============================================================

// Note: We generate field elements on-the-fly
// Field size P is passed as a template parameter to avoid
// requiring a ::characteristic() method on the field type

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

// Check if polynomial is zero
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
export 
template<field_element_c_weak K, std::size_t N>
constexpr auto poly_gcd(polynomial_nttp<K, N> a, polynomial_nttp<K, N> b) -> polynomial_nttp<K, N>
{
  // Euclidean algorithm
  while (!is_zero_poly(b))
  {
    // Compute a mod b using division
    // We need runtime division here since degrees vary
    auto deg_a = effective_degree(a);
    auto deg_b = effective_degree(b);

    if (deg_a < deg_b)
      std::swap(a, b);

    // Manual polynomial division for remainder
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

  // Make monic (leading coefficient = 1)
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

// Compute x^exp mod f(x) using repeated squaring
template<field_element_c_weak K, std::size_t N>
constexpr auto power_mod(const polynomial_nttp<K, N>& base, std::size_t exp, const polynomial_nttp<K, N>& modulus)
  -> polynomial_nttp<K, N>
{
  // Helper lambda to reduce a polynomial mod modulus
  auto reduce_mod = [&modulus](const auto& poly) -> polynomial_nttp<K, N> {
    polynomial_nttp<K, N> result{};
    
    // Copy lower-degree coefficients
    auto poly_size = poly.coefficients.size();
    for (std::size_t i = 0; i < poly_size && i <= N; ++i)
      result.coefficients[i] = poly[i];
    
    K lead_mod = modulus[effective_degree(modulus)];
    auto deg_mod = effective_degree(modulus);
    
    // First reduce any coefficients above N (from the 2N product)
    for (std::size_t i = poly_size - 1; i > N; --i) {
      if (is_negligible(poly[i])) continue;
      // Reduce x^i by substituting x^(deg_mod) = -(other terms)/lead
      K coeff = poly[i] / lead_mod;
      std::size_t shift = i - deg_mod;
      for (std::size_t j = 0; j < deg_mod; ++j) {
        if (j + shift <= N)
          result.coefficients[j + shift] = result[j + shift] - coeff * modulus[j];
      }
    }
    
    // Now reduce within result
    while (effective_degree(result) >= deg_mod && !is_zero_poly(result))
    {
      auto deg_r = effective_degree(result);
      if (deg_r < deg_mod)
        break;

      K coeff = result[deg_r] / lead_mod;
      std::size_t shift = deg_r - deg_mod;

      for (std::size_t i = 0; i <= deg_mod; ++i)
        if (i + shift <= N)
          result.coefficients[i + shift] = result[i + shift] - coeff * modulus[i];
    }
    
    return result;
  };

  polynomial_nttp<K, N> result{};
  result.coefficients[0] = K(1); // Start with 1

  polynomial_nttp<K, N> current = base;

  while (exp > 0)
  {
    if (exp & 1)
    {
      // result = result * current mod modulus
      auto product = result * current;  // This is polynomial_nttp<K, 2N>
      result = reduce_mod(product);
    }

    // current = current^2 mod modulus
    {
      auto product = current * current;  // This is polynomial_nttp<K, 2N>
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
export 
template<field_element_c_weak K, std::size_t P, std::size_t N>
constexpr auto build_berlekamp_matrix(const polynomial_nttp<K, N>& f) -> std::array<std::array<K, N>, N>
{
  std::array<std::array<K, N>, N> B{};

  // x as a polynomial
  polynomial_nttp<K, N> x_poly{};
  x_poly.coefficients[1] = K(1);

  for (std::size_t i = 0; i < N; ++i)
  {
    // Compute x^(i*P) mod f
    auto x_ip = power_mod(x_poly, i * P, f);

    // Store coefficients as row i
    for (std::size_t j = 0; j < N; ++j)
      B[i][j] = x_ip[j];
  }

  return B;
}

// ============================================================
// Gaussian Elimination for Null Space
// ============================================================

// Compute null space basis of (B - I)
// Returns: (basis vectors, dimension of null space)
export 
template<field_element_c_weak K, std::size_t N>
constexpr auto berlekamp_null_space(std::array<std::array<K, N>, N> B)
  -> std::pair<std::array<std::array<K, N>, N>, std::size_t>
{
  // Subtract identity matrix
  for (std::size_t i = 0; i < N; ++i)
    B[i][i] = B[i][i] - K(1);

  // Gaussian elimination to reduced row echelon form
  std::array<std::size_t, N> pivot_col{};
  std::size_t rank = 0;

  for (std::size_t col = 0; col < N && rank < N; ++col)
  {
    // Find pivot row
    std::size_t pivot_row = rank;
    while (pivot_row < N && is_negligible(B[pivot_row][col]))
      ++pivot_row;

    if (pivot_row == N)
      continue; // No pivot in this column, it's a free variable

    // Swap rows
    if (pivot_row != rank)
      std::swap(B[rank], B[pivot_row]);

    pivot_col[rank] = col;

    // Scale pivot row to make leading coefficient 1
    K inv = K(1) / B[rank][col];
    for (std::size_t j = col; j < N; ++j)
      B[rank][j] = B[rank][j] * inv;

    // Eliminate all other entries in this column
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

  // Null space dimension = N - rank
  std::size_t null_dim = N - rank;

  // Extract null space basis from free columns
  std::array<std::array<K, N>, N> null_basis{};
  std::size_t basis_idx = 0;

  // First basis vector is always (1, 0, 0, ...) corresponding to constant poly
  null_basis[0][0] = K(1);
  if (null_dim > 0)
    basis_idx = 1;

  // Find free variables and construct basis vectors
  // Skip free_col=0 since null_basis[0] already covers it
  std::array<bool, N> is_pivot{};
  for (std::size_t i = 0; i < rank; ++i)
    is_pivot[pivot_col[i]] = true;

  for (std::size_t free_col = 1; free_col < N && basis_idx < null_dim; ++free_col)
  {
    if (is_pivot[free_col])
      continue;

    // Construct basis vector for this free variable
    std::array<K, N> basis_vec{};
    basis_vec[free_col] = K(1);

    // Back-substitute to find dependent variables
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

// P = field size (prime characteristic)
// zero and one are used to construct field elements (avoids requiring K(int))
export 
template<field_element_c_weak K, std::size_t P, std::size_t N>
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

  // Build Berlekamp matrix
  auto B = build_berlekamp_matrix<K, P, N>(f);

  // Find null space of (B - I)
  auto [null_basis, null_dim] = berlekamp_null_space(B);

  if (null_dim == 1)
  { // f is irreducible
    factors[0] = f;
    return {factors, 1};
  }

  // Work queue: polynomials still to factor
  std::array<polynomial_nttp<K, N>, N> to_factor{};
  std::size_t queue_size = 1;
  to_factor[0] = f;

  // For each null space basis vector (starting from index 1)
  for (std::size_t k = 1; k < null_dim && queue_size > 0; ++k)
  { // Convert basis vector to polynomial h(x)
    polynomial_nttp<K, N> h{};
    for (std::size_t i = 0; i < N; ++i)
      h.coefficients[i] = null_basis[k][i];

    // Try to split each polynomial in queue
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

      // Try gcd(g, h - c) for each c in the field
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

          // Compute g / d
          // For simplicity, divide using quotient
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
        c = c + one;  // Increment to next field element
      }

      if (!split)
        new_queue[new_queue_size++] = g;
    }

    to_factor = new_queue;
    queue_size = new_queue_size;
  }

  // Move remaining queue items to factors
  for (std::size_t i = 0; i < queue_size; ++i)
    factors[factor_count++] = to_factor[i];

  return {factors, factor_count};
}

// ============================================================
// Find Roots via Berlekamp Factorization
// ============================================================

// P = field size (prime characteristic)
// zero and one are used to construct field elements
export 
template<field_element_c_weak K, std::size_t P, std::size_t N>
constexpr auto roots_berlekamp(const polynomial_nttp<K, N>& f, K zero, K one) -> roots_result<K, N>
{
  roots_result<K, N> result;

  auto [factors, count] = berlekamp_factor<K, P, N>(f, zero, one);
  std::println("DEBUG: berlekamp_factor found {} factors", count);

  for (std::size_t i = 0; i < count; ++i)
  {
    auto& factor = factors[i];
    if (effective_degree(factor) == 1)
    {
      // Linear factor: a*x + b = 0 => x = -b/a
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
