/*
 *  polynomial_nttp-univariate-interpolation.cppm – Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp:univariate.interpolation provides Lagrange interpolation
 */

export module lam.polynomial_nttp:univariate.interpolation;

import std;
import lam.concepts;
import :univariate.structure;
import :univariate.algebra;

namespace lam::polynomial::univariate
{

/**
 * @brief Lagrange interpolation
 *
 * Given a set of points (x_i, y_i), find the unique polynomial P(x) of degree at most N
 * such that P(x_i) = y_i for all i.
 *
 * @tparam K Field element type (satisfying field_element_c_weak)
 * @tparam N Degree of the resulting polynomial (number of points is N + 1)
 * @param x_vals Array of x coordinates
 * @param y_vals Array of y coordinates
 * @return The interpolating polynomial
 */
export template<lam::concepts::experimental::field_element_c_weak K, std::size_t N>
constexpr auto lagrange_interpolate(const std::array<K, N + 1>& x_vals, const std::array<K, N + 1>& y_vals)
  -> polynomial_nttp<K, N>
{
  // Result polynomial, initially zero
  polynomial_nttp<K, N> result{};

  // For each point j
  for (std::size_t j = 0; j <= N; ++j)
  {
    // Compute basis polynomial L_j(x) = product_{i != j} (x - x_i) / (x_j - x_i)

    // Denominator: product_{i != j} (x_j - x_i)
    K denominator = K(1);
    for (std::size_t i = 0; i <= N; ++i)
    {
      if (i == j)
        continue;
      denominator = denominator * (x_vals[j] - x_vals[i]);
    }

    // Numerator: product_{i != j} (x - x_i)
    // We start with 1 (degree 0 polynomial)
    polynomial_nttp<K, N> basis{K(1)};

    for (std::size_t i = 0; i <= N; ++i)
    {
      if (i == j)
        continue;

      // In-place multiplication logic for (x - c):
      // P(x) = p_k x^k + ... + p_0
      // P(x) * (x - c) = P(x)*x - c*P(x)

      K c = x_vals[i];
      for (std::size_t k = N; k > 0; --k)
      {
        basis.coefficients[k] = basis.coefficients[k - 1] - basis.coefficients[k] * c;
      }
      basis.coefficients[0] = -(basis.coefficients[0] * c); // binary minus
    }

    // Now basis is L_j(x) * denominator
    // Scale by y_vals[j] / denominator
    K scalar = y_vals[j] / denominator;

    // result += basis * scalar
    for (std::size_t k = 0; k <= N; ++k)
    {
      result.coefficients[k] = result.coefficients[k] + basis.coefficients[k] * scalar;
    }
  }

  return result;
}

} // namespace lam::polynomial::univariate

namespace lam
{
export using polynomial::univariate::lagrange_interpolate;
}
