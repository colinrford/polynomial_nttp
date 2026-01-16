/*
 *  binary_multiply_fft.cpp
 *  Unit tests for FFT-based polynomial multiplication using polynomial_nttp.
 */

import std;
import lam.polynomial_nttp;

using namespace lam;
using namespace lam::polynomial::univariate;

constexpr double tolerance = 1e-10;

constexpr bool within_abs(double a, double b, double tol)
{
  double diff = a - b;
  return diff > -tol && diff < tol;
}

// Helper to create a polynomial with all coefficients = 1.0
// Degree 35 polynomial: 1 + x + x^2 + ... + x^35
template<std::size_t N>
consteval auto make_ones_polynomial()
{
  std::array<double, N + 1> coeffs{};
  for (std::size_t i = 0; i <= N; ++i)
    coeffs[i] = 1.0;
  return polynomial_nttp<double, N>(coeffs);
}

// Test polynomial multiplication triggering FFT path
// Degree 35 + Degree 35 = Degree 70 (> 64 threshold)
consteval bool test_polynomial_multiply_fft()
{
  constexpr auto p1 = make_ones_polynomial<35>();
  constexpr auto p2 = make_ones_polynomial<35>();

  // Multiply: (1 + x + ... + x^35)^2
  // Result coefficients form a triangular pattern: 1, 2, 3, ..., 36, 35, ..., 1
  auto p3 = p1 * p2;

  // Check degree
  if (p3.degree != 70)
    return false;

  // Check key coefficients
  // p3[0] = 1 (only 1*1)
  // p3[1] = 2 (1*x + x*1)
  // p3[35] = 36 (peak)
  // p3[70] = 1 (only x^35 * x^35)
  if (!within_abs(p3[0], 1.0, tolerance))
    return false;
  if (!within_abs(p3[1], 2.0, tolerance))
    return false;
  if (!within_abs(p3[35], 36.0, tolerance))
    return false;
  if (!within_abs(p3[70], 1.0, tolerance))
    return false;

  return true;
}

// Test direct FFT implementation: (1+x)(1-x) = 1-x^2
consteval bool test_direct_fft()
{
  // Signals: [1, 1], [1, -1] -> Pad to 4
  std::vector<std::complex<double>> signal1 = {{1, 0}, {1, 0}, {0, 0}, {0, 0}};
  std::vector<std::complex<double>> signal2 = {{1, 0}, {-1, 0}, {0, 0}, {0, 0}};

  // Transform
  auto f1 = fft::fft(signal1, false);
  auto f2 = fft::fft(signal2, false);

  // Multiply Pointwise
  std::vector<std::complex<double>> f3(4);
  for (int i = 0; i < 4; ++i)
    f3[i] = f1[i] * f2[i];

  // Inverse
  auto result = fft::fft(f3, true);

  // Expected: 1, 0, -1, 0
  if (!within_abs(result[0].real(), 1.0, tolerance))
    return false;
  if (!within_abs(result[1].real(), 0.0, tolerance))
    return false;
  if (!within_abs(result[2].real(), -1.0, tolerance))
    return false;
  if (!within_abs(result[3].real(), 0.0, tolerance))
    return false;

  return true;
}

int main()
{
  constexpr bool multiply_ok = test_polynomial_multiply_fft();
  constexpr bool direct_fft_ok = test_direct_fft();

  if constexpr (multiply_ok && direct_fft_ok)
  {
    return 0; // pass
  }
  else
  {
    return 1; // fail
  }
}
