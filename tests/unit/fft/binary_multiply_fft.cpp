/*
 *  binary_multiply_fft.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *    and for more info
 *
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

// Helper to check runtime conditions
void check(bool condition, const char* msg)
{
  if (!condition)
  {
    std::print("FAIL: {}\n", msg);
    std::exit(1);
  }
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

// -----------------------------------------------------------------------------
// Runtime Test (Verifies FFT Backend Dispatch)
// -----------------------------------------------------------------------------
void test_runtime_fft_dispatch()
{
  // Runtime values (prevent constexpr)
  volatile double val = 1.0;
  polynomial_nttp<double, 40> p1;
  polynomial_nttp<double, 40> p2;

  std::fill(p1.begin(), p1.end(), val);
  std::fill(p2.begin(), p2.end(), val);

  auto p3 = p1 * p2;

  check(p3.degree == 80, "Degree mismatch in runtime FFT");

  // Verify coefficients (Degree 80 triangular pattern)
  // p3[0] = 1, p3[1] = 2 ... p3[40] = 41 ...
  check(within_abs(p3[0], 1.0, tolerance), "p3[0] mismatch");
  check(within_abs(p3[40], 41.0, tolerance), "p3[40] mismatch");
  check(within_abs(p3[80], 1.0, tolerance), "p3[80] mismatch");

  std::print("  -> Verified Runtime FFT Dispatch (N=80)\n");
}

void test_fft_runtime_paths()
{
  std::print("Testing FFT Runtime Paths...\n");

  // 1. Exact power of 2 (32)
  std::vector<std::complex<double>> input_32(32, {1.0, 0.0});
  // Just testing it runs without crash and covers code
  fft::fft_transform(input_32, false);

  // 2. Non-power of 2 multiplication
  // p1 deg 20, p2 deg 20 -> result deg 40.
  // Next power of 2 is 64.
  // This forces "Runtime: Use FFT" path in algebra.cppm
  // which handles padding to 64.
  polynomial_nttp<double, 20> p1;
  polynomial_nttp<double, 20> p2;
  // volatile to ensure runtime
  volatile double v = 1.0;
  std::fill(p1.begin(), p1.end(), v);
  std::fill(p2.begin(), p2.end(), v);

  auto p3 = p1 * p2;
  // Check peak index 20 (coeff should be 21)
  check(within_abs(p3[20], 21.0, tolerance), "p3[20] mismatch (deg 20*20)");
}

int main()
{
  // static_assert(test_polynomial_multiply_fft(), "Constexpr multiplication failed");
  // static_assert(test_direct_fft(), "Constexpr direct FFT failed");

  std::print("Running binary_multiply_fft tests...\n");
  test_runtime_fft_dispatch();
  test_fft_runtime_paths();
  std::print("All tests passed.\n");

  return 0;
}
