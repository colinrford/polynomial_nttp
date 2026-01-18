/*
 *  binary_multiply_polynomials.cpp - written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Unit test for polynomial * polynomial multiplication.
 *    Checks both naive and FFT paths using return codes.
 */
import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

bool test_naive()
{
  polynomial_nttp<double, 2> p{{1.0, 2.0, 3.0}};
  polynomial_nttp<double, 2> q{{4.0, 5.0, 6.0}};
  // (1 + 2x + 3x^2)(4 + 5x + 6x^2)
  // = 4 + 13x + 28x^2 + 27x^3 + 18x^4
  auto result = p * q;

  auto is_close = [](double a, double b) { return std::abs(a - b) < 1e-9; };

  if (!is_close(result[0], 4.0))
    return false;
  if (!is_close(result[1], 13.0))
    return false;
  if (!is_close(result[2], 28.0))
    return false;
  if (!is_close(result[3], 27.0))
    return false;
  if (!is_close(result[4], 18.0))
    return false;
  return true;
}

bool test_fft_trigger()
{
  // Trigger FFT (N >= 32)
  // Check correctness of split method (complex) and regular (real)
  constexpr std::size_t N = 32;
  polynomial_nttp<double, N> p{};
  p.coefficients[N] = 1.0;
  polynomial_nttp<double, N> q{};
  q.coefficients[N] = 1.0;

  auto result = p * q;

  // result should be x^64, so coeff[64] = 1, others 0
  if (std::abs(result[2 * N] - 1.0) > 1e-9)
    return false;
  if (std::abs(result[0]) > 1e-9)
    return false;

  // Complex Case (vDSP Split Trick)
  polynomial_nttp<std::complex<double>, N> pc{};
  pc.coefficients[N] = {1.0, 0.0};
  polynomial_nttp<std::complex<double>, N> qc{};
  qc.coefficients[N] = {0.0, 1.0}; // i*x^32

  auto result_c = pc * qc; // i*x^64

  auto diff = result_c[2 * N] - std::complex<double>(0.0, 1.0);
  if (std::abs(diff) > 1e-9)
    return false;

  return true;
}

// -----------------------------------------------------------------------------
// Compile-Time Tests
// -----------------------------------------------------------------------------
constexpr bool test_constexpr_naive()
{
  polynomial_nttp<int, 2> p{{1, 2, 3}};
  polynomial_nttp<int, 2> q{{4, 5, 6}};
  auto res = p * q;
  // 1*4 = 4
  if (res[0] != 4)
    return false;
  // 1*5 + 2*4 = 13
  if (res[1] != 13)
    return false;
  // 1*6 + 2*5 + 3*4 = 6 + 10 + 12 = 28
  if (res[2] != 28)
    return false;
  return true;
}
static_assert(test_constexpr_naive());

constexpr bool test_constexpr_fft()
{
  // N=16, M=16 -> M+N=32 (FFT Threshold)
  constexpr std::size_t N = 16;
  polynomial_nttp<double, N> p{};
  p.coefficients[N] = 1.0; // x^16
  polynomial_nttp<double, N> q{};
  q.coefficients[N] = 2.0; // 2x^16

  // Result = 2x^32
  auto res = p * q; // Should trigger compile-time FFT path

  if (!is_approx_equal(res[2 * N], 2.0, 1e-9))
    return false;
  if (!is_approx_equal(res[0], 0.0, 1e-9))
    return false;
  return true;
}
// Note: This requires -fconstexpr-steps bump if N is huge, but N=16 should be fine for simple FFT
static_assert(test_constexpr_fft());

// -----------------------------------------------------------------------------
// Fallback Tests (Non-Optimized Types)
// -----------------------------------------------------------------------------
bool test_large_degree_fallback()
{
  // int is not optimized, so even if N >= 32, it should use Naive
  constexpr std::size_t N = 40;
  polynomial_nttp<int, N> p;
  p.coefficients[0] = 1;
  p.coefficients[N] = 1; // 1 + x^40

  auto res = p * p; // (1 + x^40)^2 = 1 + 2x^40 + x^80

  if (res[0] != 1)
    return false;
  if (res[N] != 2)
    return false;
  if (res[2 * N] != 1)
    return false;
  if (res[N - 1] != 0)
    return false;

  return true;
}

int main()
{
  if (!test_naive())
    return 1;
  if (!test_fft_trigger())
    return 1;
  if (!test_large_degree_fallback())
    return 1;

  return 0; // All passed
}
