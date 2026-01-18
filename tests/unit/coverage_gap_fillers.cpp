/*
 *  coverage_gap_fillers.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  Tests designed to hit specific coverage gaps (constexpr fallbacks, unused operators)
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Helper for approximate equality
bool is_close(double a, double b, double tol = 1e-9) { return std::abs(a - b) < tol; }

void terminate_on_failure(bool condition, const char* msg)
{
  if (!condition)
  {
    std::print("FAIL: {}\n", msg);
    std::exit(1);
  }
}

// A. Math Module Tests
void test_runtime_math()
{
  std::print("Testing Runtime Math...\n");

  // volatiles force runtime execution, bypassing consteval
  volatile double val_sin = 0.5;
  volatile double val_cos = 0.5;
  volatile double val_sqrt = 2.0;

  double res_sin = math::sin(val_sin);
  terminate_on_failure(is_close(res_sin, std::sin(val_sin)), "math::sin failed");

  double res_cos = math::cos(val_cos);
  terminate_on_failure(is_close(res_cos, std::cos(val_cos)), "math::cos failed");

  double res_sqrt = math::sqrt(val_sqrt);
  terminate_on_failure(is_close(res_sqrt, std::sqrt(val_sqrt)), "math::sqrt failed");
}

// B. Algebra Module Scalar Operators
void test_algebra_scalar_ops()
{
  std::print("Testing Algebra Scalar Ops...\n");

  polynomial_nttp<double, 2> p;
  p.coefficients = {1.0, 2.0, 3.0}; // 1 + 2x + 3x^2
  double scalar = 5.0;

  // 1. Scalar + Polynomial (r + p) -> coefficients[0] += r
  auto p1 = scalar + p;
  terminate_on_failure(is_close(p1[0], 6.0), "scalar + p failed");
  terminate_on_failure(is_close(p1[1], 2.0), "scalar + p failed (idx 1)");

  // 2. Polynomial + Scalar (p + r) -> coefficients[0] += r
  auto p2 = p + scalar;
  terminate_on_failure(is_close(p2[0], 6.0), "p + scalar failed");

  // 3. Scalar - Polynomial (r - p) -> r_minus_p = -p; r_minus_p[0] += r
  // -p is {-1, -2, -3}. +5 -> {4, -2, -3}
  auto p3 = scalar - p;
  terminate_on_failure(is_close(p3[0], 4.0), "scalar - p failed (idx 0)");
  terminate_on_failure(is_close(p3[1], -2.0), "scalar - p failed (idx 1)");

  // 4. Polynomial - Scalar (p - r) -> p_minus_r[0] -= r
  // {1, 2, 3} - 5 -> {-4, 2, 3}
  auto p4 = p - scalar;
  terminate_on_failure(is_close(p4[0], -4.0), "p - scalar failed (idx 0)");
  terminate_on_failure(is_close(p4[1], 2.0), "p - scalar failed (idx 1)");

  // 5. Scalar * Polynomial (r * p) -> all * r
  auto p5 = scalar * p;
  terminate_on_failure(is_close(p5[0], 5.0), "scalar * p failed (idx 0)");

  // 6. Polynomial * Scalar (p * r) -> all * r
  auto p6 = p * scalar;
  terminate_on_failure(is_close(p6[0], 5.0), "p * scalar failed (idx 0)");

  // 7. Polynomial / Scalar (p / r) -> all / r
  auto p7 = p / scalar;
  terminate_on_failure(is_close(p7[0], 0.2), "p / scalar failed (idx 0)");
}

// C. Structure & Roots Helpers
void test_roots_result_helpers()
{
  std::print("Testing Roots Result Helpers...\n");

  roots_result<double, 5> rr;
  rr.push(1.0);
  rr.push(2.0);

  // Test to_vector
  auto vec = rr.to_vector();
  terminate_on_failure(vec.size() == 2, "to_vector size wrong");
  terminate_on_failure(is_close(vec[0].value, 1.0), "to_vector val 0 wrong");

  // Test values()
  auto vals = rr.values();
  terminate_on_failure(is_close(vals[0], 1.0), "values() val 0 wrong");
  terminate_on_failure(is_close(vals[1], 2.0), "values() val 1 wrong");

  // Test append
  roots_result<double, 5> rr2;
  rr2.push(3.0);

  rr.append(rr2);
  terminate_on_failure(rr.count == 3, "append count wrong");
  terminate_on_failure(is_close(rr.data[2].value, 3.0), "append value wrong");
}
// D. Compatibility Module Tests
void test_compat_module()
{
  std::print("Testing Compatibility Module...\n");

  using namespace lam::polynomial::univariate::compat;

  // Check capabilities
  static_assert(type_capabilities<double>::is_approximate, "double should be approximate");
  static_assert(type_capabilities<int>::is_exact, "int should be exact");
  static_assert(type_capabilities<double>::can_sqrt, "double should have sqrt");

  // verify_approximate_structure
  verify_approximate_structure<double>(1.0, 2.0, 3.0);

  // verify_approximate_algebra
  verify_approximate_algebra<double>(1.0);

  // verify_approximate_field_algebra
  verify_approximate_field_algebra<double>(1.0, 2.0);

  // verify_approximate_poly_evaluation
  verify_approximate_poly_evaluation<double>(0.0, 1.0, 2.0, 3.0);

  // verify_approximate_poly_mult_coeffs
  verify_approximate_poly_mult_coeffs<double>(0.0, 1.0);

  // verify_approximate_algebraic_identity
  verify_approximate_algebraic_identity<double>(0.0, 1.0);

  // verify_approximate_higher_degree
  verify_approximate_higher_degree<double>(0.0, 1.0, 2.0, 3.0);

  // verify_calculus
  verify_calculus<double>(0.0, 1.0);

  // verify_approximate_linear_roots: 2x - 4 = 0 -> x = 2
  auto result = verify_approximate_linear_roots<double>(2.0, -4.0, 2.0);
  terminate_on_failure(result, "verify_approximate_linear_roots failed");
}

// E. Complex FFT Multiplication
void test_complex_fft_multiply()
{
  std::print("Testing Complex FFT Multiplication...\n");

  using C = std::complex<double>;

  // Degree 32 each -> sum = 64 >= threshold, triggers FFT path
  polynomial_nttp<C, 32> p1;
  polynomial_nttp<C, 32> p2;

  // Use volatile doubles to prevent consteval, then construct complex
  volatile double re = 1.0;
  volatile double im = 0.5;
  C val{re, im};
  std::fill(p1.begin(), p1.end(), val);
  std::fill(p2.begin(), p2.end(), val);

  auto p3 = p1 * p2;

  // Verify degree
  terminate_on_failure(p3.degree == 64, "complex FFT degree mismatch");

  // Verify non-zero result
  terminate_on_failure(std::abs(p3[0]) > 0.1, "complex FFT coeff[0] too small");
}

// F. Naive Multiplication with int (non-FFT type)
void test_naive_multiply_int()
{
  std::print("Testing Naive Multiplication (int)...\n");

  // int is not FFT-compatible, so this forces the naive O(n^2) path
  polynomial_nttp<int, 5> p1{{1, 2, 3, 4, 5, 6}};
  polynomial_nttp<int, 5> p2{{1, 1, 1, 1, 1, 1}};

  auto p3 = p1 * p2;

  // (1+2x+3x^2+4x^3+5x^4+6x^5) * (1+x+x^2+x^3+x^4+x^5)
  // Result degree = 10
  terminate_on_failure(p3.degree == 10, "naive int mult degree mismatch");

  // p3[0] = 1*1 = 1
  terminate_on_failure(p3[0] == 1, "naive int mult coeff[0] wrong");

  // p3[10] = 6*1 = 6
  terminate_on_failure(p3[10] == 6, "naive int mult coeff[10] wrong");
}
// G. Software FFT Fallback
void test_soft_fft()
{
  std::print("Testing Software FFT (soft_fft)...\n");

  // Create a power-of-2 sized vector
  std::vector<std::complex<double>> data(64, {1.0, 0.0});

  // Call soft_fft directly to test the software fallback path
  fft::soft_fft(data, false);

  // Verify it ran (first element should be sum = 64)
  terminate_on_failure(std::abs(data[0].real() - 64.0) < 1e-9, "soft_fft forward failed");

  // Now inverse
  fft::soft_fft(data, true);

  // After forward + inverse, should be back to original (within tolerance)
  terminate_on_failure(std::abs(data[0].real() - 1.0) < 1e-9, "soft_fft inverse failed");
}

// H. Real-Input FFT Overload
void test_fft_real_input()
{
  std::print("Testing FFT with Real Input...\n");

  // Create real input vector (non-power-of-2 to test padding)
  std::vector<double> real_data = {1.0, 2.0, 3.0, 4.0, 5.0};

  // This calls the fft(const std::vector<double>&, bool) overload
  auto result = fft::fft(real_data, false);

  // Should be padded to 8 (next power of 2)
  terminate_on_failure(result.size() == 8, "real-input fft size wrong");

  // Sum of inputs is 15, so DC component (idx 0) should be 15
  terminate_on_failure(std::abs(result[0].real() - 15.0) < 1e-9, "real-input fft DC wrong");
}

// I. Float Evaluation (Accelerate/BLAS checks)
void test_float_evaluation()
{
  std::print("Testing Float Evaluation...\n");

  constexpr int N = 100; // Degree >= 80 to trigger accelerated path
  polynomial_nttp<float, N> p;

  // p(x) = 1 + x + x^2 + ... + x^100
  std::fill(p.begin(), p.end(), 1.0f);

  // Evaluate at x = 1.0 -> sum of 101 coefficients = 101.0
  float val = p(1.0f);
  terminate_on_failure(std::abs(val - 101.0f) < 1e-4f, "float evaluation at x=1 failed");

  // Evaluate at x = 0.5
  // Sum of geometric series: (1 - x^101) / (1 - x)
  // (1 - 0.5^101) / 0.5 approx 1 / 0.5 = 2.0
  float val2 = p(0.5f);
  terminate_on_failure(std::abs(val2 - 2.0f) < 1e-4f, "float evaluation at x=0.5 failed");
}

// J. Complex Float Evaluation
void test_complex_float_evaluation()
{
  std::print("Testing Complex Float Evaluation...\n");

  constexpr int N = 100; // Degree >= 80
  using C = std::complex<float>;
  polynomial_nttp<C, N> p;

  // p(x) = (1+i) + (1+i)x + ...
  std::fill(p.begin(), p.end(), C{1.0f, 1.0f});

  // Evaluate at x = 1.0
  // Sum = 101 * (1+i)
  C val = p(C{1.0f, 0.0f});
  terminate_on_failure(std::abs(val.real() - 101.0f) < 1e-4f, "complex float real part failed");
  terminate_on_failure(std::abs(val.imag() - 101.0f) < 1e-4f, "complex float imag part failed");
}

int main()
{
  test_runtime_math();
  test_algebra_scalar_ops();
  test_roots_result_helpers();
  test_compat_module();
  test_complex_fft_multiply();
  test_naive_multiply_int();
  test_soft_fft();
  test_fft_real_input();
  test_float_evaluation();
  test_complex_float_evaluation();

  std::print("All coverage gap fillers passed.\n");
  return 0;
}
