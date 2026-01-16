/*
 *  parentheses_operator.cpp
 *  Unit test for polynomial_nttp parentheses operator() evaluation.
 *  Covers:
 *   - Compile-time evaluation (constexpr)
 *   - Runtime evaluation (Horner)
 *   - Runtime optimized paths (Accelerate/BLAS) for N >= 80
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

// -----------------------------------------------------------------------------
// Reference Implementation (Simple Horner's Method)
// -----------------------------------------------------------------------------
template<typename T, std::size_t N>
T ref_horner(const polynomial_nttp<T, N>& poly, T x)
{
  T result = poly.coefficients[N];
  for (std::size_t i = N; i > 0; --i)
  {
    result = result * x + poly.coefficients[i - 1];
  }
  return result;
}

// -----------------------------------------------------------------------------
// Test Helpers
// -----------------------------------------------------------------------------
void check(bool condition, const char* msg)
{
  if (!condition)
  {
    std::print("FAIL: {}\n", msg);
    std::exit(1);
  }
}

template<typename T>
void check_approx(T val, T expected, const char* msg)
{
  if (!is_approx_equal(val, expected, 1e-9))
  { // Loose tolerance for floating point accumulation
    // std::print can't easily print complex without formatting helpers, so we keep it simple
    if constexpr (std::is_same_v<T, double>)
      std::print("FAIL: {} (Got: {}, Expected: {})\n", msg, val, expected);
    else
      std::print("FAIL: {}\n", msg);
    std::exit(1);
  }
}

// -----------------------------------------------------------------------------
// 1. Compile-Time Tests (Static Assertions)
// -----------------------------------------------------------------------------
constexpr bool test_constexpr()
{
  // Monomial: x^2
  constexpr polynomial_nttp<int, 2> p2{0, 0, 1};
  static_assert(p2(2) == 4);
  static_assert(p2(-3) == 9);

  // Linear: 2x + 1
  constexpr polynomial_nttp<int, 1> p1{1, 2};
  static_assert(p1(0) == 1);
  static_assert(p1(10) == 21);

  // Constant: 42
  constexpr polynomial_nttp<int, 0> p0{42};
  static_assert(p0(999) == 42);

  return true;
}
static_assert(test_constexpr());

// -----------------------------------------------------------------------------
// 2. Runtime Tests
// -----------------------------------------------------------------------------
int main()
{
  std::print("Running parentheses_operator tests...\n");

  // A. Small Degree (Horner Path)
  {
    polynomial_nttp<double, 3> cubic{1.0, -2.0, 0.0, 1.0}; // x^3 - 2x + 1
    check_approx(cubic(0.0), 1.0, "Small degree at 0");
    check_approx(cubic(1.0), 0.0, "Small degree at 1");
    check_approx(cubic(2.0), 5.0, "Small degree at 2");
  }

  // B. Large Degree (N >= 80) -> Targets Accelerate / BLAS paths
  // We construct a polynomial where P(x) = \sum x^i = (x^{N+1}-1)/(x-1)
  // Actually, simpler to just use random coefficients and compare vs reference.
  {
    constexpr std::size_t N = 100;
    polynomial_nttp<double, N> poly_large;

    // Fill with alternating pattern
    for (std::size_t i = 0; i <= N; ++i)
    {
      poly_large.coefficients[i] = (i % 2 == 0) ? 1.0 : -0.5;
    }

    double x = 0.9; // |x| < 1 avoids explosion

    double val_optimized = poly_large(x);
    double val_ref = ref_horner(poly_large, x);

    // Check against reference implementation
    check_approx(val_optimized, val_ref, "Large Degree Double (Optimization Check)");

    // Explicitly check configuration status
    if constexpr (lam::polynomial::config::use_accelerate)
    {
      std::print("  -> Verified using Apple Accelerate path for N=100\n");
    }
    else if constexpr (lam::polynomial::config::use_blas)
    {
      std::print("  -> Verified using BLAS path for N=100\n");
    }
    else
    {
      std::print("  -> Verified using Generic Horner path for N=100\n");
    }
  }

  // C. Complex Large Degree
  {
    constexpr std::size_t N = 85;
    polynomial_nttp<std::complex<double>, N> poly_complex;

    for (std::size_t i = 0; i <= N; ++i)
    {
      poly_complex.coefficients[i] = {static_cast<double>(i), 1.0};
    }

    std::complex<double> x{0.5, 0.5};

    auto val_optimized = poly_complex(x);
    auto val_ref = ref_horner(poly_complex, x);

    check_approx(val_optimized, val_ref, "Large Degree Complex (Optimization Check)");

    if constexpr (lam::polynomial::config::use_accelerate)
    {
      std::print("  -> Verified using Apple Accelerate path (Complex) for N=85\n");
    }
    else if constexpr (lam::polynomial::config::use_blas)
    {
      std::print("  -> Verified using BLAS Complex path for N=85\n");
    }
    else
    {
      std::print("  -> Verified using Generic Complex Horner path for N=85\n");
    }
  }

  std::print("All tests passed.\n");
  return 0;
}
