/*
 *  division_complex.cpp - written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp polynomial division with structural complex-like coefficients
 *  (std::complex is not structural, so we mock it to test the generic logic)
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

// A structural type that mimics std::complex for testing purposes
template<typename T>
struct structural_complex
{
  T r; // real part
  T i; // imaginary part

  using value_type = T;

  constexpr T real() const { return r; }
  constexpr T imag() const { return i; }

  // Arithmetic logic needed for division
  constexpr structural_complex operator+(const structural_complex& other) const { return {r + other.r, i + other.i}; }
  constexpr structural_complex operator-(const structural_complex& other) const { return {r - other.r, i - other.i}; }
  constexpr structural_complex operator-() const { return {-r, -i}; }
  constexpr structural_complex operator*(const structural_complex& other) const
  { return {r * other.r - i * other.i, r * other.i + i * other.r}; }
  constexpr structural_complex operator/(const structural_complex& other) const
  {
    T d = other.r * other.r + other.i * other.i;
    return {(r * other.r + i * other.i) / d, (i * other.r - r * other.i) / d};
  }
  // needed for equality check in is_approx_equal (fallback) and non-is_approx_equal paths
  constexpr bool operator==(const structural_complex& other) const = default;

  // Construct from scalar (value_type) - implicit or explicit
  constexpr structural_complex(T v) : r(v), i(0) {}
  constexpr structural_complex(T re, T im) : r(re), i(im) {}
  constexpr structural_complex() : r(0), i(0) {}

  // Satisfy multiplicative_group_element_c_weak requirement for inversion/identity
  static constexpr structural_complex one() { return {static_cast<T>(1), static_cast<T>(0)}; }
};

// Helper helper for approx check of our mock
constexpr bool complex_approx_eq(structural_complex<double> a, structural_complex<double> b)
{
  constexpr double tol = 1e-7;
  auto diff = a - b;
  auto real_diff = diff.r;
  auto imag_diff = diff.i;
  if (real_diff < 0.0)
    real_diff = -real_diff;
  if (imag_diff < 0.0)
    imag_diff = -imag_diff;
  return (real_diff < tol) && (imag_diff < tol);
}

consteval bool test_complex_division()
{
  using Complex = structural_complex<double>;

  // P(x) = (1+i)x^2 + (2+2i)x + (3+3i)
  // Coefficients are stored [x^0, x^1, x^2...]
  constexpr polynomial_nttp<Complex, 2> p{
    Complex{3.0, 3.0}, // x^0
    Complex{2.0, 2.0}, // x^1
    Complex{1.0, 1.0}  // x^2
  };

  // D(x) = (1+i)x
  // Coefficients: 0*x^0 + (1+i)*x^1
  constexpr polynomial_nttp<Complex, 1> d{
    Complex{0.0, 0.0}, // x^0
    Complex{1.0, 1.0}  // x^1
  };

  // Expected Quotient: x + 2
  constexpr auto result = division_prototype<Complex, 2, p, 1, d>();
  constexpr auto q = result.first;
  constexpr auto r = result.second;

  // q = {2, 1} -> 2*x^0 + 1*x^1 = 2 + x
  bool q_ok = complex_approx_eq(q[0], Complex{2.0, 0.0}) && complex_approx_eq(q[1], Complex{1.0, 0.0});

  // Remainder should be ~ 3+3i (constant term)
  // The sizes might vary, but let's check the constant coeff
  bool r_ok = false;
  // Assuming standard storage: last coeff is constant term
  if (r.degree >= 0)
  {
    // Find constant term (last element in coefficients array if it were full,
    // but polynomial_nttp might be sparse or not?
    // structure.cppm: "index 0 is x^N"
    // so index N is x^0 constant.
    // For remainder 'r', degree is usually < deg(d)=1, so degree 0.
    // So r[0] is strictly x^deg, wait.

    // Let's rely on iterating or fixed indexing.
    // r should be degree 0 (or -1/empty if 0 polynomial, but here it's 3+3i).
    // If degree=0, r[0] is x^0
    r_ok = complex_approx_eq(r[r.degree], Complex{3.0, 3.0});
  }

  return q_ok && r_ok;
}

int main()
{
  static_assert(test_complex_division()); // Ensure compile-time validity
  return 0;
}
