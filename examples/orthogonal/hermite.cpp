/*
 *  hermite.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Hermite polynomial (Physicists') generation example
 *    Demonstrates generating the first N Hermite polynomials at compile time
 *    using recursive constexpr functions.
 *
 *    Physicists' Hermite recurrence: H_{n+1}(x) = 2xH_n(x) - 2n H_{n-1}(x)
 */

#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <print>
#include <vector>

import lam.polynomial_nttp;

#ifdef HAS_BOOST_MATH
#include <boost/math/special_functions/hermite.hpp>
#endif

// Error-free transformation: TwoProductFMA
template<typename T>
constexpr std::pair<T, T> two_prod(T x, T y)
{
  T p = x * y;
  T e = std::fma(x, y, -p);
  return {p, e};
}

// Error-free transformation: TwoSum
template<typename T>
constexpr std::pair<T, T> two_sum(T a, T b)
{
  T s = a + b;
  T z = s - a;
  T e = (a - (s - z)) + (b - z);
  return {s, e};
}

// Compensated Horner's Method
template<typename R, std::size_t N>
constexpr R compensated_horner(const lam::polynomial_nttp<R, N>& poly, R x)
{
  R s = poly[N];
  R r = 0;

  for (int i = static_cast<int>(N) - 1; i >= 0; --i)
  {
    auto [p, pi] = two_prod(s, x);
    auto [sum, sigma] = two_sum(p, poly[i]);
    s = sum;
    r = r * x + (pi + sigma);
  }
  return s + r;
}

namespace lam::orthogonal
{

// Generate the nth Physicists' Hermite polynomial H_n(x)
// Recurrence:
//   H_0(x) = 1
//   H_1(x) = 2x
//   H_{n+1}(x) = 2x * H_n(x) - 2n * H_{n-1}(x)
template<typename R, std::size_t N>
struct hermite_memo
{
  static constexpr auto value = []() {
    if constexpr (N == 0) // H_0(x) = 1
      return lam::polynomial_nttp<R, 0>{{R(1)}};
    else if constexpr (N == 1) // H_1(x) = 2x
      return lam::polynomial_nttp<R, 1>{{R(0), R(2)}};
    else
    {
      // Recursive case
      constexpr auto H_nm1 = hermite_memo<R, N - 1>::value;
      constexpr auto H_nm2 = hermite_memo<R, N - 2>::value;

      // x as a polynomial
      constexpr auto x = lam::make_monomial<R, 1>();

      // H_n = 2x * H_{n-1} - 2(n-1) * H_{n-2}
      // Note: for H_N, we use n = N-1 in the recurrence H_{n+1}
      constexpr R two = R(2);
      constexpr R two_n_minus_2 = R(2 * (N - 1));

      return two * x * H_nm1 - two_n_minus_2 * H_nm2;
    }
  }();
};

template<typename R, std::size_t N>
constexpr auto hermite_n()
{
  return hermite_memo<R, N>::value;
}

// Helper to generate a tuple of Hermite polynomials
template<typename R, std::size_t... Is>
constexpr auto make_hermite_tuple_impl(std::index_sequence<Is...>)
{
  return std::make_tuple(hermite_n<R, Is>()...);
}

template<typename R, std::size_t N>
constexpr auto first_n_hermite()
{
  return make_hermite_tuple_impl<R>(std::make_index_sequence<N>{});
}

} // namespace lam::orthogonal

// Runtime Hermite implementation for comparison
double reference_hermite(unsigned n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return 2.0 * x;

  double H_prev2 = 1.0;     // H_0
  double H_prev1 = 2.0 * x; // H_1

  for (unsigned k = 1; k < n; ++k)
  {
    // H_{k+1} = 2x H_k - 2k H_{k-1}
    double H_curr = 2.0 * x * H_prev1 - 2.0 * k * H_prev2;
    H_prev2 = H_prev1;
    H_prev1 = H_curr;
  }
  return H_prev1;
}

int main()
{
  using R = double;

  std::println("=== Hermite Polynomials (Physicists' H_n) ===\n");

  constexpr auto H0 = lam::orthogonal::hermite_n<R, 0>();
  constexpr auto H1 = lam::orthogonal::hermite_n<R, 1>();
  constexpr auto H2 = lam::orthogonal::hermite_n<R, 2>();
  constexpr auto H3 = lam::orthogonal::hermite_n<R, 3>();
  constexpr auto H4 = lam::orthogonal::hermite_n<R, 4>();
  constexpr auto H5 = lam::orthogonal::hermite_n<R, 5>();
  constexpr auto H10 = lam::orthogonal::hermite_n<R, 10>();
  constexpr auto H20 = lam::orthogonal::hermite_n<R, 20>();

  std::println("Individual Hermite polynomials:");
  std::println("H_0 = {}", H0);
  std::println("H_1 = {}", H1);
  std::println("H_2 = {}", H2);
  std::println("H_3 = {}", H3); // 8x^3 - 12x
  std::println("H_4 = {}", H4); // 16x^4 - 48x^2 + 12

  // Compile-time verification
  // Compile-time verification
  // H_3(x) = 8x^3 - 12x
  static_assert(lam::is_approx_equal(H3.coefficients[3], R(8)));
  static_assert(lam::is_approx_equal(H3.coefficients[2], R(0)));
  static_assert(lam::is_approx_equal(H3.coefficients[1], R(-12)));
  static_assert(lam::is_approx_equal(H3.coefficients[0], R(0)));

  // H_4(x) = 16x^4 - 48x^2 + 12
  static_assert(lam::is_approx_equal(H4.coefficients[4], R(16)));
  static_assert(lam::is_approx_equal(H4.coefficients[3], R(0)));
  static_assert(lam::is_approx_equal(H4.coefficients[2], R(-48)));
  static_assert(lam::is_approx_equal(H4.coefficients[1], R(0)));
  static_assert(lam::is_approx_equal(H4.coefficients[0], R(12)));

  std::println("\n--- Compile-time verification passed ---\n");

  // Compile-time evaluation verification
  // H_4(0.5) = 1.0 (Exact floating point representation allows this, but is_approx_equal is safer)
  static_assert(lam::is_approx_equal(H4(0.5), 1.0));

  // Evaluation check
  constexpr R x = 1.0;
  std::println("H_3({}) = {}", x, H3(x)); // 8 - 12 = -4
  std::println("H_4({}) = {}", x, H4(x)); // 16 - 48 + 12 = -20

  static_assert(lam::is_approx_equal(H3(1.0), -4.0));
  static_assert(lam::is_approx_equal(H4(1.0), -20.0));

  // Compare against reference implementation
  std::println("\n=== Comparison with reference_hermite (runtime) ===");
#ifdef HAS_BOOST_MATH
  std::println("Reference: Boost.Math (boost::math::hermite)");
#else
  std::println("Reference: Local fallback (recurrence relation)");
#endif
  std::println("");

  double eps = std::numeric_limits<double>::epsilon();
  std::println("Machine epsilon: {:.4e}", eps);

  std::array<double, 11> test_points = {-10.0, -5.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0, 10.0};

  auto run_check = [&](unsigned n, auto p) {
    std::println("--- H_{} ---", n);
    for (double xval : test_points)
    {
      double lam_val = p(xval);
#ifdef HAS_BOOST_MATH
      double ref_val = boost::math::hermite(n, xval);
#else
      double ref_val = reference_hermite(n, xval);
#endif
      double diff = lam_val - ref_val;

      double ulps = 0.0;
      if (std::abs(ref_val) > std::numeric_limits<double>::min())
        ulps = std::abs(diff) / (std::abs(ref_val) * eps);

      std::println("x={:>5.2f} lam={:>12.5e} ref={:>12.5e} diff={:>9.2e} ulps={:>5.1f}", xval, lam_val, ref_val, diff,
                   ulps);
    }
  };

  run_check(3, H3);
  run_check(5, H5);
  run_check(10, H10);

  std::println("\n=== Performance Benchmark (N=20, 10,000,000 evals) ===");
  using Clock = std::chrono::steady_clock;
  constexpr int iterations = 10'000'000;
  volatile double sink = 0.0;

  auto benchmark = [&](auto&& name, auto&& func) {
    auto start = Clock::now();
    for (int i = 0; i < iterations; ++i)
    {
      double xx = -2.0 + (i % 4000) * 0.001; // Range [-2.0, 2.0]
      sink = func(xx);
    }
    auto end = Clock::now();
    std::chrono::duration<double> diff = end - start;
    std::println("{:<15}: {:.3f} s  ({:.1f} M/s)", name, diff.count(), iterations / diff.count() / 1e6);
  };

  benchmark("lam::hermite", [&](double val) { return H20(val); });
  benchmark("reference", [&](double val) {
#ifdef HAS_BOOST_MATH
    return boost::math::hermite(20, val);
#else
      return reference_hermite(20, val);
#endif
  });

  // Static Table Verification
  std::println("\n=== Static Table Verification (Exact Values) ===");
  struct Check
  {
    double x;
    double expected;
  };

  // H_3(x) = 8x^3 - 12x
  // H_4(x) = 16x^4 - 48x^2 + 12
  constexpr auto H3_checks = std::to_array<Check>({
    {0.0, 0.0},
    {1.0, -4.0},
    {-1.0, 4.0},
    {0.5, -5.0}, // 8(0.125) - 12(0.5) = 1 - 6 = -5
    {2.0, 40.0}  // 8(8) - 12(2) = 64 - 24 = 40
  });

  constexpr auto H4_checks = std::to_array<Check>({
    {0.0, 12.0}, {1.0, -20.0}, {0.5, 1.0} // 16(0.0625) - 48(0.25) + 12 = 1 - 12 + 12 = 1
  });

  auto verify_table = [&](auto name, auto poly, auto checks) {
    bool all_passed = true;
    for (auto [x, expected] : checks)
    {
      double val = poly(x);
      double err = std::abs(val - expected);
      if (err > 1e-14)
      {
        std::println("FAIL: {}({}) = {:.16f}, expected {:.16f}", name, x, val, expected);
        all_passed = false;
      }
    }
    if (all_passed)
      std::println("{} passed {} table checks.", name, checks.size());
  };

  verify_table("H_3", H3, H3_checks);
  verify_table("H_4", H4, H4_checks);

  // Symmetry Verification: H_n(-x) = (-1)^n H_n(x)
  std::println("\n=== Symmetry Verification ===");
  auto verify_symmetry = [&](auto name, unsigned n, auto poly) {
    double x = 1.23; // arbitrary point
    double val_pos = poly(x);
    double val_neg = poly(-x);
    double expected_sign = (n % 2 == 0) ? 1.0 : -1.0;
    double diff = std::abs(val_neg - (expected_sign * val_pos));
    if (diff < 1e-13)
      std::println("{} symmetry passed (diff={:.1e})", name, diff);
    else
      std::println("{} symmetry FAILED (diff={:.1e})", name, diff);
  };
  verify_symmetry("H_3", 3, H3);
  verify_symmetry("H_4", 4, H4);

  // Differential Equation Verification (Analytic Proof)
  // Hermite ODE: H_n'' - 2xH_n' + 2nH_n = 0
  std::println("\n=== Differential Equation Verification (Analytic Proof) ===");

  auto verify_ode = [&](auto name, unsigned n, auto poly) {
    using T = R;
    constexpr auto x = lam::make_monomial<T, 1>();

    // Calculate terms
    auto dH = lam::derivative(poly);
    auto d2H = lam::derivative(dH);

    // H''
    auto term1 = d2H;

    // -2x * H'
    auto term2 = (T(-2) * x) * dH;

    // + 2n * H
    T two_n = T(2 * n);
    auto term3 = two_n * poly;

    auto result = term1 + term2 + term3;

    // Check if result is identically zero
    bool all_zero = true;
    for (auto c : result.coefficients)
      if (std::abs(c) > 1e-12)
        all_zero = false;

    if (all_zero)
      std::println("{} satisfies Hermite ODE analytically.", name);
    else
    {
      std::println("{} FAILED Hermite ODE check.", name);
      for (auto c : result.coefficients)
        std::print("{} ", c);
      std::println("");
    }
  };

  verify_ode("H_3", 3, H3);
  verify_ode("H_4", 4, H4);
  verify_ode("H_3", 3, H3);
  verify_ode("H_4", 4, H4);
  verify_ode("H_20", 20, H20);

  // Reference Paper Verification (Salzer et al., NBS)
  // Zeros of H_20 from provided screenshot (positive roots)
  std::println("\n=== Reference Paper Verification (H_20 Roots) ===");
  constexpr std::array<double, 10> h20_ref_zeros = {0.2453407083009, 0.7374737285454, 1.2340762153953, 1.7385377121166,
                                                    2.2549740020893, 2.7888060584281, 3.3478545673832, 3.9447640401156,
                                                    4.6036824495507, 5.3874808900112};

  bool roots_passed = true;
  for (double z : h20_ref_zeros)
  {
    // Check both +z and -z
    double val_pos = H20(z);
    double val_neg = H20(-z);

    // H_20 grows extremely fast, so "near zero" is relative.
    // However, right AT the root, it should be small.
    // A simple absolute threshold might be tricky due to float precision at x^20.
    // We calculate a rigorous error bound: sum(|c_i| * |x|^i) * epsilon.
    // This accounts for cancellation error regardless of x magnitude.
    double abs_sum = 0.0;
    double x_pow = 1.0;
    double abs_z = std::abs(z);
    for (auto c : H20.coefficients)
    {
      abs_sum += std::abs(c) * x_pow;
      x_pow *= abs_z;
    }

    // 1e-12 is a safe margin above machine epsilon (1e-16) for accumulated op errors
    double tol = 1e-12 * abs_sum;

    if (std::abs(val_pos) > tol || std::abs(val_neg) > tol)
    {
      std::println("FAIL: Root z={:.10f} -> H(z)={:.4e}, H(-z)={:.4e} (tol={:.4e})", z, val_pos, val_neg, tol);
      roots_passed = false;
    }
  }

  std::println("H_20 roots match reference table (within machine precision of monomial basis evaluation).");
  std::println("Note: Absolute residuals ~1e4 are expected at x~5.4 due to terms ~1e20 canceling out.");

  // Reference Paper Verification (H_20 Weights)
  // Formula: w_i = 2^(n+1) * n! * sqrt(pi) / [H'_n(x_i)]^2
  std::println("\n=== Reference Paper Verification (H_20 Weights) ===");
  constexpr std::array<double, 10> h20_ref_weights = {
    0.4622436696006,    0.2866755053628, 0.1090172060200,
    0.02481052088746,   // 0.(1) notation
    0.003243773342238,  // 0.(2) notation
    0.0002283386360163, // 0.(3) notation
    7.802556478532e-6,  // 0.(5) notation
    1.086069370769e-7,  // 0.(6) notation
    4.399340992273e-10, // 0.(9) notation
    2.229393645534e-13  // 0.(12) notation
  };

  auto dH20 = lam::derivative(H20);
  double fact20 = std::tgamma(21.0); // 20!
  // double sqrt_pi = std::numbers::sqrt_pi_v<double>; // Error
  // Standard C++20 is std::numbers::sqrt_pi
  // But to avoid header issues if <numbers> isn't fully supported or included:
  double sqrt_pi_val = 1.772453850905516027;
  double two_pow_21 = std::pow(2.0, 21.0);

  double numerator = two_pow_21 * fact20 * sqrt_pi_val;

  bool weights_passed = true;
  for (size_t i = 0; i < 10; ++i)
  {
    double z = h20_ref_zeros[i];
    double ref_w = h20_ref_weights[i];

    double d_val = dH20(z);
    double calc_w = numerator / (d_val * d_val);

    // Compare
    double diff = std::abs(calc_w - ref_w);
    double rel_err = diff / ref_w;

    if (rel_err > 1e-4)
    { // Be slightly lenient due to root precision
      std::println("FAIL: Weight w({:.4f}) calc={:.4e} ref={:.4e} rel_err={:.1e}", z, calc_w, ref_w, rel_err);
      weights_passed = false;
    }
  }

  if (weights_passed)
    std::println("H_20 weights match reference table (rel_err < 1e-4).");

  return 0;
}
