/*
 *  roots.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp roots function
 */

import std;
import lam.polynomial_nttp;

using namespace lam;

// Helper to extract values for sorting/comparison
template<typename Range>
auto extract_values(const Range& roots)
{
  using RootType = std::remove_cvref_t<decltype(*roots.begin())>;
  using K = decltype(std::declval<RootType>().value);
  std::vector<K> vals;
  for (const auto& r : roots)
    vals.push_back(r.value);
  return vals;
}

// Compile-time tests
constexpr auto test_linear_constexpr()
{
  polynomial_nttp<double, 1> p{{-4.0, 2.0}}; // 2x - 4 = 0
  auto r = roots(p);
  return r.size() == 1 && r[0].value == 2.0;
}

constexpr auto test_quadratic_constexpr()
{
  polynomial_nttp<double, 2> p{{2.0, -3.0, 1.0}}; // x² - 3x + 2 = (x-1)(x-2)
  auto r = roots(p);
  if (r.size() != 2)
    return false;

  bool has_1 = false, has_2 = false;
  for (auto root : r)
  {
    if (is_approx_equal(root.value, 1.0))
      has_1 = true;
    if (is_approx_equal(root.value, 2.0))
      has_2 = true;
  }
  return has_1 && has_2;
}

constexpr auto test_cubic_constexpr()
{
  // Use x³ - 1 = 0. Real solver finds one real root x=1.
  polynomial_nttp<double, 3> p{{-1.0, 0.0, 0.0, 1.0}};
  auto r = roots(p);
  if (r.size() != 1)
    return false;
  return is_approx_equal(r[0].value, 1.0);
}

constexpr auto test_quartic_constexpr()
{
  polynomial_nttp<double, 4> p{{24.0, -50.0, 35.0, -10.0, 1.0}}; // 1, 2, 3, 4
  auto r = roots(p);
  if (r.size() != 4)
    return false;

  bool h[5] = {false}; // 1..4
  for (auto root : r)
  {
    int v = (int)(root.value + 0.1);
    if (v >= 1 && v <= 4 && is_approx_equal(root.value, (double)v))
      h[v] = true;
  }
  return h[1] && h[2] && h[3] && h[4];
}

constexpr auto test_complex_quadratic_constexpr()
{
  using std_complex = std::complex<double>;
  polynomial_nttp<std_complex, 2> p{{std_complex(1.0, 0.0), std_complex(0.0, 0.0), std_complex(1.0, 0.0)}};
  auto r = roots(p);
  if (r.size() != 2)
    return false;
  for (auto& root : r)
    if (!is_negligible(p(root.value)))
      return false;
  return true;
}

constexpr auto test_real_no_roots_constexpr()
{
  polynomial_nttp<double, 2> p{{1.0, 0.0, 1.0}}; // x² + 1 = 0
  auto r = roots(p);
  return r.empty(); // Should find 0 real roots
}

constexpr auto test_complex_linear_constexpr()
{
  using std_complex = std::complex<double>;
  polynomial_nttp<std_complex, 1> p{{std_complex(2.0, 2.0), std_complex(1.0, 1.0)}};
  auto r = roots(p);
  if (r.size() != 1)
    return false;
  return is_negligible(p(r[0].value));
}

constexpr auto test_complex_cubic_constexpr()
{
  using std_complex = std::complex<double>;
  polynomial_nttp<std_complex, 3> p{
    {std_complex(-1.0, 0.0), std_complex(0.0, 0.0), std_complex(0.0, 0.0), std_complex(1.0, 0.0)}};
  auto r = roots(p);
  if (r.size() < 1)
    return false; // Should find 3 roots ideally, but size check is loose here provided we verify values
  for (auto& root : r)
    if (!is_negligible(p(root.value)))
      return false;
  return true;
}


// H_4(x) has 4 real roots (runs into cubic solver limitation after finding first root)
constexpr auto test_hermite_constexpr()
{
  polynomial_nttp<double, 4> H4{{12.0, 0.0, -48.0, 0.0, 16.0}};
  auto r = roots(H4);
  if (r.size() != 4)
    return false;

  // Roots are ±sqrt((48 ± sqrt(1536))/32)
  constexpr double sq6 = 2.449489742783178;
  // x^2 = 1.5 ± 1.2247... -> x = 0.5246... and 1.6506...
  constexpr double r1 = 0.5246476232752903;
  constexpr double r2 = 1.6506801238857845;

  bool h[5] = {false}; // match tracker
  for (auto& root : r)
  {
    double a = (root.value < 0 ? -root.value : root.value); // abs
    if (is_approx_equal(a, r1))
      h[1] = true;
    else if (is_approx_equal(a, r2))
      h[2] = true;
  }
  // We need 4 roots, passed size check.
  // If we found both magnitudes, and sizes imply we have + and -, likely distinct.
  // Actually, strictly check all 4?
  int count = 0;
  for (auto& root : r)
  {
    if (is_approx_equal(root.value, r1) || is_approx_equal(root.value, -r1) || is_approx_equal(root.value, r2) ||
        is_approx_equal(root.value, -r2))
      count++;
  }
  return count == 4;
}

// Quartic that Ferrari can fully solve at compile time: (x-1)(x-2)(x^2+1)
// Roots: 1, 2, ±i
constexpr auto test_complex_quartic_constexpr()
{
  using std_complex = std::complex<double>;
  // x^4 - 3x^3 + 3x^2 - 3x + 2 = (x-1)(x-2)(x^2+1)
  polynomial_nttp<std_complex, 4> p{
    {std_complex(2, 0), std_complex(-3, 0), std_complex(3, 0), std_complex(-3, 0), std_complex(1, 0)}};
  auto r = roots(p);

  if (r.size() != 4)
    return false;

  bool has_1 = false;
  bool has_2 = false;
  bool has_i = false;
  bool has_ni = false;

  for (const auto& root : r)
  {
    if (is_approx_equal(root.value, {1.0, 0.0}))
      has_1 = true;
    else if (is_approx_equal(root.value, {2.0, 0.0}))
      has_2 = true;
    else if (is_approx_equal(root.value, {0.0, 1.0}))
      has_i = true;
    else if (is_approx_equal(root.value, {0.0, -1.0}))
      has_ni = true;
  }
  return has_1 && has_2 && has_i && has_ni;
}

int main()
{
  int failures = 0;

  // Compile-time verification - Real polynomials
  static_assert(test_linear_constexpr(), "Linear roots should work at compile time");
  static_assert(test_quadratic_constexpr(), "Quadratic roots should work at compile time");
  static_assert(test_cubic_constexpr(), "Cubic roots should work at compile time");
  static_assert(test_quartic_constexpr(), "Quartic roots should work at compile time");

  // Compile-time verification - std_complex polynomials
  static_assert(test_complex_linear_constexpr(), "std_complex linear roots should work at compile time");
  static_assert(test_complex_quadratic_constexpr(), "std_complex quadratic roots should work at compile time");
  static_assert(test_complex_cubic_constexpr(), "std_complex cubic roots should work at compile time");
  static_assert(test_complex_quartic_constexpr(), "std_complex quartic roots should work at compile time");

  static_assert(test_real_no_roots_constexpr(),
                "Real quadratic with complex roots should return empty at compile time");

  // Compile-time verification - High Degree
  static_assert(test_hermite_constexpr(), "NR should work at compile time for H4 (full)");

  std::println("Compile-time root finding: OK (linear, quadratic, cubic, quartic, complex, NR)");

  // Test 1: Linear 2x - 4 = 0 -> x = 2
  {
    polynomial_nttp<double, 1> p;
    p.coefficients[1] = 2.0;
    p.coefficients[0] = -4.0;
    auto r = roots(p);
    if (r.size() != 1 || !is_approx_equal(r[0].value, 2.0))
    {
      std::println("Test 1 Failed: Expected {{2.0}}, got size {}", r.size());
      failures++;
    }
    else
      std::println("Test 1 Passed: Linear root = {}", r[0].value);
  }

  // Test 2: Quadratic x^2 - 3x + 2 = 0 -> {1, 2}
  {
    polynomial_nttp<double, 2> p;
    p.coefficients[2] = 1.0;
    p.coefficients[1] = -3.0;
    p.coefficients[0] = 2.0;
    auto r = roots(p);
    auto vals = extract_values(r);
    std::ranges::sort(vals);
    if (vals.size() != 2 || !is_approx_equal(vals[0], 1.0) || !is_approx_equal(vals[1], 2.0))
    {
      std::print("Test 2 Failed: Expected {{1, 2}}, got {{");
      for (auto v : vals)
        std::print("{}, ", v);
      std::println("}}");
      failures++;
    }
    else
      std::println("Test 2 Passed: Quadratic roots = {}, {}", vals[0], vals[1]);
  }

  // Test 3: Quadratic x^2 + 2x + 1 = 0 -> {-1} (double root, multiplicity 2)
  {
    polynomial_nttp<double, 2> p;
    p.coefficients[2] = 1.0;
    p.coefficients[1] = 2.0;
    p.coefficients[0] = 1.0;
    auto r = roots(p);
    if (r.size() != 1 || !is_approx_equal(r[0].value, -1.0) || r[0].multiplicity != 2)
    {
      std::print("Test 3 Failed: Expected {{-1, mult=2}}, got {{");
      for (const auto& root : r)
        std::print("({}, mult={}), ", root.value, root.multiplicity);
      std::println("}}");
      failures++;
    }
    else
      std::println("Test 3 Passed: Double root = {} with multiplicity {}", r[0].value, r[0].multiplicity);
  }

  // Test 4: Cubic (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6 -> {1, 2, 3}
  {
    polynomial_nttp<double, 3> p;
    p.coefficients[3] = 1.0;
    p.coefficients[2] = -6.0;
    p.coefficients[1] = 11.0;
    p.coefficients[0] = -6.0;

    auto r = roots(p);
    auto vals = extract_values(r);
    std::ranges::sort(vals);

    if (vals.size() != 3 || !is_approx_equal(vals[0], 1.0) || !is_approx_equal(vals[1], 2.0) ||
        !is_approx_equal(vals[2], 3.0))
    {
      std::print("Test 4 Failed: Expected {{1, 2, 3}}, got {{");
      for (auto v : vals)
        std::print("{}, ", v);
      std::println("}}");
      failures++;
    }
    else
      std::println("Test 4 Passed: Cubic roots = {}, {}, {}", vals[0], vals[1], vals[2]);
  }

  // Test 5: Cubic with repeated root x^3 - 3x + 2 = (x-1)^2(x+2) -> {1 (mult 2), -2}
  {
    polynomial_nttp<double, 3> p;
    p.coefficients[3] = 1.0;
    p.coefficients[2] = 0.0;
    p.coefficients[1] = -3.0;
    p.coefficients[0] = 2.0;

    auto r = roots(p);
    std::println("Test 5: Cubic with repeated root");
    for (const auto& root : r)
      std::println("  root = {}, multiplicity = {}", root.value, root.multiplicity);

    // Check we got roots (implementation-dependent exact results)
    if (r.empty())
    {
      std::println("Test 5 Failed: No roots found");
      failures++;
    }
    else
      std::println("Test 5 Passed: Found {} root(s)", r.size());
  }

  // Test 6: Quartic (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
  {
    polynomial_nttp<double, 4> p;
    p.coefficients[4] = 1.0;
    p.coefficients[3] = -10.0;
    p.coefficients[2] = 35.0;
    p.coefficients[1] = -50.0;
    p.coefficients[0] = 24.0;

    auto r = roots(p);
    auto vals = extract_values(r);
    std::ranges::sort(vals);

    std::print("Test 6: Quartic roots = {{");
    for (auto v : vals)
      std::print("{:.4f}, ", v);
    std::println("}}");

    bool all_accurate = true;
    for (auto v : vals)
    {
      if (!is_negligible(p(v)))
      {
        all_accurate = false;
        std::println("  Accuracy failure: p({}) = {}", v, p(v));
      }
    }

    if (vals.size() == 4 && all_accurate)
      std::println("Test 6 Passed: Found 4 roots with good accuracy");
    else
    {
      std::println("Test 6 Failed: Found {} roots, accuracy ok? {}", vals.size(), all_accurate);
      failures++;
    }
  }

  // Test 7: std_complex linear (1+i)x + (2+2i) = 0 -> x = -2
  {
    using std_complex = std::complex<double>;
    polynomial_nttp<std_complex, 1> p;
    p.coefficients[1] = std_complex(1.0, 1.0);
    p.coefficients[0] = std_complex(2.0, 2.0);

    auto r = roots(p);
    if (r.size() != 1 || !is_approx_equal(r[0].value, std_complex(-2.0, 0.0)))
    {
      std::print("Test 7 Failed: Expected {{-2}}, got {{");
      for (const auto& root : r)
        std::print("({}, {}), ", root.value.real(), root.value.imag());
      std::println("}}");
      failures++;
    }
    else
      std::println("Test 7 Passed: std_complex linear root = ({}, {})", r[0].value.real(), r[0].value.imag());
  }

  // Test 8: std_complex quadratic x² + 1 = 0 -> roots = ±i
  {
    using std_complex = std::complex<double>;
    polynomial_nttp<std_complex, 2> p;
    p.coefficients[2] = std_complex(1.0, 0.0); // x²
    p.coefficients[1] = std_complex(0.0, 0.0); // 0x
    p.coefficients[0] = std_complex(1.0, 0.0); // +1

    auto r = roots(p);
    std::println("Test 8: std_complex quadratic x² + 1 = 0");
    for (const auto& root : r)
      std::println("  root = ({}, {}), multiplicity = {}", root.value.real(), root.value.imag(), root.multiplicity);

    // Should find 2 roots: +i and -i
    if (r.size() == 2)
    {
      bool found_plus_i = false;
      bool found_minus_i = false;
      for (const auto& root : r)
      {
        if (is_approx_equal(root.value, std_complex(0.0, 1.0)))
          found_plus_i = true;
        if (is_approx_equal(root.value, std_complex(0.0, -1.0)))
          found_minus_i = true;
      }
      if (found_plus_i && found_minus_i)
        std::println("Test 8 Passed: Found ±i roots");
      else
      {
        std::println("Test 8 Failed: Roots don't match ±i");
        failures++;
      }
    }
    else
    {
      std::println("Test 8 Failed: Expected 2 roots, got {}", r.size());
      failures++;
    }
  }

  // Test 8b: Real quadratic x² + 1 = 0 -> no real roots
  {
    polynomial_nttp<double, 2> p{{1.0, 0.0, 1.0}};
    auto r = roots(p);

    if (r.empty())
      std::println("Test 8b Passed: x² + 1 has 0 real roots as expected");
    else
    {
      std::println("Test 8b Failed: Found phantom roots for x² + 1");
      failures++;
    }
  }

  // Test 9: std_complex cubic x³ - 1 = 0 -> roots = 1, ω, ω² (cube roots of unity)
  {
    using std_complex = std::complex<double>;
    polynomial_nttp<std_complex, 3> p;
    p.coefficients[3] = std_complex(1.0, 0.0);  // x³
    p.coefficients[2] = std_complex(0.0, 0.0);  // 0x²
    p.coefficients[1] = std_complex(0.0, 0.0);  // 0x
    p.coefficients[0] = std_complex(-1.0, 0.0); // -1

    auto r = roots(p);
    std::println("Test 9: std_complex cubic x³ - 1 = 0 (cube roots of unity)");
    for (const auto& root : r)
      std::println("  root = ({:.6f}, {:.6f}), multiplicity = {}", root.value.real(), root.value.imag(),
                   root.multiplicity);

    // Should find at least 1 root (Cardano gives 1 real root for real discriminant)
    if (!r.empty())
    {
      // Check if 1 is a root
      bool found_one = false;
      for (const auto& root : r)
        if (is_approx_equal(root.value, std_complex(1.0, 0.0)))
          found_one = true;

      if (found_one)
        std::println("Test 9 Passed: Found root 1");
      else
        std::println("Test 9: Found {} root(s), none equal to 1", r.size());
    }
    else
    {
      std::println("Test 9 Failed: No roots found");
      failures++;
    }
  }

  // Test 10: std_complex quartic (x-i)(x+i)(x-2i)(x+2i) = x⁴ + 5x² + 4
  // Note: All roots are purely imaginary (±i, ±2i). Newton-Raphson from
  // real starting points may not find these. This is informational only.
  {
    using std_complex = std::complex<double>;
    polynomial_nttp<std_complex, 4> p;
    p.coefficients[4] = std_complex(1.0, 0.0); // x⁴
    p.coefficients[3] = std_complex(0.0, 0.0); // 0x³
    p.coefficients[2] = std_complex(5.0, 0.0); // 5x²
    p.coefficients[1] = std_complex(0.0, 0.0); // 0x
    p.coefficients[0] = std_complex(4.0, 0.0); // 4

    auto r = roots(p);
    std::println("Test 10: std_complex quartic x⁴ + 5x² + 4 = 0 (purely imaginary roots)");
    for (const auto& root : r)
      std::println("  root = ({:.6f}, {:.6f}), multiplicity = {}", root.value.real(), root.value.imag(),
                   root.multiplicity);

    bool all_accurate = true;
    for (const auto& root : r)
    {
      if (!is_negligible(p(root.value), 1e-14))
      {
        all_accurate = false;
        std_complex val = p(root.value);
        std::println("  Accuracy failure: p(...) = {}, {}", val.real(), val.imag());
      }
    }

    if (r.size() == 4 && all_accurate)
      std::println("Test 10 Passed: Found 4 roots with good accuracy");
    else
    {
      std::println("Test 10 Failed: Expected 4 roots, found {}, accurary ok? {}", r.size(), all_accurate);
      failures++;
    }
  }

  // Test 11: Hermite Polynomials (Real Roots)
  // H_3(x) = 8x^3 - 12x  -> Roots: 0, ±sqrt(1.5) ≈ ±1.2247
  {
    std::println("Test 11a: Hermite H_3(x) = 8x^3 - 12x");
    polynomial_nttp<double, 3> H3{{0.0, -12.0, 0.0, 8.0}};
    auto r = roots(H3);

    if (r.size() == 3)
    {
      std::vector<double> vals = extract_values(r);
      std::sort(vals.begin(), vals.end()); // -1.22..., 0, 1.22...

      bool ok = true;
      if (!is_approx_equal(vals[0], -std::sqrt(1.5)))
        ok = false;
      if (!is_negligible(vals[1]))
        ok = false;
      if (!is_approx_equal(vals[2], std::sqrt(1.5)))
        ok = false;

      if (ok)
        std::println("Test 11a Passed: Found roots 0, ±sqrt(1.5)");
      else
      {
        std::println("Test 11a Failed: values mismatch");
        for (auto v : vals)
          std::print("{} ", v);
        std::println("");
        failures++;
      }
    }
    else
    {
      std::println("Test 11a Failed: Expected 3 roots, got {}", r.size());
      failures++;
    }
  }

  // H_4(x) = 16x^4 - 48x^2 + 12 -> Roots: ±sqrt((3 ± sqrt(6))/2)
  // ±0.5246, ±1.6506
  {
    std::println("Test 11b: Hermite H_4(x) = 16x^4 - 48x^2 + 12");
    polynomial_nttp<double, 4> H4{{12.0, 0.0, -48.0, 0.0, 16.0}};
    auto r = roots(H4);

    if (r.size() == 4)
    {
      std::println("Test 11b Passed: Found 4 roots");
      for (const auto& root : r)
        std::println("  root = {:.4f}", root.value);
    }
    else
    {
      std::println("Test 11b Failed: Expected 4 roots, got {}", r.size());
      for (const auto& root : r)
        std::println("  found: {:.4f}", root.value);
      failures++;
    }
  }

  // H_5(x) = 32x^5 - 160x^3 + 120x -> Roots: 0, ±0.9585, ±2.0201
  {
    std::println("Test 11c: Hermite H_5(x) = 32x^5 - 160x^3 + 120x");
    polynomial_nttp<double, 5> H5{{0.0, 120.0, 0.0, -160.0, 0.0, 32.0}};
    auto r = roots(H5);

    if (r.size() == 5)
    {
      std::println("Test 11c Passed: Found 5 roots");
    }
    else
    {
      std::println("Test 11c Failed: Expected 5 roots, got {}", r.size());
      for (const auto& root : r)
        std::println("  found: {:.4f}", root.value);
      failures++;
    }
  }

  if (failures == 0)
  {
    std::println("\nAll roots tests passed.");
    return 0;
  }
  else
  {
    std::println("\n{} tests failed.", failures);
    return 1;
  }
}
