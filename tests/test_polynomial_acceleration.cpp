/*
 *  test_polynomial_acceleration.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 */

import std;
import lam.concepts;
import lam.polynomial_nttp;

using namespace lam;

// Helper to fill a polynomial with data
template<typename R, std::size_t N>
void fill_poly(lam::polynomial_nttp<R, N>& p, R start_val)
{
  for (std::size_t i = 0; i <= N; ++i)
  {
    if constexpr (std::is_same_v<R, std::complex<double>> || std::is_same_v<R, std::complex<float>>)
    {
      p.coefficients[i] = {static_cast<typename R::value_type>(start_val + i), 0};
    }
    else
    {
      p.coefficients[i] = static_cast<R>(start_val + i);
    }
  }
}

template<typename R, std::size_t N>
void verify_poly(const lam::polynomial_nttp<R, N>& p, R start_val, std::string_view name)
{
  for (std::size_t i = 0; i <= N; ++i)
  {
    R expected;
    if constexpr (std::is_same_v<R, std::complex<double>> || std::is_same_v<R, std::complex<float>>)
    {
      expected = {static_cast<typename R::value_type>(start_val + i), 0};
    }
    else
    {
      expected = static_cast<R>(start_val + i);
    }

    if (!lam::is_approx_equal(p[i], expected))
    {
      if constexpr (std::is_same_v<R, std::complex<double>> || std::is_same_v<R, std::complex<float>>)
      {
        std::println("Mismatch in {} at index {}: got ({}, {}), expected ({}, {})", name, i, p[i].real(), p[i].imag(),
                     expected.real(), expected.imag());
      }
      else
      {
        std::println("Mismatch in {} at index {}: got {}, expected {}", name, i, p[i], expected);
      }
      std::exit(1);
    }
  }
  std::println("{} verified successfully.", name);
}

void test_large_addition()
{
  std::println("Testing Large Addition (N=100)...");
  constexpr std::size_t N = 100;
  lam::polynomial_nttp<double, N> p;
  lam::polynomial_nttp<double, N> q;

  fill_poly(p, 1.0); // p[i] = 1 + i
  fill_poly(q, 2.0); // q[i] = 2 + i

  // Expected: (p+q)[i] = 3 + 2*i
  // Calculation:
  // p[0]=1, q[0]=2 -> sum=3 = 3 + 0
  // p[1]=2, q[1]=3 -> sum=5 = 3 + 2
  // ...

  auto sum = p + q;

  for (std::size_t i = 0; i <= N; ++i)
  {
    double expected = 3.0 + 2.0 * i;
    if (!lam::is_approx_equal(sum[i], expected))
    {
      std::println("Mismatch in Addition at index {}: got {}, expected {}", i, sum[i], expected);
      std::exit(1);
    }
  }
  std::println("Large Addition Verified.");
}

void test_large_subtraction()
{
  std::println("Testing Large Subtraction (N=100)...");
  constexpr std::size_t N = 100;
  lam::polynomial_nttp<double, N> p;
  lam::polynomial_nttp<double, N> q;

  fill_poly(p, 10.0); // p[i] = 10 + i
  fill_poly(q, 1.0);  // q[i] = 1 + i

  // Expected: (p-q)[i] = 9

  auto diff = p - q;

  for (std::size_t i = 0; i <= N; ++i)
  {
    double expected = 9.0;
    if (!lam::is_approx_equal(diff[i], expected))
    {
      std::println("Mismatch in Subtraction at index {}: got {}, expected {}", i, diff[i], expected);
      std::exit(1);
    }
  }
  std::println("Large Subtraction Verified.");
}

void test_complex_addition()
{
  std::println("Testing Large Complex Addition (N=100)...");
  constexpr std::size_t N = 100;
  lam::polynomial_nttp<std::complex<double>, N> p;
  lam::polynomial_nttp<std::complex<double>, N> q;

  for (std::size_t i = 0; i <= N; ++i)
  {
    p.coefficients[i] = {static_cast<double>(i), 1.0};
    q.coefficients[i] = {static_cast<double>(2 * i), 2.0};
  }

  auto sum = p + q; // Expected: {3*i, 3.0}

  for (std::size_t i = 0; i <= N; ++i)
  {
    std::complex<double> expected = {static_cast<double>(3 * i), 3.0};
    if (!lam::is_approx_equal(sum[i], expected))
    {
      std::println("Mismatch in Complex Addition at index {}: got ({},{}), expected ({},{})", i, sum[i].real(),
                   sum[i].imag(), expected.real(), expected.imag());
      std::exit(1);
    }
  }
  std::println("Large Complex Addition Verified.");
}

void test_large_int_threading()
{
  std::println("Testing Large Int Addition (Threading Path, N=60000)...");
  // int is not accelerated, so this triggers TBB or std::jthread
  constexpr std::size_t N = 60000;
  lam::polynomial_nttp<int, N> p;
  lam::polynomial_nttp<int, N> q;

  // Fill
  for (std::size_t i = 0; i <= N; ++i)
  {
    p.coefficients[i] = 1;
    q.coefficients[i] = 2;
  }

  auto sum = p + q;

  for (std::size_t i = 0; i <= N; ++i)
  {
    if (sum[i] != 3)
    {
      std::println("Mismatch in Threaded Int Addition at index {}: got {}, expected 3", i, sum[i]);
      std::exit(1);
    }
  }
  std::println("Large Int Threading Verified.");

  auto diff = p - q;
  for (std::size_t i = 0; i <= N; ++i)
  {
    if (diff[i] != -1)
    {
      std::println("Mismatch in Threaded Int Subtraction at index {}: got {}, expected -1", i, diff[i]);
      std::exit(1);
    }
  }
  std::println("Large Int Subtraction Threading Verified.");
}

void test_complex_subtraction()
{
  std::println("Testing Large Complex Subtraction (N=100)...");
  constexpr std::size_t N = 100;
  lam::polynomial_nttp<std::complex<double>, N> p;
  lam::polynomial_nttp<std::complex<double>, N> q;

  // p = {2i, 2}, q = {i, 1} -> p-q = {i, 1}
  for (std::size_t i = 0; i <= N; ++i)
  {
    p.coefficients[i] = {static_cast<double>(2 * i), 2.0};
    q.coefficients[i] = {static_cast<double>(i), 1.0};
  }

  auto diff = p - q;

  for (std::size_t i = 0; i <= N; ++i)
  {
    std::complex<double> expected = {static_cast<double>(i), 1.0};
    if (!lam::is_approx_equal(diff[i], expected))
    {
      std::println("Mismatch in Complex Subtraction at index {}: got ({},{}), expected ({},{})", i, diff[i].real(),
                   diff[i].imag(), expected.real(), expected.imag());
      std::exit(1);
    }
  }
  std::println("Large Complex Subtraction Verified.");
}

void test_small_serial_fallback()
{
  std::println("Testing Small Serial Fallback (N=5)...");
  constexpr std::size_t N = 5;
  // int type, small size -> Should NOT trigger acceleration or threading.
  // Falls back to std::ranges::transform.
  lam::polynomial_nttp<int, N> p{1, 2, 3, 4, 5, 6};
  lam::polynomial_nttp<int, N> q{1, 1, 1, 1, 1, 1};

  // Sum: 2, 3, 4, 5, 6, 7
  auto sum = p + q;
  for (std::size_t i = 0; i <= N; ++i)
  {
    if (sum[i] != static_cast<int>(i + 2))
    {
      std::println("Mismatch in Serial Fallback at index {}: got {}, expected {}", i, sum[i], i + 2);
      std::exit(1);
    }
  }

  // Diff: 0, 1, 2, 3, 4, 5
  auto diff = p - q;
  for (std::size_t i = 0; i <= N; ++i)
  {
    if (diff[i] != static_cast<int>(i))
    {
      std::println("Mismatch in Serial Fallback Diff at index {}: got {}, expected {}", i, diff[i], i);
      std::exit(1);
    }
  }

  std::println("Small Serial Fallback Verified.");
}

int main()
{
  test_large_addition();
  test_large_subtraction();
  test_complex_addition();
  test_complex_subtraction();
  test_large_int_threading();
  test_small_serial_fallback();
  std::println("All acceleration tests passed!");
  return 0;
}
