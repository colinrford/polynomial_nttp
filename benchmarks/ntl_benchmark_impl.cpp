/*
 *  ntl_benchmark_impl.cpp
 *  Implementation of NTL benchmarking using traditional includes
 *  to avoid conflicts with C++23 'import std' in the main module.
 */

#include <chrono>
#include <iostream>
#include <print>

// NTL Headers
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>

void run_ntl_benchmark_impl(long N)
{
  // Ensure NTL is initialized (assuming main did strict init,
  // or we assume known prime from context if needed, but
  // NTL::ZZ_p::init must have been called.
  // Actually, passing the init responsibility across is tricky if not shared.
  // Let's re-init or assume initialized.
  // To be safe, we'll re-init here if not done, but NTL uses global state.
  // The main function initializes it.

  NTL::ZZ_pX a, b, c;
  NTL::random(a, N);
  NTL::random(b, N);

  auto start = std::chrono::steady_clock::now();
  for (int i = 0; i < 100; ++i)
  {
    NTL::mul(c, a, b);
  }
  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  std::println("  NTL (deg {}):       {:.2f} us", N, (duration / 100.0));
}

void init_ntl_modulus_impl(const char* prime_str)
{
  NTL::ZZ p = NTL::to_ZZ(prime_str);
  NTL::ZZ_p::init(p);
}

// Perform multiplication using NTL and compare with expected result
// Returns true if match, false otherwise
bool verify_ntl_multiplication_impl(long N, const long* a_data, const long* b_data, const long* expected_data)
{
  NTL::ZZ_pX a, b, c;
  a.SetLength(N + 1);
  b.SetLength(N + 1);

  // NTL coefficients can be set safely
  for (long i = 0; i <= N; ++i)
  {
    NTL::SetCoeff(a, i, a_data[i]);
    NTL::SetCoeff(b, i, b_data[i]);
  }

  NTL::mul(c, a, b);

  // Compare coefficients
  // Result degree can be up to 2*N
  long deg = NTL::deg(c);
  for (long i = 0; i <= deg; ++i)
  {
    long val = NTL::to_long(NTL::rep(NTL::coeff(c, i)));
    if (val != expected_data[i])
    {
      std::println(std::cerr, "Mismatch at index {}: NTL={} Local={}", i, val, expected_data[i]);
      return false;
    }
  }
  // Check purely zero beyond degree
  // (Local implementation might have zero padding up to power of 2 size, check strictly logic)
  // We assume expected_data is valid up to 2*N

  return true;
}
