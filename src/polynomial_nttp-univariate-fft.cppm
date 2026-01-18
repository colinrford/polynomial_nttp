/*
 *  polynomial_nttp-univariate-fft.cppm
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    FFT Module Interface:
 *    - Constexpr FFT implementation (Pure C++)
 *    - Runtime Dispatch (vDSP / Software Fallback)
 */

module;

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <numbers>
#include <span>
#include <vector>

export module lam.polynomial_nttp:univariate.fft;

// Import the Math module for constexpr sin/cos
import :univariate.math;

// Forward declare implementation functions (defined in .cpp)
extern "C++"
{
  void fft_accelerate_impl(std::span<std::complex<double>> data, bool inverse);
  void fft_software_impl(std::span<std::complex<double>> data, bool inverse);
}

export namespace lam::polynomial::univariate::fft
{

// Helper: Bit Reversal Permutation (Compile-Time Friendly)
constexpr void bit_reverse_permutation(std::span<std::complex<double>> a)
{
  std::size_t n = a.size();
  std::size_t j = 0;
  for (std::size_t i = 1; i < n; i++)
  {
    std::size_t bit = n >> 1;
    while (j & bit)
    {
      j ^= bit;
      bit >>= 1;
    }
    j ^= bit;
    if (i < j)
    {
      // std::swap is constexpr in C++20
      std::swap(a[i], a[j]);
    }
  }
}

// Pure C++ Constexpr FFT Implementation (Cooley-Tukey)
constexpr void fft_constexpr_impl(std::span<std::complex<double>> a, bool inverse)
{
  std::size_t n = a.size();
  if (n <= 1)
    return;

  bit_reverse_permutation(a);

  // Precompute direction
  // For standard FFT: exp(-2pi * i / N)
  // For Inverse FFT:  exp( 2pi * i / N)
  double angle_sign = inverse ? 1.0 : -1.0;

  for (std::size_t len = 2; len <= n; len <<= 1)
  {
    double angle = angle_sign * 2.0 * std::numbers::pi_v<double> / len;

    for (std::size_t i = 0; i < n; i += len)
    {
      for (std::size_t j = 0; j < len / 2; j++)
      {
        double current_angle = angle * j;
        std::complex<double> w = {math::cos(current_angle), math::sin(current_angle)};

        std::complex<double> u = a[i + j];
        std::complex<double> v = a[i + j + len / 2] * w;

        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
      }
    }
  }

  if (inverse)
  {
    double inv_n = 1.0 / static_cast<double>(n);
    for (auto& x : a)
      x *= inv_n;
  }
}

// Main In-Place Transform Interface (Span-based)
constexpr void fft_transform(std::span<std::complex<double>> data, bool inverse = false)
{
  if consteval
  {
    fft_constexpr_impl(data, inverse);
  }
  else
  {
    // Runtime Dispatch
    // If linked with Accelerate, this calls the vDSP version.
    // Otherwise it calls the generic software version.
    fft_accelerate_impl(data, inverse);
  }
}

// Vector Wrapper (Handles allocation and padding)
constexpr std::vector<std::complex<double>> fft(std::vector<std::complex<double>> data, bool inverse = false)
{
  // Pad to power of 2
  std::size_t n = data.size();
  if (std::popcount(n) != 1)
  {
    std::size_t next_pow2 = std::bit_ceil(n);
    data.resize(next_pow2, {0, 0});
  }

  fft_transform(data, inverse);

  return data; // Return by value (Move semantics)
}

// Convenience overload for real input
constexpr std::vector<std::complex<double>> fft(const std::vector<double>& real_data, bool inverse = false)
{
  std::vector<std::complex<double>> complex_data(real_data.size());
  for (size_t i = 0; i < real_data.size(); ++i)
  {
    complex_data[i] = {real_data[i], 0.0};
  }
  return fft(std::move(complex_data), inverse);
}

// Public exposure of generic software FFT for benchmarks (legacy wrapper)
void soft_fft(std::vector<std::complex<double>>& data, bool inverse) { fft_software_impl(data, inverse); }

} // namespace lam::polynomial::univariate::fft
