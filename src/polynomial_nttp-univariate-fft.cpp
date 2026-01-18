/*
 *  polynomial_nttp-univariate-fft.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *    Implementation of Fast Fourier Transform (FFT) using:
 *      - Apple vDSP (Accelerate) if available (Zero-Copy Stride Trick)
 *      - Generic Fallback (Cooley-Tukey)
 */

#include <bit>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <span>
#include <vector>

#if defined(__APPLE__) && defined(LAM_USE_ACCELERATE)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Welaborated-enum-base"
#include <Accelerate/Accelerate.h>
#pragma clang diagnostic pop
#endif

// Forward declarations to match extern "C++" in module
void fft_accelerate_impl(std::span<std::complex<double>> data, bool inverse);
void fft_software_impl(std::span<std::complex<double>> data, bool inverse);

#if defined(__APPLE__) && defined(LAM_USE_ACCELERATE)
static std::mutex setup_cache_mutex;
static std::map<vDSP_Length, FFTSetupD> setup_cache;

FFTSetupD get_fft_setup(vDSP_Length log2n)
{
  std::lock_guard<std::mutex> lock(setup_cache_mutex);
  auto it = setup_cache.find(log2n);
  if (it != setup_cache.end())
  {
    return it->second;
  }
  FFTSetupD setup = vDSP_create_fftsetupD(log2n, FFT_RADIX2);
  setup_cache[log2n] = setup;
  return setup;
}
#endif

void fft_accelerate_impl(std::span<std::complex<double>> data, bool inverse)
{
#if defined(__APPLE__) && defined(LAM_USE_ACCELERATE)
  std::size_t n = data.size();
  if (n == 0)
    return;
  vDSP_Length log2n = static_cast<vDSP_Length>(std::bit_width(n) - 1);

  FFTSetupD setup = get_fft_setup(log2n);
  if (!setup)
  {
    fft_software_impl(data, inverse);
    return;
  }

  // Zero-Copy Stride Trick:
  // std::complex<double> standard layout [R, I, R, I...] matches
  // vDSP split complex with stride 2.
  auto* raw_ptr = reinterpret_cast<double*>(data.data());
  DSPDoubleSplitComplex split;
  split.realp = raw_ptr;
  split.imagp = raw_ptr + 1;

  constexpr vDSP_Stride stride = 2;

  // vDSP_fft_zipD is in-place complex-to-complex
  vDSP_fft_zipD(setup, &split, stride, log2n, inverse ? FFT_INVERSE : FFT_FORWARD);

  if (inverse)
  {
    // vDSP results are unscaled; apply 1/N normalization
    double scale = 1.0 / static_cast<double>(n);
    vDSP_vsmulD(raw_ptr, 1, &scale, raw_ptr, 1, n * 2);
  }

#else
  fft_software_impl(data, inverse);
#endif
}

// --------------------------------------------------------------------------
// Software Implementation (Generic C++)
// --------------------------------------------------------------------------

static void bit_reverse_permutation(std::span<std::complex<double>> a)
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
      std::swap(a[i], a[j]);
    }
  }
}

void fft_software_impl(std::span<std::complex<double>> data, bool inverse)
{
  std::size_t n = data.size();
  if (n <= 1)
    return;

  bit_reverse_permutation(data);

  // ACCURACY IMPROVEMENT: Use 'long double' for intermediate trig calculations
  long double angle_sign = inverse ? 1.0L : -1.0L;

  for (std::size_t len = 2; len <= n; len <<= 1)
  {
    long double angle = angle_sign * 2.0L * std::numbers::pi_v<long double> / static_cast<long double>(len);

    for (std::size_t i = 0; i < n; i += len)
    {
      for (std::size_t j = 0; j < len / 2; j++)
      {
        long double current_angle = angle * static_cast<long double>(j);
        // Compute in high precision, then cast to double
        std::complex<double> w = {static_cast<double>(std::cos(current_angle)),
                                  static_cast<double>(std::sin(current_angle))};

        std::complex<double> u = data[i + j];
        std::complex<double> v = data[i + j + len / 2] * w;

        data[i + j] = u + v;
        data[i + j + len / 2] = u - v;
      }
    }
  }

  if (inverse)
  {
    double inv_n = 1.0 / static_cast<double>(n);
    for (auto& x : data)
      x *= inv_n;
  }
}
