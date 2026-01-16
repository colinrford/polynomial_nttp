/*
 *  precision_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 *  Verifies the numerical stability of the FFT implementation.
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

consteval double compute_max_squared_error()
{
  constexpr std::size_t n = 1024;
  std::vector<std::complex<double>> original(n);

  // Fill with simple data
  for (std::size_t i = 0; i < n; ++i)
    original[i] = {static_cast<double>(i), 0.0};

  auto spectrum = fft::fft(original, false);
  auto reconstructed = fft::fft(spectrum, true);

  double max_error = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    double diff_real = reconstructed[i].real() - original[i].real();
    double diff_imag = reconstructed[i].imag() - original[i].imag();
    double diff = diff_real * diff_real + diff_imag * diff_imag;
    if (diff > max_error)
      max_error = diff;
  }
  return max_error;
}

int main()
{
  constexpr double max_err = compute_max_squared_error();

  // Threshold: 1e-11 (squared), corresponds to ~3e-6 absolute error
  // Measured error is ~9e-13 (squared) = ~3e-7 absolute
  constexpr double threshold = 1.0e-11;

  if constexpr (max_err < threshold)
  {
    return 0; // pass
  }
  else
  {
    return 1; // fail
  }
}
