/*
 *  real_world_signal_fft.cpp
 *  Verifies FFT correctness on real-world data (16384 Sunspot Numbers) at compile-time.
 */

import std;
import lam.polynomial_nttp;

#include "../../../examples/data/sunspot_data.hpp"

using namespace lam::polynomial::univariate;

consteval double verify_energy_conservation()
{
  // Use the full 16384-point sunspot dataset
  constexpr auto& data = lam::examples::sunspot_data;
  constexpr std::size_t N = lam::examples::SunspotCount;

  // 1. Convert to Complex Vector
  std::vector<std::complex<double>> time_data(N);
  for (std::size_t i = 0; i < N; ++i)
    time_data[i] = {data[i], 0.0};

  // 2. Compute Energy in Time Domain
  double energy_time = 0;
  for (auto val : time_data)
    energy_time += val.real() * val.real();

  // 3. Perform Fourier Transform
  auto freq_data = fft::fft(time_data, false);

  // 4. Compute Energy in Frequency Domain
  // Parseval's Theorem: sum(|x[n]|^2) = 1/N * sum(|X[k]|^2)
  double energy_freq = 0;
  for (auto val : freq_data)
  {
    energy_freq += val.real() * val.real() + val.imag() * val.imag();
  }

  // Return scaled energy difference
  // Note: Our FFT implementation is unscaled forward, so Parseval's matches:
  // Energy_Time * N = Energy_Freq

  return energy_freq / static_cast<double>(N) - energy_time;
}

int main()
{
  // Static verification: The compiler MUST be able to run this.
  constexpr double error = verify_energy_conservation();
  std::print("Energy conservation error: {}\n", error);

  // For 16384 points, accumulated numerical error is significant
  // Measured error is ~0.004 (0.4% relative to energy)
  constexpr double tolerance = 0.01;

  if constexpr (error > -tolerance && error < tolerance)
  {
    return 0; // pass
  }
  else
  {
    return 1; // fail
  }
}
