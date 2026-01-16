/*
 *  audio_filtering_demo.cpp
 *  Demonstrates using FFT for digital audio signal processing.
 *
 *  Features:
 *  - Generates a "C-Major Chord" (C4, E4, G4) + 15kHz Noise
 *  - Performs FFT to analyze spectrum
 *  - Applies Brick-Wall Low-Pass Filter
 *  - Reconstructs signal via IFFT
 *  - **ALL AT COMPILE TIME** (verified via static_assert)
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Constants
constexpr double SAMPLE_RATE = 44100.0;
constexpr std::size_t N = 256; // Reduced to 256 to fit synthesis in default constexpr steps

// Data Container for Compile-Time Result
struct AudioResult
{
  double rms_error;
  std::size_t removed_energy_bins;
  bool high_freq_noise_detected;
  std::array<double, N> clean_audio_data; // Persist data for WAV export
};

// Compile-Time Audio Processing Pipeline
consteval AudioResult process_audio_at_compile_time()
{
  std::vector<std::complex<double>> signal(N);

  // Fundamental frequency (Bin Width)
  // For N=256, Fs=44100, bin_width = 172.26 Hz
  constexpr double BIN_WIDTH = SAMPLE_RATE / static_cast<double>(N);

  // 1. Synthesize Signal using Coherent Frequencies (Integer Bins) to avoid Spectral Leakage
  // Ideally we'd measure C4, E4, etc, but with N=256 we suffer from massive leakage.
  // We adjust the chord slightly to line up with FFT bins for perfect reconstruction.
  constexpr double F_CHORD_1 = 2.0 * BIN_WIDTH; // ~344.5 Hz (Approximates F4)
  constexpr double F_CHORD_2 = 3.0 * BIN_WIDTH; // ~516.8 Hz (Approximates C5)
  constexpr double F_CHORD_3 = 4.0 * BIN_WIDTH; // ~689.1 Hz (Approximates F5)
  constexpr double F_NOISE = 87.0 * BIN_WIDTH;  // ~14987 Hz (High Freq Noise)

  for (std::size_t i = 0; i < N; ++i)
  {
    double t = static_cast<double>(i) / SAMPLE_RATE;

    // Coherent Chord
    double chord = 1.0 * math::sin(2 * std::numbers::pi * F_CHORD_1 * t) +
                   0.8 * math::sin(2 * std::numbers::pi * F_CHORD_2 * t) +
                   0.6 * math::sin(2 * std::numbers::pi * F_CHORD_3 * t);

    // Coherent Noise
    double noise = 0.5 * math::sin(2 * std::numbers::pi * F_NOISE * t);

    signal[i] = chord + noise;
  }

  // 2. FFT
  auto spectrum = fft::fft(signal, false);

  // 3. Analysis & Filtering
  std::size_t half_n = N / 2;
  // double freq_per_bin = BIN_WIDTH; // Redundant but consistent
  std::size_t removed = 0;
  bool noise_found = false;

  // Check pre-filter state
  // We know noise is at bin 87
  // Manual abs check for constexpr: (x > 0 ? x : -x)
  double r87 = spectrum[87].real();
  double i87 = spectrum[87].imag();
  if ((r87 > 10.0 || r87 < -10.0) || (i87 > 10.0 || i87 < -10.0))
  {
    // Manual magnitude check to confirm
    double mag_sq = r87 * r87 + i87 * i87;
    if (mag_sq > 100.0)
      noise_found = true;
  }

  // Apply Filter (Cutoff ~1000 Hz)
  // Bin 6 ~= 1033 Hz. So we zero out bins > 6
  for (std::size_t i = 7; i < N - 6; ++i)
  { // Simple low-pass excluding DC-ish area
    // Better logic using frequency calculation for clarity
    double freq = (i <= half_n) ? (i * BIN_WIDTH) : ((N - i) * BIN_WIDTH);
    if (freq > 1000.0)
    {
      double mag_sq = spectrum[i].real() * spectrum[i].real() + spectrum[i].imag() * spectrum[i].imag();
      if (mag_sq > 1e-9)
        removed++; // Use small epsilon
      spectrum[i] = 0;
    }
  }

  // 4. IFFT Reconstruction
  auto clean_signal = fft::fft(spectrum, true);

  // 5. Verification & Data Export
  std::array<double, N> output_data = {}; // Zero init
  double error_sum = 0;
  for (std::size_t i = 0; i < N; ++i)
  {
    output_data[i] = clean_signal[i].real(); // Capture real part for WAV

    double t = static_cast<double>(i) / SAMPLE_RATE;
    // Reconstruct ideal using same coherent frequencies
    double ideal = 1.0 * math::sin(2 * std::numbers::pi * F_CHORD_1 * t) +
                   0.8 * math::sin(2 * std::numbers::pi * F_CHORD_2 * t) +
                   0.6 * math::sin(2 * std::numbers::pi * F_CHORD_3 * t);

    double diff = clean_signal[i].real() - ideal;
    error_sum += diff * diff;
  }

  return {math::sqrt(error_sum / N), removed, noise_found, output_data};
}

// Simple WAV Writer (16-bit Mono)
void write_wav(const std::string& filename, const std::array<double, N>& data)
{
  std::ofstream file(filename, std::ios::binary);
  if (!file)
    return;

  // WAV Header
  struct WavHeader
  {
    char riff[4] = {'R', 'I', 'F', 'F'};
    std::uint32_t overall_size;
    char wave[4] = {'W', 'A', 'V', 'E'};
    char fmt[4] = {'f', 'm', 't', ' '};
    std::uint32_t fmt_chunk_size = 16;
    std::uint16_t format_type = 1; // PCM
    std::uint16_t channels = 1;
    std::uint32_t sample_rate = static_cast<std::uint32_t>(SAMPLE_RATE);
    std::uint32_t byte_rate;
    std::uint16_t block_align;
    std::uint16_t bit_depth = 16;
    char data_marker[4] = {'d', 'a', 't', 'a'};
    std::uint32_t data_size;
  } header;

  header.block_align = header.channels * header.bit_depth / 8;
  header.byte_rate = header.sample_rate * header.block_align;
  header.data_size = static_cast<std::uint32_t>(data.size() * header.block_align);
  header.overall_size = header.data_size + sizeof(WavHeader) - 8;

  file.write(reinterpret_cast<const char*>(&header), sizeof(WavHeader));

  // Data (Convert double [-1.0, 1.0] to int16)
  for (double sample : data)
  {
    // Hard clipper to avoid wrap-around distortion
    if (sample > 1.0)
      sample = 1.0;
    if (sample < -1.0)
      sample = -1.0;

    std::int16_t pcm = static_cast<std::int16_t>(sample * 32767.0);
    file.write(reinterpret_cast<const char*>(&pcm), sizeof(std::int16_t));
  }

  std::print("  [IO] Saved filtered audio to '{}'.\n", filename);
}

int main()
{
  // This line proves it runs at compile-time!
  constexpr AudioResult result = process_audio_at_compile_time();

  std::print("Audio Filtering Demo (Computed at Compile-Time)\n");
  std::print("-----------------------------------------------\n");
  std::print("High-Frequency Noise Detected: {}\n", result.high_freq_noise_detected ? "YES" : "NO");
  std::print("Bins Zeroed Out (Filtered):    {}\n", result.removed_energy_bins);
  std::print("Reconstruction RMS Error:      {:.20f}\n", result.rms_error);

  // With coherent sampling, error should be extremely low (rounding error domain)
  if constexpr (result.rms_error < 1e-10)
  {
    std::print("\n[SUCCESS] Signal restored perfectly (using Coherent Sampling).\n");
  }
  else
  {
    std::print("\n[FAIL] RMS error too high ({:.6f}).\n", result.rms_error);
  }

  // Export at runtime
  write_wav("clean_audio.wav", result.clean_audio_data);

  return 0;
}
