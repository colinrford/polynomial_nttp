/*
 *  embed_wav_demo.cpp
 *  Demonstrates loading a WAV file at compile-time using #embed and processing it.
 */

import std;
import lam.polynomial_nttp;

using namespace lam::polynomial::univariate;

// Helper: Read u32 little endian from bytes
consteval std::uint32_t read_u32(const unsigned char* bytes, std::size_t offset)
{
  return static_cast<std::uint32_t>(bytes[offset]) | (static_cast<std::uint32_t>(bytes[offset + 1]) << 8) |
         (static_cast<std::uint32_t>(bytes[offset + 2]) << 16) | (static_cast<std::uint32_t>(bytes[offset + 3]) << 24);
}

// Helper: Read u16 little endian
consteval std::uint16_t read_u16(const unsigned char* bytes, std::size_t offset)
{
  return static_cast<std::uint16_t>(bytes[offset]) | (static_cast<std::uint16_t>(bytes[offset + 1]) << 8);
}

// Helper: Read i16 little endian
consteval std::int16_t read_i16(const unsigned char* bytes, std::size_t offset)
{
  // Cast to unsigned first to define shift behavior, then cast back
  std::uint16_t val = static_cast<std::uint16_t>(bytes[offset]) | (static_cast<std::uint16_t>(bytes[offset + 1]) << 8);
  return std::bit_cast<std::int16_t>(val);
}

// -----------------------------------------------------------------------------
// The "Magic" Link to the File
// -----------------------------------------------------------------------------


// Load the exact bytes of the WAV file into the compiler's memory
// We use a lambda to process it immediately to avoid keeping the huge array in symbol table
consteval auto parse_embedded_wav()
{
  // Clang Extension: #embed
  // Note: We suppress warning via pragma in real projects, here we just show it.
  constexpr unsigned char wav_bytes[] = {
#embed "clean_audio.wav"
  };

  // 1. Parse Header
  // RIFF is at 0, WAVE at 8, fmt at 12
  // Data marker usually at 36, but can vary if extra chunks exist.
  // For our specific generated file, we know the layout is standard PCM.

  // Header size for standard PCM is 44 bytes
  constexpr std::size_t HEADER_SIZE = 44;

  // Safety Check: Verify it's a minimal valid WAV
  // In a real parser we'd check "RIFF", "WAVE" magic bytes.
  if (wav_bytes[0] != 'R' || wav_bytes[1] != 'I')
    throw "Invalid WAV";

  constexpr std::uint32_t data_size = read_u32(wav_bytes, 40);
  constexpr std::uint32_t sample_count = data_size / 2; // 16-bit

  // 2. Parse Samples
  std::vector<std::complex<double>> signal;
  signal.reserve(sample_count);

  for (std::size_t i = 0; i < sample_count; ++i)
  {
    std::size_t offset = HEADER_SIZE + i * 2;
    if (offset + 1 >= sizeof(wav_bytes))
      break; // Safety

    std::int16_t sample_i16 = read_i16(wav_bytes, offset);

    // Normalize int16 [-32768, 32767] to double [-1.0, 1.0]
    double sample_norm = static_cast<double>(sample_i16) / 32768.0;
    signal.push_back({sample_norm, 0.0});
  }

  return signal;
}

// Main Compile-Time Processing
consteval double calculate_signal_energy_at_compile_time()
{
  // STEP 1: LOAD
  auto signal = parse_embedded_wav();

  // STEP 2: FFT (Analysis)
  auto spectrum = fft::fft(signal, false);

  // STEP 3: COMPUTE STATS (Parseval's Thm)
  // Energy in time domain == Energy in freq domain / N
  double energy_freq = 0;
  for (auto val : spectrum)
  {
    energy_freq += val.real() * val.real() + val.imag() * val.imag();
  }

  double N = static_cast<double>(spectrum.size());
  return energy_freq / N;
}

int main()
{
  // This executes the entire pipeline:
  // 1. Embed "clean_audio.wav"
  // 2. Parse it
  // 3. FFT it
  // 4. Sum energy
  // ... all before the program starts!
  constexpr double total_energy = calculate_signal_energy_at_compile_time();

  std::print("Compile-Time WAV Analysis\n");
  std::print("-------------------------\n");
  std::print("Loaded 'clean_audio.wav' via #embed.\n");
  std::print("Total Signal Energy (calculated by compiler): {:.6f}\n", total_energy);
}
