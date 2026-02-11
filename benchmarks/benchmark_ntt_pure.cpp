/*
 *  benchmark_ntt_pure.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */
import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Solinas-like Prime: 4179340454199820289 = 29 * 2^57 + 1
constexpr auto ntt_prime = 4179340454199820289_Z;
using field = decltype(lam::cbn::Zq(ntt_prime));
using wrapper = lam::polynomial::univariate::finite_field_traits<field>;

// Must specialize finite_field_traits so NTT engine recognizes ZqElement
namespace lam::polynomial::univariate
{
template<typename T, T... Modulus>
struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
{
  static constexpr bool is_finite_field = true;
  static constexpr T modulus = []() {
    constexpr T mods[] = {Modulus...};
    return mods[0];
  }();
};
} // namespace lam::polynomial::univariate

using namespace lam::polynomial::univariate::ntt;

template<std::size_t N>
void run_benchmark()
{
  std::vector<field> data(N);
  for (std::size_t i = 0; i < N; ++i)
  {
    data[i] = field(static_cast<long>((i * 12345) ^ 0xDEADBEEF));
  }

  auto start = std::chrono::steady_clock::now();
  constexpr int iterations = 1000;

  // Warmup
  ntt_transform(data, false);
  ntt_transform(data, true);

  for (int i = 0; i < iterations; ++i)
  {
    ntt_transform(data, false);
    // Don't inverse every time to keep it hot, just measure forward transform density
    // But to prevent overflow/garbage, maybe reset? No, modular arithmetic is fine.
    // Actually, let's measure a round trip pair to be realistic.
    ntt_transform(data, true);
  }

  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  double avg_us = static_cast<double>(duration) / iterations; // microseconds per Round Trip
  double avg_ns_per_element = (avg_us * 1000.0) / N;          // nanoseconds per scalar

  std::println("N={:<5} | Round-Trip: {:<8.3f} us | Throughput: {:<6.1f} ns/elem", N, avg_us, avg_ns_per_element);
}

int main()
{
  std::println("Benchmarking Pure NTT (Forward + Inverse)");
  std::println("Prime P = 29*2^57 + 1 (64-bit)");
  std::println("---------------------------------------------------------------");

  run_benchmark<64>();
  run_benchmark<128>();
  run_benchmark<256>();
  run_benchmark<512>();
  run_benchmark<1024>();
  run_benchmark<2048>();
  run_benchmark<4096>();
  run_benchmark<8192>();

  return 0;
}
