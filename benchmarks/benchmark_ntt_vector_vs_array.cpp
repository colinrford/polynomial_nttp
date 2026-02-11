/*
 *  benchmark_ntt_vector_vs_array.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */
import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Solinas-like Prime: 4179340454199820289 = 29 * 2^57 + 1
constexpr auto local_prime = 4179340454199820289_Z;
using field = decltype(lam::cbn::Zq(local_prime));

// Traits specialization
namespace lam::polynomial::univariate
{
template<typename T, T... Modulus>
struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
{
  using K = lam::cbn::ZqElement<T, Modulus...>;
  static constexpr bool is_finite_field = true;
  static constexpr T modulus = []() {
    constexpr T mods[] = {Modulus...};
    return mods[0];
  }();
  static constexpr K mul(const K& a, const K& b) { return a * b; }
  static constexpr K add(const K& a, const K& b) { return a + b; }
  static constexpr K sub(const K& a, const K& b) { return a - b; }
};
} // namespace lam::polynomial::univariate

using namespace lam::polynomial::univariate::ntt;

template<std::size_t N>
void run_comparison()
{
  std::println("\n--- Degree {} ---", N);

  // 1. Vector Benchmark
  std::vector<field> vec_data(N);
  for (std::size_t i = 0; i < N; ++i)
    vec_data[i] = field(i);

  // Warmup
  ntt_transform(vec_data, false);
  ntt_transform(vec_data, true);

  auto start_vec = std::chrono::steady_clock::now();
  constexpr int iterations = 1000;
  for (int i = 0; i < iterations; ++i)
  {
    ntt_transform(vec_data, false);
    ntt_transform(vec_data, true);
    lam::polynomial::is_negligible(vec_data[0]);
  }
  auto end_vec = std::chrono::steady_clock::now();
  auto dur_vec = std::chrono::duration_cast<std::chrono::microseconds>(end_vec - start_vec).count();

  // 2. Array Benchmark
  static std::array<field, N> arr_data; // Static to avoid stack overflow for large N
  for (std::size_t i = 0; i < N; ++i)
    arr_data[i] = field(i);

  // Warmup
  ntt_transform(arr_data, false);
  ntt_transform(arr_data, true);

  auto start_arr = std::chrono::steady_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    ntt_transform(arr_data, false);
    ntt_transform(arr_data, true);
    lam::polynomial::is_negligible(arr_data[0]);
  }
  auto end_arr = std::chrono::steady_clock::now();
  auto dur_arr = std::chrono::duration_cast<std::chrono::microseconds>(end_arr - start_arr).count();

  double avg_vec = static_cast<double>(dur_vec) / iterations;
  double avg_arr = static_cast<double>(dur_arr) / iterations;

  std::println("  Vector: {:>8.2f} us", avg_vec);
  std::println("  Array:  {:>8.2f} us", avg_arr);
  std::println("  Speedup: {:>8.2f}x", avg_vec / avg_arr);
}

int main()
{
  std::println("Benchmarking NTT: std::vector vs std::array");

  run_comparison<64>();
  run_comparison<256>();
  run_comparison<1024>();
  run_comparison<4096>();
  // run_comparison<8192>(); // Might be large for stack if not static

  return 0;
}
