/*
 *  benchmark_ntt_vs_ntl.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Solinas Prime: 18446744069414584321 = 2^64 - 2^32 + 1
constexpr auto local_prime = 18446744069414584321ULL;
using field = lam::cbn::ZqElement<std::uint64_t, local_prime>;

// Forward declarations of external C++ implementation
extern "C++" void run_ntl_benchmark_impl(long N);
extern "C++" void init_ntl_modulus_impl(const char* prime_str);
extern "C++" bool verify_ntl_multiplication_impl(long N, const long* a_data, const long* b_data,
                                                 const long* expected_data);

// ... (existing includes and formatting) ...

template<std::size_t N>
void run_benchmark_local()
{
  using poly = lam::polynomial::univariate::polynomial_nttp<field, N>;
  poly a, b;
  // Fill with some data
  // Use long for easy bridging to NTL
  std::vector<long> a_vec(N + 1), b_vec(N + 1);

  for (std::size_t i = 0; i <= N; ++i)
  {
    a.coefficients[i] = field(i + 1);
    b.coefficients[i] = field(N - i);
    a_vec[i] = static_cast<long>(i + 1);
    b_vec[i] = static_cast<long>(N - i);
  }

  // Warmup & Correctness Check
  auto c = a * b;

  std::vector<long> res_vec(2 * N + 1);
  for (std::size_t i = 0; i <= 2 * N; ++i)
  {
    // Handle potential indexing depending on result size type
    if (i < c.coefficients.size())
      res_vec[i] = static_cast<long>(c.coefficients[i].data[0]);
    else
      res_vec[i] = 0;
  }

  bool match = verify_ntl_multiplication_impl(N, a_vec.data(), b_vec.data(), res_vec.data());
  if (match)
  {
    std::println("  [Correctness] Local result MATCHES NTL result.");
  }
  else
  {
    std::println("  [Incorrect] Local result DOES NOT MATCH NTL result!");
    std::exit(1);
  }

  // Performance Benchmark
  auto start = std::chrono::steady_clock::now();
  constexpr int iterations = 100;
  for (int i = 0; i < iterations; ++i)
  {
    auto res = a * b;
    lam::polynomial::is_negligible(res[0]); // Prevent optimize out
  }
  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  std::println("  Local (deg {}):     {:>8.2f} us", N, static_cast<double>(duration) / iterations);
}
namespace lam::polynomial::univariate
{
template<typename T, T... Modulus>
struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
{
  using K = lam::cbn::ZqElement<T, Modulus...>;
  static constexpr bool is_finite_field = true;
  static constexpr std::uint64_t modulus = 18446744069414584321ULL;

  static constexpr K mul(const K& a, const K& b) {
      unsigned __int128 prod = static_cast<unsigned __int128>(a.data[0]) * b.data[0];
      return K(lam::polynomial::univariate::ntt::reduce_solinas(prod));
  }
  static constexpr K add(const K& a, const K& b) { return a + b; }
  static constexpr K sub(const K& a, const K& b) { return a - b; }
};
} // namespace lam::polynomial::univariate


int main()
{
  std::println("Benchmarking NTT Multiplication: Local vs NTL");

  // Initialize NTL Modulus via wrapper (Solinas Prime)
  init_ntl_modulus_impl("18446744069414584321");

  std::println("\n--- Degree 1024 ---");
  run_benchmark_local<1024>();
  run_ntl_benchmark_impl(1024);

  std::println("\n--- Degree 2048 ---");
  run_benchmark_local<2048>();
  run_ntl_benchmark_impl(2048);

  std::println("\n--- Degree 4096 ---");
  run_benchmark_local<4096>();
  run_ntl_benchmark_impl(4096);

  std::println("\n--- Degree 8192 ---");
  run_benchmark_local<8192>();
  run_ntl_benchmark_impl(8192);

  std::println("\n--- Degree 16384 (Parallel Crossover) ---");
  run_benchmark_local<16384>();
  run_ntl_benchmark_impl(16384);

  std::println("\n--- Degree 32768 ---");
  run_benchmark_local<32768>();
  run_ntl_benchmark_impl(32768);

  std::println("\n--- Degree 65536 ---");
  run_benchmark_local<65536>();
  run_ntl_benchmark_impl(65536);

  return 0;
}
