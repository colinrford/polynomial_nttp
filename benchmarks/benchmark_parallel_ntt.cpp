/*
 *  benchmark_parallel_ntt.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License
 *
 */
 
import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;

using namespace lam::polynomial;

// Define field K
using field = lam::cbn::ZqElement<std::uint64_t, 18446744069414584321ULL>; // Solinas Prime 2^64 - 2^32 + 1

void benchmark(std::size_t N)
{
  if (N == 0)
    return;

  // Setup Data
  std::vector<field> data(N);
  std::mt19937_64 rng(12345);
  std::uniform_int_distribution<std::uint64_t> dist;
  for (auto& x : data)
    x = field(dist(rng));

  // Warmup
  auto data_copy = data;
  lam::polynomial::univariate::ntt::ntt_transform(data_copy, false);

  // Benchmark
  constexpr int ITERATIONS = 100;
  auto start = std::chrono::steady_clock::now();

  for (int i = 0; i < ITERATIONS; ++i)
  {
    data_copy = data; // Reset
    lam::polynomial::univariate::ntt::ntt_transform(data_copy, false);
  }

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double, std::nano> diff = end - start;
  double avg_ns = diff.count() / ITERATIONS;
  double ns_per_n = avg_ns / (N * std::log2(N)); // Normalized metric

  std::println("N={:<8} Time={:<10.0f} ns   Norm={:<8.2f} ns/NlogN", N, avg_ns, ns_per_n);
}

int main()
{
  std::println("Benchmarking Parallel NTT Thresholds...");

  // Powers of 2 around the thresholds (16384, 65536)
  std::vector<std::size_t> sizes = {4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576};

  for (auto n : sizes)
  {
    benchmark(n);
  }

  return 0;
}
