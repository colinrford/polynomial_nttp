/*
 *  benchmark_threading_crossover.cpp
 *
 *    Determines the crossover point where threaded addition (std::jthread)
 *    becomes faster than naive serial addition for int types.
 */
#include <algorithm>
#include <chrono>
#include <numeric>
#include <print>
#include <thread>
#include <vector>

// Naive Serial Addition
template<typename T>
void naive_add(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& res)
{
  std::size_t n = a.size();
  for (std::size_t i = 0; i < n; ++i)
  {
    res[i] = a[i] + b[i];
  }
}

// Threaded Addition (Copy of implementation logic)
template<typename T>
void threaded_add(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& res)
{
  std::size_t N_size = a.size();
  std::size_t num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0)
    num_threads = 2;

  std::size_t chunk_size = N_size / num_threads;
  std::vector<std::jthread> threads;
  threads.reserve(num_threads);

  auto worker = [&](std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; ++i)
    {
      res[i] = a[i] + b[i];
    }
  };

  for (std::size_t t = 0; t < num_threads; ++t)
  {
    std::size_t start = t * chunk_size;
    std::size_t end = (t == num_threads - 1) ? N_size : start + chunk_size;
    threads.emplace_back(worker, start, end);
  }
}

// TBB Addition
#ifdef LAM_USE_TBB
#include <tbb/parallel_for.h>
template<typename T>
void tbb_add(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& res)
{
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, a.size()), [&](const tbb::blocked_range<std::size_t>& r) {
    for (std::size_t i = r.begin(); i != r.end(); ++i)
    {
      res[i] = a[i] + b[i];
    }
  });
}
#endif

void run_benchmark(std::size_t N)
{
  using T = int;
  std::vector<T> a(N, 1), b(N, 2), res(N);

  // Warmup & Iteration Count Heuristic
  int iterations = 1000;
  if (N > 10000)
    iterations = 100;
  if (N > 100000)
    iterations = 10;
  if (N < 1000)
    iterations = 10000;

  // Measure Naive
  auto start_naive = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    naive_add(a, b, res);
    volatile auto sink = res[0];
  }
  auto end_naive = std::chrono::high_resolution_clock::now();

  auto dur_naive =
    std::chrono::duration_cast<std::chrono::nanoseconds>(end_naive - start_naive).count() / (double)iterations;

  // Measure Threaded (jthread)
  auto start_thread = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    threaded_add(a, b, res);
    volatile auto sink = res[0];
  }
  auto end_thread = std::chrono::high_resolution_clock::now();

  auto dur_thread =
    std::chrono::duration_cast<std::chrono::nanoseconds>(end_thread - start_thread).count() / (double)iterations;

#ifdef LAM_USE_TBB
  // Measure TBB
  auto start_tbb = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < iterations; ++i)
  {
    tbb_add(a, b, res);
    volatile auto sink = res[0];
  }
  auto end_tbb = std::chrono::high_resolution_clock::now();

  auto dur_tbb = std::chrono::duration_cast<std::chrono::nanoseconds>(end_tbb - start_tbb).count() / (double)iterations;

  std::print("N={:<7} | Naive: {:>8.0f} ns | jThread: {:>8.0f} ns | TBB: {:>8.0f} ns | Win: {}\n", N, dur_naive,
             dur_thread, dur_tbb,
             (dur_tbb < dur_naive && dur_tbb < dur_thread) ? "TBB" : (dur_thread < dur_naive ? "jThread" : "Naive"));
#else
  std::print("N={:<7} | Naive: {:>8.0f} ns | jThread: {:>8.0f} ns | Win: {}\n", N, dur_naive, dur_thread,
             (dur_thread < dur_naive ? "jThread" : "Naive"));
#endif
}

int main()
{
  std::print("=== Threading Crossover Benchmark ===\n");
  std::vector<size_t> sizes = {1000, 5000, 10000, 20000, 50000, 100000, 500000, 1000000};

  for (auto n : sizes)
  {
    run_benchmark(n);
  }
  return 0;
}
