/*
 *  root_finding_demo.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Demonstrates compile-time (constexpr) and runtime root finding
 *    for Number Theoretic Transforms using ctbignum interop.
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::polynomial::univariate::ntt;

constexpr std::uint64_t P = 18446744069414584321ULL;
using Zq = lam::cbn::ZqElement<std::uint64_t, P>;

constexpr std::uint64_t N_compile = 4096;
constexpr std::uint64_t root_compile = find_nth_root_of_unity(N_compile, P);

static_assert(root_compile != 0, "Root must exist for Solinas prime");
static_assert(power(root_compile, N_compile, P) == 1);
static_assert(power(root_compile, N_compile / 2, P) != 1);

void demo_runtime_root_finding()
{
  std::println("--- Runtime Root Finding Demo ---");

  std::size_t N = 1048576;
  std::println("Finding {}-th root of unity (mod {})...", N, P);

  auto start = std::chrono::high_resolution_clock::now();
  std::uint64_t root = find_nth_root_of_unity(static_cast<std::uint64_t>(N), P);
  auto end = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

  if (root != 0)
  {
    std::println("  Found root: {:x}", root);
    std::println("  Time taken: {} us", duration.count());

    Zq root_zq(root);
    std::println("  Verified in Zq: root^N = {:x}", (std::uint64_t)(power(root, static_cast<std::uint64_t>(N), P)));
  }
  else
  {
    std::println("  No root found for N={}!", N);
  }
}

int main()
{
  std::println("Compile-time root finding for N=4096: SUCCESS");
  std::println("  Root: {:x}", root_compile);

  demo_runtime_root_finding();

  return 0;
}
