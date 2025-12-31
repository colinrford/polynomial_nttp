/*
 *  type_erased.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  not sure if this is interesting or not but 'twas an idea
 */

import std;
import lam.polynomial_nttp;

namespace stdr = std::ranges;
namespace stdv = std::views;
// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
const auto indexing_set = [](auto n) {
  return stdv::iota(static_cast<decltype(n)>(0), n);
};

const auto indexing_set_from_to = [](auto m, auto n) {
  return stdv::iota(static_cast<decltype(n)>(m), n);
};

constexpr double factorial(std::size_t n)
{ return n == 0 ? 1. : static_cast<double>(n) * factorial(n - 1); }

constexpr double neg_one_multiplier(std::size_t index_i)
{ return index_i == 0 ? 1. : -1. * neg_one_multiplier(index_i - 1); }

constexpr auto sin_impl_test()
{
  return []() {
    constexpr std::size_t n = 12;
    constexpr std::size_t two_n_plus_1 = 2 * n + 1;
    lam::polynomial_nttp<double, two_n_plus_1> power_series{};
    for (std::size_t&& i : indexing_set(n))
    {
      double neg_one_i = neg_one_multiplier(i);
      double two_i_plus_one_fact = factorial(2 * i + 1);
      power_series.coefficients.at(2 * i + 1) = neg_one_i / two_i_plus_one_fact;
    }
    return power_series;
  }();
}

constexpr double sin_test(double x)
{ return sin_impl_test()(x); }

constexpr auto cos_impl_test()
{
  return []() {
    constexpr std::size_t n = 12;
    constexpr std::size_t two_n = 2 * n;
    lam::polynomial_nttp<double, two_n> power_series{};
    for (std::size_t&& i : indexing_set(n))
    {
      double neg_one_i = neg_one_multiplier(i);
      double two_i_fact = factorial(2 * i);
      power_series.coefficients.at(2 * i) = neg_one_i / two_i_fact;
    }
    return power_series;
  }();
}

constexpr double cos_test(double x)
{ return cos_impl_test()(x); }

constexpr auto exp_impl_test()
{
  return []() {
    constexpr std::size_t n = 16;
    constexpr std::size_t two_n = 30;
    lam::polynomial_nttp<double, two_n> power_series{};
    for (std::size_t&& i : indexing_set(two_n))
    {
      double one_over_n_fact = 1 / factorial(i);
      power_series.coefficients.at(i) = one_over_n_fact;
    }
    return power_series;
  }();
}

constexpr double exp_test(double x)
{ return exp_impl_test()(x); }

constexpr double sqrt_iterator(double of, double x_k)
{
  constexpr auto one_half = 0.5;
  return one_half * (x_k + of / x_k);
}

constexpr double sqrt_test(double of)
{
  if (of < 0.)
    return std::nan("");
  if (of == 0.)
    return 0.;
  double running_root = sqrt_iterator(of, 1.);
  double prev = 1.;
  // Iterate until we reach the Newton-Raphson fixed point (no change)
  while (running_root != prev)
  {
    prev = running_root;
    running_root = sqrt_iterator(of, running_root);
  }
  // If the result squared is less than 'of', it's an underestimate.
  // Bump up by 1 ULP to match correctly-rounded behavior.
  if (running_root * running_root < of)
  {
    auto bits = std::bit_cast<std::uint64_t>(running_root);
    return std::bit_cast<double>(bits + 1);
  }
  return running_root;
}

// Fast sqrt using bit manipulation for initial guess.
// sqrt(x) has roughly half the exponent of x in IEEE 754 representation.
// This gives us a ~12-bit accurate initial guess, so we only need ~4 iterations.
constexpr double sqrt_fast(double of)
{
  if (of < 0.)
    return std::nan("");
  if (of == 0.)
    return 0.;
  
  // Magic constant for double precision sqrt initial guess.
  // The idea: for IEEE 754 doubles, (bits >> 1) + magic ≈ sqrt.
  // 0x1FF7A3BEA91D9B1B is tuned for minimal initial error.
  constexpr std::uint64_t magic = 0x1FF7A3BEA91D9B1BULL;
  auto bits = std::bit_cast<std::uint64_t>(of);
  double guess = std::bit_cast<double>((bits >> 1) + magic);
  
  // Only ~4 Newton-Raphson iterations needed with this good initial guess
  guess = sqrt_iterator(of, guess);
  guess = sqrt_iterator(of, guess);
  guess = sqrt_iterator(of, guess);
  guess = sqrt_iterator(of, guess);
  
  // Final ULP refinement for correct rounding
  if (guess * guess < of)
  {
    auto result_bits = std::bit_cast<std::uint64_t>(guess);
    return std::bit_cast<double>(result_bits + 1);
  }
  return guess;
}

std::function<double(double)> which_type_erased_polynomial(int index)
{
  switch (index)
  {
    case 0:
      return [](double d) { return std::sin<double>(d); };
    case 1:
      return [](double d) { return std::cos<double>(d); };
    case 2:
      return [](double d) { return std::exp<double>(d); };
    case 3:
      return [](double d) { return std::sqrt<double>(d); };
    default:
      return {};
  }
}

int main()
{
  constexpr double sin_of_two = sin_test(2.);
  std::println("power series of sin:");
  std::println("\tsin_test(2) = {}", sin_of_two);
  std::println("\tstd::sin(2) = {}", std::sin(2.));
  std::println("\tdifference = {}", sin_of_two - std::sin(2.));

  constexpr double cos_of_two = cos_test(2.);
  std::println("power series of cos:");
  std::println("\tcos_test(2) = {}", cos_of_two);
  std::println("\tstd::cos(2) = {}", std::cos(2.));
  std::println("\tdifference = {}", cos_of_two - std::cos(2.));

  constexpr double exp_val = 1.;
  constexpr double exp_of_val = exp_test(exp_val);
  std::println("power series of exp:");
  std::println("\texp_test({}) = {}", exp_val, exp_of_val);
  std::println("\tstd::exp({}) = {}", exp_val, std::exp(exp_val));
  std::println("\tdifference = {}", exp_of_val - std::exp(exp_val));

  constexpr double sqrt_val = 2.;
  constexpr double sqrt_of_val = sqrt_test(sqrt_val);
  std::println("iterative sqrt:");
  std::println("\tsqrt_test({}) = {}", sqrt_val, sqrt_of_val);
  std::println("\tstd::sqrt({}) = {}", sqrt_val, std::sqrt(sqrt_val));
  std::println("\tdifference = {}", sqrt_of_val - std::sqrt(sqrt_val));

  std::vector<std::function<double(double)>> type_erased_polynomials{};
  type_erased_polynomials.push_back(sin_test);
  type_erased_polynomials.push_back(cos_test);
  type_erased_polynomials.push_back(exp_test);
  type_erased_polynomials.push_back(sqrt_test); // not actually implemented as
                                                // polynomial / power series

  int type_erased_counter = 0;
  constexpr int last = 100;
  constexpr double one_half = 0.5;
  for (auto&& p : type_erased_polynomials)
  {
    std::println("type_erased_counter =\t{}", type_erased_counter);
    const auto current_test_function = which_type_erased_polynomial(
                                         type_erased_counter
                                       );
    for (int i = 0; i < last; ++i)
    {
      double increment = -one_half
                        + static_cast<double>(last - i)
                        / static_cast<double>(last);
      if (std::isnan(p(increment)))
      {
        std::print("\tencountered NaN at {}, ", i);
        std::println("breaking iteration over interval");
        break;
      }
      auto diff = std::abs(current_test_function(increment) - p(increment));
      if (diff > 0) // only print nonzero differences
      {
        std::println("\tpolynomial: {};\tposition: {};\tdifference: {}",
                      type_erased_counter,
                      increment,
                      diff);
      }
    }
    ++type_erased_counter;
    if (type_erased_counter != type_erased_polynomials.size())
      std::println("\tswitch function");
  }

  // Quick test: compare sqrt_test vs std::sqrt for 1000 rational numbers in [1, 50]
  constexpr int N = 1000;
  constexpr double range_start = 1.0;
  constexpr double range_end = 50.0;
  int exact_matches = 0;
  int ulp_1_off = 0;
  double max_diff = 0.0;
  double max_diff_at = 0.0;
  
  std::println("\n--- sqrt comparison test ({} rationals in [{}, {}]) ---", N, range_start, range_end);
  for (int i = 0; i < N; ++i)
  {
    double val = range_start + (range_end - range_start) * static_cast<double>(i) / static_cast<double>(N - 1);
    double mine = sqrt_test(val);
    double standard = std::sqrt(val);
    double diff = std::abs(mine - standard);
    
    if (diff == 0.0)
      ++exact_matches;
    else if (diff <= std::numeric_limits<double>::epsilon() * standard)
      ++ulp_1_off;
    
    if (diff > max_diff)
    {
      max_diff = diff;
      max_diff_at = val;
    }
  }
  
  std::println("Exact matches: {}/{}", exact_matches, N);
  std::println("Within 1 ULP:  {}/{}", exact_matches + ulp_1_off, N);
  std::println("Max difference: {} (at sqrt({}))", max_diff, max_diff_at);

  // Speed comparison: time both implementations over many iterations
  constexpr int iterations = 1'000'000;
  constexpr int num_values = 100;
  std::array<double, num_values> test_values{};
  for (int i = 0; i < num_values; ++i)
    test_values[i] = 1.0 + 49.0 * static_cast<double>(i) / static_cast<double>(num_values - 1);

  volatile double sink = 0.0; // prevent optimizer from eliminating the loop

  std::println("\n--- sqrt speed comparison ({} iterations x {} values) ---", iterations, num_values);

  // Time sqrt_test
  auto start_mine = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& val : test_values)
      sink = sqrt_test(val);
  auto end_mine = std::chrono::steady_clock::now();
  auto duration_mine = std::chrono::duration_cast<std::chrono::milliseconds>(end_mine - start_mine);

  // Time sqrt_fast
  auto start_fast = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& val : test_values)
      sink = sqrt_fast(val);
  auto end_fast = std::chrono::steady_clock::now();
  auto duration_fast = std::chrono::duration_cast<std::chrono::milliseconds>(end_fast - start_fast);

  // Time std::sqrt
  auto start_std = std::chrono::steady_clock::now();
  for (int iter = 0; iter < iterations; ++iter)
    for (const auto& val : test_values)
      sink = std::sqrt(val);
  auto end_std = std::chrono::steady_clock::now();
  auto duration_std = std::chrono::duration_cast<std::chrono::milliseconds>(end_std - start_std);

  std::println("sqrt_test: {} ms", duration_mine.count());
  std::println("sqrt_fast: {} ms", duration_fast.count());
  std::println("std::sqrt: {} ms", duration_std.count());
  std::println("Ratio sqrt_test/std: {:.2f}x", 
               static_cast<double>(duration_mine.count()) / static_cast<double>(duration_std.count()));
  std::println("Ratio sqrt_fast/std: {:.2f}x", 
               static_cast<double>(duration_fast.count()) / static_cast<double>(duration_std.count()));

  // Quick accuracy check for sqrt_fast
  int fast_exact = 0;
  for (int i = 0; i < N; ++i)
  {
    double val = range_start + (range_end - range_start) * static_cast<double>(i) / static_cast<double>(N - 1);
    if (sqrt_fast(val) == std::sqrt(val))
      ++fast_exact;
  }
  std::println("\nsqrt_fast exact matches: {}/{}", fast_exact, N);
} // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
