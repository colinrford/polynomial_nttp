/*
 *  type_erased.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  not sure if this is interesting or not but 'twas an idea
 */

import std;
import polynomial_nttp;

namespace stdr = std::ranges;
namespace stdv = std::views;

auto indexing_set = [](auto n) {
  return stdv::iota(decltype(n)(0), n);
};

auto indexing_set_from_to = [](auto m, auto n) {
  return stdv::iota(decltype(n)(m), n);
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
    math_nttp::polynomial_nttp<double, two_n_plus_1> power_series{};
    for (std::size_t&& i : indexing_set(n))
    {
      double neg_one_i = neg_one_multiplier(i);
      double two_i_plus_one_fact = factorial(2 * i + 1);
      power_series.coefficients[2 * i + 1] = neg_one_i / two_i_plus_one_fact;
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
    math_nttp::polynomial_nttp<double, two_n> power_series{};
    for (std::size_t&& i : indexing_set(n))
    {
      double neg_one_i = neg_one_multiplier(i);
      double two_i_fact = factorial(2 * i);
      power_series.coefficients[2 * i] = neg_one_i / two_i_fact;
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
    math_nttp::polynomial_nttp<double, two_n> power_series{};
    for (std::size_t&& i : indexing_set(two_n))
    {
      double one_over_n_fact = 1 / factorial(i);
      power_series.coefficients[i] = one_over_n_fact;
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
  constexpr double epsilon = std::numeric_limits<double>::epsilon();
  double running_root = 1.;
  double most_recent_running_root = 1.;
  auto abs_val = [](double val) { return val > 0. ? val : -val; };
  do {
    most_recent_running_root = running_root;
    running_root = sqrt_iterator(of, running_root);
  } while (abs_val(running_root - most_recent_running_root) > epsilon);
  return running_root;
}

std::function<double(double)> which_type_erased_polynomial(int index)
{
  if (index == 0)
    return [](double d) { return std::sin<double>(d); };
  else if (index == 1)
    return [](double d) { return std::cos<double>(d); };
  else if (index == 2)
    return [](double d) { return std::exp<double>(d); };
  else if (index == 3)
    return [](double d) { return std::sqrt<double>(d); };
  else return {};
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
}
