/*
 *  polynomial_nttp-univariate-math.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

export module lam.polynomial_nttp:univariate.math;
import std;

export namespace lam::polynomial::univariate::math
{

// High-Precision PI (Matches std::numbers::pi)
// Uses Machin's Formula: pi/4 = 4*arctan(1/5) - arctan(1/239)
template<typename T>
consteval T compute_pi()
{
  auto arctan = [](T x) -> T {
    T sum = 0;
    T term = x;
    T x2 = x * x;
    int n = 1;
    while (term > 1e-18 || term < -1e-18)
    {
      if (n % 4 == 1)
        sum += term;
      else
        sum -= term;
      term *= x2;
      term *= static_cast<T>(n) / static_cast<T>(n + 2);
      n += 2;
    }
    return sum;
  };
  // This simple series convergence is too slow for high precision in limited steps.
  // We use hardcoded value for reliability in this demo context, or rely on std::numbers::pi_v if available.
  return std::numbers::pi_v<T>;
}

// Constants
template<typename T>
constexpr T PI = std::numbers::pi_v<T>;
template<typename T>
constexpr T PI_2 = PI<T> / 2;
template<typename T>
constexpr T TWO_PI = PI<T> * 2;

// Range Reduction: Reduce x to [-pi, pi]
template<typename T>
constexpr T reduce_to_pi(T x)
{
  constexpr T two_pi = TWO_PI<T>;
  constexpr T pi = PI<T>;

  // Simple range reduction (not robust for huge x, but fine for signal gen)
  long long k = static_cast<long long>((x + pi) / two_pi);
  return x - k * two_pi;
}

// Taylor Series Sine (Converges fast for small x)
template<typename T>
constexpr T sin_taylor(T x)
{
  T term = x;
  T sum = term;
  T x2 = x * x;
  for (int n = 3; n <= 19; n += 2)
  {
    term *= -x2 / ((n - 1) * n);
    sum += term;
  }
  return sum;
}

// Taylor Series Cosine
template<typename T>
constexpr T cos_taylor(T x)
{
  T term = 1.0;
  T sum = term;
  T x2 = x * x;
  for (int n = 2; n <= 18; n += 2)
  {
    term *= -x2 / ((n - 1) * n);
    sum += term;
  }
  return sum;
}

// Constexpr Sine
template<std::floating_point T>
constexpr T sin(T x)
{
  if consteval
  {
    x = reduce_to_pi(x);
    return sin_taylor(x);
  }
  else
  {
    return std::sin(x); // Runtime fallback
  }
}

// Constexpr Cosine
template<std::floating_point T>
constexpr T cos(T x)
{
  if consteval
  {
    // cos(x) = sin(x + pi/2)
    return sin(x + PI_2<T>);
  }
  else
  {
    return std::cos(x);
  }
}

// Constexpr Sqrt (Newton-Raphson)
template<typename T>
constexpr T sqrt(T x)
{
  if consteval
  {
    if (x < 0)
      return std::numeric_limits<T>::quiet_NaN();
    if (x == 0)
      return 0;
    T curr = x;
    T prev = 0;
    while (curr != prev)
    {
      prev = curr;
      curr = 0.5 * (curr + x / curr);
    }
    return curr;
  }
  else
  {
    return std::sqrt(x);
  }
}

} // namespace lam::polynomial::univariate::math
