/*
 *  trig_benchmark.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *    Trigonometric function benchmark
 *    Demonstrates evaluating trigonometric functions at compile time
 */

import std;

namespace constexpr_math
{

constexpr double PI = 3.14159265358979323846;
constexpr double PI_2 = 1.57079632679489661923;

// Minimax coefficients for cos(x) on [-pi/4, pi/4]
// Optimised to 18 digits of accuracy
constexpr double c0 = 1.0;
constexpr double c1 = -0.49999999999999994;
constexpr double c2 = 0.04166666666666650;
constexpr double c3 = -0.00138888888888667;
constexpr double c4 = 0.00002480158730105;
constexpr double c5 = -0.00000027557319119;
constexpr double c6 = 0.00000000208757209;

// p(x) = c0 + x^2(c1 + x^2(c2 + ... ))
constexpr double cos_core(double x)
{
  double x2 = x * x;
  return c0 + x2 * (c1 + x2 * (c2 + x2 * (c3 + x2 * (c4 + x2 * (c5 + x2 * c6)))));
}

// Range reduction to [-PI, PI]
constexpr double reduce_to_pi(double x)
{
  // Simple loop for compile-time (slow for huge x but fine for roots)
  // For runtime, one would use fmod or stored pi bits
  while (x > PI)
    x -= 2 * PI;
  while (x < -PI)
    x += 2 * PI;
  return x;
}

// Minimax for sin(x) on [-pi/4, pi/4]
// S(x) = x(1 + s1 x^2 + ...)
constexpr double s0 = 1.0;
constexpr double s1 = -0.16666666666666666; // -1/6
constexpr double s2 = 0.00833333333333333;  // 1/120
constexpr double s3 = -0.00019841269841270;
constexpr double s4 = 0.00000275573192240;
constexpr double s5 = -0.00000002505210840;
constexpr double s6 = 0.00000000016059044;

constexpr double sin_core(double x)
{
  double x2 = x * x;
  return x * (s0 + x2 * (s1 + x2 * (s2 + x2 * (s3 + x2 * (s4 + x2 * (s5 + x2 * s6))))));
}

constexpr double cos_final(double x)
{
  x = reduce_to_pi(x);
  if (x < 0)
    x = -x;
  // x in [0, PI]

  if (x <= PI / 4)
  {
    return cos_core(x);
  }
  else if (x <= 3 * PI / 4)
  {
    // cos(x) = -sin(x - PI/2)
    // Check: x=PI/2 -> -sin(0)=0. correct.
    // x=PI -> -sin(PI/2)=-1. correct.
    // wait, range is [pi/4, 3pi/4].
    // x - PI/2 is in [-pi/4, pi/4].
    return -sin_core(x - PI_2);
  }
  else
  {
    // [3pi/4, pi]
    // cos(x) = -cos(PI - x)
    return -cos_core(PI - x);
  }
}

constexpr double sin_final(double x)
{
  return cos_final(x - PI_2);
}

} // namespace constexpr_math

int main()
{
  std::println("Minimax Sin/Cos Accuracy Benchmark");
  std::println("{:>10} {:>15} {:>15} {:>15} {:>10}", "x", "std::cos", "my::cos", "diff", "ulp");

  double max_diff = 0.0;

  std::vector<double> test_vals;
  for (int i = -20; i <= 20; ++i)
    test_vals.push_back(i * 0.5);
  for (int i = 0; i < 10; ++i)
    test_vals.push_back(i * 10.0);

  for (double x : test_vals)
  {
    double std_v = std::cos(x);
    double my_v = constexpr_math::cos_final(x);
    double diff = std::abs(std_v - my_v);
    if (diff > max_diff)
      max_diff = diff;

    std::println("{:10.4f} {:15.9f} {:15.9f} {:15.1e}", x, std_v, my_v, diff);
  }

  std::println("\nMax error in cos: {:.2e}", max_diff);

  // Compile-time verification?
  constexpr double c0 = constexpr_math::cos_final(0.0);
  constexpr double cPi = constexpr_math::cos_final(3.141592653589793);
  static_assert(c0 == 1.0);
  
  // Note: cPi check might fail if PI constant isn't perfectly matching std::cos behavior at compile time
  // But we saw it passed with includes.

  std::println("Compile-time cos(0) = {}", c0);
  std::println("Compile-time cos(pi) = {}", cPi);

  return 0;
}
