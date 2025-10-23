/*
 *  division.cpp – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp polynomial division
 *    need to return (placeholder circa Oct. 8, 2025)
 */

import polynomial_nttp;

using namespace math_nttp;
// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
constexpr bool comparitore(double a, double b)
{
  constexpr auto abs_val = [](double val) { return val > 0. ? val : -val; };
  if (abs_val(a - b) < 1e-7)
    return true;
  else
    return false;
}

consteval auto division_by_zero_result()
{
  constexpr polynomial_nttp<double, 2> a{0., 1., 2.};
  constexpr polynomial_nttp<double, 1> b{0., 0.};
  constexpr auto a_divided_by_b =
                        division_prototype<double, 2, a, 1, b>();
  return a_divided_by_b;
}

consteval bool divide_by_zero()
{
  constexpr auto q_and_r = division_by_zero_result();
  constexpr auto q_of_x = q_and_r.first;
  constexpr auto r_of_x = q_and_r.second;
  if (norm(q_of_x) != 0 || !comparitore(q_of_x[0], 0.))
    return false;
  else if (2 != norm(r_of_x))
    return false;
  else if (!comparitore(0., r_of_x[0])
        || !comparitore(1., r_of_x[1])
        || !comparitore(2., r_of_x[2]))
    return false;
  else
    return true;
}

consteval auto division_by_greater_degree_result()
{
  constexpr polynomial_nttp<double, 2> a{-1.2, 1.5, 2.2};
  constexpr polynomial_nttp<double, 1> b{1.5, 2.5};
  constexpr auto b_divided_by_a =
                        division_prototype<double, 1, b, 2, a>();
  return b_divided_by_a;
}

consteval bool divide_by_greater_degree()
{
  constexpr auto q_and_r = division_by_greater_degree_result();
  constexpr auto q_of_x = q_and_r.first;
  constexpr auto r_of_x = q_and_r.second;
  if (norm(q_of_x) != 0 || !comparitore(q_of_x[0], 0.))
    return false;
  else if (1 != norm(r_of_x))
    return false;
  else if (!comparitore(1.5, r_of_x[0])
        || !comparitore(2.5, r_of_x[1]))
    return false;
  else
    return true;
}

/*
 *  a polynomial divided by itself should have quotient = 1, remainder = 0
 */
consteval bool divide_by_self()
{
  constexpr polynomial_nttp<double, 2> quadratic{0., 1., 2.};
  constexpr auto q_and_r = division_prototype<double,
                                              2,
                                              quadratic,
                                              2,
                                              quadratic>();
  constexpr auto q_of_x = q_and_r.first;
  constexpr auto r_of_x = q_and_r.second;
  if (norm(r_of_x) != 0 || !comparitore(r_of_x[0], 0.))
    return false;
  else if (norm(q_of_x) != 0)
    return false;
  else if (!comparitore(1., q_of_x[0]))
    return false;
  else
    return true;
}

consteval bool self_squared_divided_by_self()
{
  constexpr polynomial_nttp<double, 2> quadratic{0., 1., 2.};
  constexpr auto quartic = quadratic * quadratic;
  constexpr auto q_and_r = division_prototype<double,
                                              4,
                                              quartic,
                                              2,
                                              quadratic>();
  constexpr auto q_of_x = q_and_r.first;
  constexpr auto r_of_x = q_and_r.second;
  if (norm(r_of_x) != 0 || !comparitore(r_of_x[0], 0.))
    return false;
  else if (norm(q_of_x) != 2)
    return false;
  else if (!comparitore(quadratic[0], q_of_x[0])
        || !comparitore(quadratic[1], q_of_x[1])
        || !comparitore(quadratic[2], q_of_x[2]))
    return false;
  else
    return true;
}

consteval bool reconstruct_original_polynomial()
{
  constexpr polynomial_nttp<double, 3> cubic{0., 1., 2., 3.};
  constexpr polynomial_nttp<double, 1> linear{-1., 1.};
  constexpr auto q_and_r = division_prototype<double, 3, cubic, 1, linear>();
  constexpr auto quotient = q_and_r.first;
  constexpr auto remainder = q_and_r.second;
  constexpr auto reconstructed_cubic = quotient * linear + remainder;
  if (norm(cubic) != norm(reconstructed_cubic))
    return false;
  else if (!comparitore(cubic[0], reconstructed_cubic[0])
        && !comparitore(cubic[1], reconstructed_cubic[1])
        && !comparitore(cubic[2], reconstructed_cubic[2])
        && !comparitore(cubic[3], reconstructed_cubic[3]))
    return false;
  else
    return true;
}
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

int main()
{
  constexpr bool divide_by_zero_works = divide_by_zero();
  constexpr bool divide_by_greater_degree_works = divide_by_greater_degree();
  constexpr bool divide_by_self_works = divide_by_self();
  constexpr bool self_squared_divided_by_self_works
                  = self_squared_divided_by_self();
  constexpr bool reconstruct_original_polynomial_works
                  = reconstruct_original_polynomial();

  if constexpr (divide_by_zero_works
             && divide_by_greater_degree_works
             && divide_by_self_works
             && self_squared_divided_by_self_works
             && reconstruct_original_polynomial_works)
  {
    return 0;
  } else
  {
    return 1;
  }
}
