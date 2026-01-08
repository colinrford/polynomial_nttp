/*
 *  roots_finite_field.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  unit test for polynomial_nttp roots function over Finite Fields (ctbignum)
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;
import lam.interop;

using namespace lam::cbn;
using namespace lam::cbn::literals;
using lam::polynomial::univariate::polynomial_nttp;
using lam::roots;

// Helper to check if roots match expected values (order independent)
template<typename Roots, typename Expected>
constexpr bool check_roots(const Roots& roots, const Expected& expected)
{
  if (roots.size() != expected.size())
    return false;

  auto local_roots = roots.to_vector();
  auto local_expected = expected;

  // Simple O(N^2) match since N is tiny
  for (const auto& exp_val : local_expected)
  {
    bool found = false;
    for (auto it = local_roots.begin(); it != local_roots.end(); ++it)
    {
      if (it->value == exp_val)
      {
        found = true;
        local_roots.erase(it);
        break;
      }
    }
    if (!found)
      return false;
  }
  return true;
}

// ============================================================
// Test 1: Linear ax + b = 0 in GF(11)
// ============================================================
bool test_linear()
{
  using GF = decltype(Zq(11_Z));
  // 4x + 9 = 6  -> 4x = -3 = 8 -> x = 2
  // 4*2 + 9 = 8 + 9 = 17 = 6 (mod 11). Correct.
  // Polynomial: 4x + (9-6) = 4x + 3 = 0
  constexpr GF four(4_Z);
  constexpr GF three(3_Z);

  polynomial_nttp<GF, 1> p{{three, four}};
  auto r = roots(p);

  std::vector<GF> expected{GF(2_Z)};
  if (!check_roots(r, expected))
  {
    std::println("Linear GF(11) Failed");
    return false;
  }
  std::println("Linear GF(11) Passed");
  return true;
}

// ============================================================
// Test 2: Quadratic x^2 - 3 = 0 in GF(13)
// ============================================================
// Squares in GF(13): 0, 1, 4, 9, 3, 12, 10
// 3 is a square: 4^2 = 16 = 3, 9^2 = 81 = 3.
// Roots should be 4 and 9.
bool test_quadratic()
{
  using GF = decltype(Zq(13_Z));
  constexpr GF one(1_Z);
  constexpr GF minus_three(10_Z); // -3 = 10 mod 13

  polynomial_nttp<GF, 2> p{{minus_three, GF(0_Z), one}}; // x^2 - 3
  
  // NOTE: This relies on ctbignum providing sqrt() -> std::optional
  auto r = roots(p);

  std::vector<GF> expected{GF(4_Z), GF(9_Z)};
  if (!check_roots(r, expected))
  {
    std::println("Quadratic GF(13) Failed. Found {} roots.", r.size());
    for(auto root : r) std::println("  Found: {}", root.value);
    return false;
  }
  std::println("Quadratic GF(13) Passed");
  return true;
}

bool test_cubic()
{
  using GF = decltype(Zq(7_Z));
  constexpr GF one(1_Z);
  constexpr GF minus_one(6_Z);

  polynomial_nttp<GF, 3> p{{minus_one, GF(0_Z), GF(0_Z), one}}; // x^3 - 1
  
  // Now this should automatically dispatch to roots_berlekamp!
  auto r = roots(p);

  std::vector<GF> expected{one, GF(2_Z), GF(4_Z)};
  if (!check_roots(r, expected))
  {
    std::println("Cubic GF(7) Failed. Found {} roots.", r.size());
    for(auto root : r) std::println("  Found: {}", root.value);
    return false;
  }
  std::println("Cubic GF(7) Passed");
  return true;
}

int main()
{
  bool result = true;
  result &= test_linear();
  result &= test_quadratic();
  result &= test_cubic();

  return result ? 0 : 1;
}
