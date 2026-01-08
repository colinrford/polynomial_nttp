/*
 *  roots_brute_force.cpp
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Brute force root finding for polynomial_nttp
 */
import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::polynomial::univariate;
using namespace lam::polynomial::univariate::roots;
using namespace lam::cbn;
using namespace lam::cbn::literals;

// GF(2) Test: x^2 + x = 0 (Roots 0, 1)
void test_gf2()
{
  constexpr auto mod = 2_Z;
  using GF = decltype(Zq(mod));

  constexpr GF zero = GF(0_Z);
  constexpr GF one = GF(1_Z);

  polynomial_nttp<GF, 2> p{{zero, one, one}}; // x^2 + x

  std::vector<GF> range = {zero, one};
  auto r = roots_brute_force(p, range);

  if (r.size() != 2)
  {
    std::println("GF(2) x^2+x failed: size {}", r.size());
    std::exit(1);
  }
  std::println("GF(2) x^2+x passed");
}

// GF(3) Test: x^3 - x = 0 (Roots 0, 1, 2)
void test_gf3()
{
  constexpr auto mod = 3_Z;
  using GF = decltype(Zq(mod));

  constexpr GF zero = GF(0_Z);
  constexpr GF one = GF(1_Z);
  constexpr GF two = GF(2_Z);
  // x^3 - x = x^3 + 2x
  polynomial_nttp<GF, 3> p{{zero, two, zero, one}};

  std::vector<GF> range = {zero, one, two};
  auto r = roots_brute_force(p, range);

  if (r.size() != 3)
  {
    std::println("GF(3) x^3-x failed: size {}", r.size());
    std::exit(1);
  }
  std::println("GF(3) x^3-x passed");
}

// GF(5) Test: x^2 + 1 = 0 (Roots 2, 3)
void test_gf5()
{
  constexpr auto mod = 5_Z;
  using GF = decltype(Zq(mod));

  constexpr GF zero = GF(0_Z);
  constexpr GF one = GF(1_Z);
  constexpr GF two = GF(2_Z);
  constexpr GF three = GF(3_Z);
  constexpr GF four = GF(4_Z);

  polynomial_nttp<GF, 2> p{{one, zero, one}}; // x^2 + 1

  std::vector<GF> range = {zero, one, two, three, four};
  auto r = roots_brute_force(p, range);

  if (r.size() != 2)
  {
    std::println("GF(5) x^2+1 failed: size {}", r.size());
    std::exit(1);
  }
  // Check specific roots
  bool has_2 = false, has_3 = false;
  for (auto root : r)
  {
    if (root.value == two)
      has_2 = true;
    if (root.value == three)
      has_3 = true;
  }
  if (has_2 && has_3)
    std::println("GF(5) x^2+1 passed");
  else
  {
    std::println("GF(5) x^2+1 failed values");
    std::exit(1);
  }
}

int main()
{
  test_gf2();
  test_gf3();
  test_gf5();
  std::println("All ctbignum brute force tests passed");
  return 0;
}
