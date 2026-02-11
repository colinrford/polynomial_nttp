/*
 *  ntt_constexpr_test.cpp
 *  Verifies that NTT implementation works in purely compile-time contexts.
 */

import std;
import lam.polynomial_nttp;
import lam.ctbignum;

using namespace lam::cbn::literals;

// Prime: 29 * 2^57 + 1
constexpr auto ntt_prime = 4179340454199820289_Z;
using Field = decltype(lam::cbn::Zq(ntt_prime));

namespace lam::polynomial::univariate
{
template<typename T, T... Modulus>
struct finite_field_traits<lam::cbn::ZqElement<T, Modulus...>>
{
  static constexpr bool is_finite_field = true;
  static constexpr T modulus = []() {
    constexpr T mods[] = {Modulus...};
    return mods[0];
  }();
};
} // namespace lam::polynomial::univariate

using namespace lam::polynomial::univariate::ntt;

// Test Function executed at compile time
consteval bool verify_ntt_compile_time()
{
  // 1. Array-based NTT (Stack/Static)
  std::array<Field, 4> data = {Field(1), Field(2), Field(3), Field(4)};
  auto original = data;

  // Forward Transform
  ntt_transform(data, false);

  // Check it did something different
  bool changed = false;
  for (std::size_t i = 0; i < 4; ++i)
    if (data[i] != original[i])
      changed = true;
  if (!changed)
    return false;

  // Inverse Transform
  ntt_transform(data, true);

  // Verify Round Trip
  for (std::size_t i = 0; i < 4; ++i)
  {
    // Strict equality for field
    if (data[i] != original[i])
      return false;
  }

  return true;
}

// Compile-time assertion
static_assert(verify_ntt_compile_time(), "NTT must work at compile time!");

int main()
{
  // Runtime check that confirms the binary was compiled successfully
  std::println("Compile-time NTT verification passed (static_assert).");
  return 0;
}
