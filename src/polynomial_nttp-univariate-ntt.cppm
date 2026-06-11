// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: 2025-2026 Colin Ford

/*
 *  polynomial_nttp-univariate-ntt.cppm
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *    Number Theoretic Transform (NTT) implementation
 */

module;

#ifdef LAM_USE_TBB
#include <tbb/parallel_for.h>
#endif

export module lam.polynomial.nttp:univariate.ntt;

import std;
import :univariate.structure;
import :config;

namespace lam::polynomial::nttp::univariate::ntt
{

// =============================================================================
// Modular Arithmetic Helpers (Constexpr)
// =============================================================================

namespace detail
{

struct u128_parts { std::uint64_t hi; std::uint64_t lo; };

// Schoolbook 64x64 -> 128 multiply; no 128-bit type, so MSVC-portable.
constexpr u128_parts mul_wide_portable(std::uint64_t a, std::uint64_t b)
{
  const std::uint64_t a_lo = a & 0xFFFFFFFFULL, a_hi = a >> 32;
  const std::uint64_t b_lo = b & 0xFFFFFFFFULL, b_hi = b >> 32;

  const std::uint64_t ll = a_lo * b_lo;
  const std::uint64_t lh = a_lo * b_hi;
  const std::uint64_t hl = a_hi * b_lo;
  const std::uint64_t hh = a_hi * b_hi;

  const std::uint64_t cross = (ll >> 32) + (lh & 0xFFFFFFFFULL) + (hl & 0xFFFFFFFFULL);
  const std::uint64_t lo = (ll & 0xFFFFFFFFULL) | (cross << 32);
  const std::uint64_t hi = hh + (lh >> 32) + (hl >> 32) + (cross >> 32);
  return { hi, lo };
}

// Native where available, else schoolbook. The token is preprocessor-gated:
// if constexpr cannot hide __int128 from a front-end that rejects it.
constexpr u128_parts mul_wide(std::uint64_t a, std::uint64_t b)
{
#if defined(LAM_HAS_INT128)
  const unsigned __int128 p = static_cast<unsigned __int128>(a) * b;
  return { static_cast<std::uint64_t>(p >> 64), static_cast<std::uint64_t>(p) };
#else
  return mul_wide_portable(a, b);
#endif
}

// (a * b) % m for a 64-bit modulus, no 128-bit type (double-and-add).
constexpr std::uint64_t mulmod64_portable(std::uint64_t a, std::uint64_t b,
                                          std::uint64_t m)
{
  a %= m;
  b %= m;
  std::uint64_t res = 0;
  while (b > 0)
  {
    if (b & 1ULL)
      res = (res < m - a) ? res + a : res - (m - a);
    if (b > 1)
      a = (a < m - a) ? a + a : a - (m - a);
    b >>= 1;
  }
  return res;
}

} // namespace detail

// Goldilocks reduction for P = 2^64 - 2^32 + 1 on a 128-bit value (hi, lo).
// 2^64 == 2^32 - 1 (mod P), 2^96 == -1 (mod P)  =>
//   x == lo + a2*2^32 - a2 - a3 (mod P),  a2 = hi low half, a3 = hi high half.
export
constexpr std::uint64_t reduce_solinas(std::uint64_t hi, std::uint64_t lo)
{
  constexpr std::uint64_t P = 0xFFFFFFFF00000001ULL;

  const std::uint64_t a2 = hi & 0xFFFFFFFFULL;
  const std::uint64_t a3 = hi >> 32;

  const std::uint64_t term = (a2 << 32) - a2;    // a2*(2^32 - 1), < P

  std::uint64_t r = lo + term;
  if (r < lo)                  // carry out of 64 bits
    r += 0xFFFFFFFFULL;        // 2^64 == 2^32 - 1 (mod P)

  if (r < a3)
    r += P;
  r -= a3;

  while (r >= P)
    r -= P;
  return r;
}

// Back-compat overload for callers holding a native 128-bit value.
#if defined(LAM_HAS_INT128)
export
constexpr std::uint64_t reduce_solinas(unsigned __int128 x)
{
  return reduce_solinas(static_cast<std::uint64_t>(x >> 64),
                        static_cast<std::uint64_t>(x));
}
#endif

// (a * b) % m, tiered by operand width and modulus.
export
constexpr auto mul_mod(std::unsigned_integral auto a, std::unsigned_integral auto b,
                              std::unsigned_integral auto m)
{
  using T = std::common_type_t<decltype(a), decltype(b), decltype(m)>;
  a %= m;
  b %= m;

  // Goldilocks prime -> Solinas reduction.
  if constexpr (std::is_same_v<T, std::uint64_t>)
    if (m == 0xFFFFFFFF00000001ULL)
    {
      const auto [hi, lo] = detail::mul_wide(a, b);
      return static_cast<T>(reduce_solinas(hi, lo));
    }

  if constexpr (sizeof(T) <= sizeof(std::uint32_t))
    return static_cast<T>((static_cast<std::uint64_t>(a) * b) % m);
  else
  {
#if defined(LAM_HAS_INT128)
    return static_cast<T>(static_cast<unsigned __int128>(a) * b % m);
#else
    return static_cast<T>(detail::mulmod64_portable(a, b, m));
#endif
  }
}

// Compile-time M: the Goldilocks branch folds to the Solinas fast path.
export
template<std::uint64_t M>
constexpr std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b)
{
  if constexpr (M == 0)
    return 0;
  else if constexpr (M == 0xFFFFFFFF00000001ULL)
  {
    const auto [hi, lo] = detail::mul_wide(a, b);
    return reduce_solinas(hi, lo);
  }
  else
    return mul_mod(a, b, M);
}

// Where __int128 exists, prove the portable kernels match it at compile time;
// passing here is the assurance for the cl.exe build, which has no native ref.
#if defined(LAM_HAS_INT128)
namespace detail::selfcheck
{
  constexpr bool wide_matches_native(std::uint64_t a, std::uint64_t b)
  {
    const unsigned __int128 p = static_cast<unsigned __int128>(a) * b;
    const auto w = mul_wide_portable(a, b);
    return w.hi == static_cast<std::uint64_t>(p >> 64)
        && w.lo == static_cast<std::uint64_t>(p);
  }

  constexpr bool solinas_matches_native(std::uint64_t hi, std::uint64_t lo)
  {
    constexpr unsigned __int128 P = 0xFFFFFFFF00000001ULL;
    const unsigned __int128 x = (static_cast<unsigned __int128>(hi) << 64) | lo;
    return reduce_solinas(hi, lo) == static_cast<std::uint64_t>(x % P);
  }

  constexpr bool run_all()
  {
    std::uint64_t s = 0x243F6A8885A308D3ULL;
    auto next = [&s]() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; };
    for (int i = 0; i < 256; ++i)
    {
      const std::uint64_t a = next(), b = next();
      if (!wide_matches_native(a, b))              return false;
      if (!solinas_matches_native(next(), next())) return false;
    }
    if (!wide_matches_native(~0ULL, ~0ULL))     return false;
    if (!solinas_matches_native(0, 0))          return false;
    if (!solinas_matches_native(~0ULL, ~0ULL))  return false;
    if (!solinas_matches_native(0, 1ULL << 63)) return false;
    if (!solinas_matches_native(1ULL << 32, 0)) return false;
    return true;
  }

  static_assert(run_all(), "portable 128-bit kernels disagree with native __int128");
}
#endif

// Computes (base^exp) % mod
export 
constexpr auto power(std::unsigned_integral auto base, std::unsigned_integral auto exp,
                            std::unsigned_integral auto mod)
{
  using T = decltype(base);
  T res = 1;
  base %= mod;
  while (exp > 0)
  {
    if (exp % 2 == 1)
      res = mul_mod(res, base, mod);
    base = mul_mod(base, base, mod);
    exp /= 2;
  }
  return res;
}

// Computes modular inverse using Fermat's Little Theorem (for prime mod)
// inv(n) = n^(mod-2) % mod
export 
constexpr auto mod_inverse(std::unsigned_integral auto n, std::unsigned_integral auto mod)
{ return power(n, mod - 2, mod); }

// Check if a primitive n-th root of unity exists for modulus p
// Returns the root if found, 0 otherwise.
// We search for a generator 'g' such that order of g is at least n.
// Then w = g^((p-1)/n) has order n.
export 
template<typename T>
constexpr T find_nth_root_of_unity(T n, T p)
{
  if (n == 0 || p <= 1)
    return 0;

  if ((p - 1) % n != 0)
    return 0;

  T exponent = (p - 1) / n;

  // Try small bases to find a generator.
  // For the cryptographic primes we care about (Solinas, Proth),
  // small generators (like 3, 5, 7) usually exist.
  for (T g = 2; g < 1000; ++g)
  {
    // w = g^((p-1)/n) mod p
    T w = power(g, exponent, p);

    // Check if w has order n.
    // We know w^n = (g^((p-1)/n))^n = g^(p-1) = 1 mod p (Fermat).
    // Sufficiency check: w^(n/2) != 1 mod p.
    // This is necessary and sufficient if n is a power of 2 (which is true for NTT).
    if (w != 1)
    {
      if (power(w, n / 2, p) != static_cast<T>(1))
        return w;
    }
  }
  return 0;
}

// Find a primitive root generator 'g' for prime 'p'
// Returns a generator g such that order(g) = p-1.
export 
constexpr auto find_primitive_root(std::unsigned_integral auto p)
{
  using T = decltype(p);
  // Iterating small integers is standard practice for finding primitive roots
  // for primes of cryptographic interest.
  for (T g = 2; g < 1000; ++g)
  {
    // Euler's totient phi(p) = p-1 for prime p.
    // We need g^(p-1) = 1 (always true per Fermat if g not mult of p)
    // AND g^((p-1)/q) != 1 for all prime factors q of p-1.
    // Since factoring p-1 is expensive at compile time, we rely on the
    // user or higher-level logic to provide primes where we know the structure,
    // or we rely on `find_nth_root_of_unity` which effectively checks the
    // subgroup order we actually care about (N).

    // For general usage, just returning a likely candidate or checking
    // simple divisors (2) is a reasonable default for this module level.
    // Logic: if g^((p-1)/2) != 1, it is a primitive root if p = 2q+1 (safe prime).
    if (power(g, (p - 1) / 2, p) != static_cast<T>(1))
      return g;
  }
  return T(0);
}


// =============================================================================
// NTT Transform
// =============================================================================

/*
 * In-place Number Theoretic Transform
 * data: range of coefficients (field elements)
 * inverse: true for Inverse NTT
 *
 * Requirements:
 * - data.size() must be a power of 2
 * - field K must have a modulus P where data.size() divides P-1
 */
export 
template<typename Range>
constexpr void ntt_transform(Range&& data, bool inverse)
{
  using K = std::remove_cvref_t<decltype(data[0])>;
  using traits = finite_field_traits<K>;
  // Check for field at compile time
  if constexpr (!traits::is_finite_field)
    return;
  else
  {
    using val_t = std::decay_t<decltype(traits::field_order)>;
    constexpr val_t P = static_cast<val_t>(traits::field_order);
    (void)P;

    std::size_t n = std::ranges::size(data);

    if ((n & (n - 1)) != 0 || n == 0) // Check power of 2
      return;

    // 1. Find Root of Unity
    val_t root = find_nth_root_of_unity<val_t>(static_cast<val_t>(n), P);

    if (root == 0)
    {
      return;
    }
    (void)root;

    // If inverse, use w^-1
    if (inverse)
      root = mod_inverse(root, P);

    // 2. Bit-Reversal Permutation
    std::size_t j = 0;
    for (std::size_t i = 1; i < n; i++)
    {
      std::size_t bit = n >> 1;
      while (j & bit)
      {
        j ^= bit;
        bit >>= 1;
      }
      j ^= bit;
      if (i < j)
      {
        std::swap(data[i], data[j]);
      }
    }
    // 3. Butterfly Operations
    for (std::size_t len = 2; len <= n; len <<= 1)
    {
      val_t wlen = power(root, static_cast<val_t>(n / len), P);

      // Serial Butterfly Loop (simplified for verification)
      for (std::size_t i = 0; i < n; i += len)
      {
        val_t w = 1;
        for (std::size_t j = 0; j < len / 2; j++)
        {
          K u = data[i + j];
          // Traits usage!
          K v = traits::mul(data[i + j + len / 2], K(w));

          data[i + j] = traits::add(u, v);
          data[i + j + len / 2] = traits::sub(u, v);

          w = mul_mod(w, wlen, P);
        }
      }
    }

    // 4. Scaling for Inverse
    if (inverse)
    {
      val_t inv_n = mod_inverse(static_cast<val_t>(n), P);
      K inv_n_k = K(inv_n);
      for (std::size_t i = 0; i < n; ++i)
      {
        data[i] = data[i] * inv_n_k;
      }
    }
  }
}

} // namespace lam::polynomial::nttp::univariate::ntt
