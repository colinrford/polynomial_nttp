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

export module lam.polynomial_nttp:univariate.ntt;

import std;
import :univariate.structure;
import :config;

namespace lam::polynomial::univariate::ntt
{

// =============================================================================
// Modular Arithmetic Helpers (Constexpr)
// =============================================================================

export 
constexpr std::uint64_t reduce_solinas(unsigned __int128 x)
{
  constexpr std::uint64_t P = 0xFFFFFFFF00000001ULL;

  // Split 128-bit x into four 32-bit chunks: a3, a2, a1, a0
  std::uint64_t lo = static_cast<std::uint64_t>(x);
  std::uint64_t hi = static_cast<std::uint64_t>(x >> 64);

  std::uint64_t a0 = lo & 0xFFFFFFFF;
  std::uint64_t a1 = lo >> 32;
  std::uint64_t a2 = hi & 0xFFFFFFFF;
  std::uint64_t a3 = hi >> 32;

  // Formula: Result = (a1 + a2) * 2^32 + a0 - a2 - a3
  // S = (a1 + a2) << 32 + a0
  // M = a2 + a3
  // res = S - M

  unsigned __int128 S = (static_cast<unsigned __int128>(a1 + a2) << 32) + a0;
  unsigned __int128 M = a2 + a3;

  // Handle subtraction (modulo P)
  // If S < M, we need to wrap around. Since we want mod P, we add P.
  // However, P = 2^64 - 2^32 + 1.
  // It is simpler to add a multiple of P large enough, or just use signed logic conceptually.
  // Since M is small (~66 bits max? No, 32+32=33 bits max), adding P once is enough.

  if (S < M)
  {
    S += P; // Borrow from the field modulus
  }
  S -= M;

  // Final reduction steps
  // S is now approx 65 bits. We reduce until it fits in [0, P-1].
  while (S >= P)
  {
    S -= P;
  }

  return static_cast<std::uint64_t>(S);
}

// Template-dispatched modular multiplication
// When M is known at compile-time, the compiler generates optimized Barrett reduction
export 
template<std::uint64_t M>
constexpr std::uint64_t mul_mod(std::uint64_t a, std::uint64_t b)
{
  if constexpr (M == 0)
  {
    // Handle runtime or non-field case gracefully (should not happen in NTT loop)
    return 0;
  }
  else if constexpr (M == 0xFFFFFFFF00000001ULL)
  {
    return reduce_solinas(static_cast<unsigned __int128>(a) * b);
  }
  else
  {
    return static_cast<std::uint64_t>((static_cast<unsigned __int128>(a) * b) % M);
  }
}

// Portable modular multiplication: (a * b) % m
// Uses tiered optimizations based on modulus properties
export 
constexpr auto mul_mod(std::unsigned_integral auto a, std::unsigned_integral auto b,
                              std::unsigned_integral auto m)
{
  using T = std::common_type_t<decltype(a), decltype(b), decltype(m)>;
  a %= m;
  b %= m;

  // If m is actually a compile-time constant (passed via template in ntt_transform)
  // this is effectively specialization.

  // Use specialized Solinas path if applicable
  if constexpr (std::is_same_v<T, std::uint64_t>)
  {
    if (m == 0xFFFFFFFF00000001ULL)
    {
      return static_cast<T>(reduce_solinas(static_cast<unsigned __int128>(a) * b));
    }
  }

  // Fallback: Use hardware 128-bit multiply + modulo
  if constexpr (config::has_int128)
  {
    if constexpr (sizeof(T) <= 8)
    {
      unsigned __int128 res = static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
      return static_cast<T>(res % m);
    }
  }

  // Fallback: Double-and-Add Algorithm
  if constexpr (sizeof(T) <= sizeof(std::uint32_t))
  {
    return static_cast<T>((static_cast<std::uint64_t>(a) * b) % m);
  }
  else
  {
    T res = 0;
    while (b > 0)
    {
      if (b % 2 == 1)
      {
        if (res < m - a)
          res += a;
        else
          res -= (m - a);
      }
      if (b > 1)
      {
        if (a < m - a)
          a += a;
        else
          a -= (m - a);
      }
      b /= 2;
    }
    return res;
  }
}

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

} // namespace lam::polynomial::univariate::ntt
