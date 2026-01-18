/*
 *  polynomial_nttp-univariate-algebra.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

module;

#ifdef LAM_USE_TBB
#include <tbb/parallel_for.h>
#endif

export module lam.polynomial_nttp:univariate.algebra;

import std;
import lam.concepts;
import :univariate.structure;
import :univariate.fft;
import :univariate.acceleration;
import :config;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace lam::polynomial::univariate
{

const auto indexing_set = [](auto n) { return stdv::iota(static_cast<decltype(n)>(0), n); };

// enable syntax for adding polynomials in R[X]
export template<ring_element_c_weak R = double, std::size_t M, std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, M>& p, const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);
  polynomial_nttp<R, max_degree> p_plus_q{};

  // Determine if we can use acceleration
  constexpr bool is_accelerated_type = std::is_same_v<R, double> || std::is_same_v<R, float> ||
                                       std::is_same_v<R, std::complex<double>> ||
                                       std::is_same_v<R, std::complex<float>>;
  // Benchmarking shows vDSP is faster even at N=16 (21ns vs 110ns Naive)
  constexpr std::size_t threshold = 16;

  if consteval
  {
    stdr::transform(indexing_set(max_degree + 1), p_plus_q.coefficients.begin(), [&](auto i) { return p[i] + q[i]; });
  }
  else
  {
    bool accelerated = false;
    if constexpr (is_accelerated_type && max_degree >= threshold)
    {
      auto perform_acceleration = [&](const auto& larger_poly, const auto& smaller_poly) {
        stdr::copy(larger_poly.coefficients, p_plus_q.coefficients.begin());
        std::size_t min_count = std::min(M, N) + 1;

#ifdef LAM_USE_ACCELERATE
        if constexpr (std::is_same_v<R, double>)
        {
          acceleration::vDSP_vaddD(smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1,
                                   p_plus_q.coefficients.data(), 1, min_count);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, float>)
        {
          acceleration::vDSP_vadd(smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1,
                                  p_plus_q.coefficients.data(), 1, min_count);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<double>>)
        {
          acceleration::vDSP_vaddD(reinterpret_cast<const double*>(smaller_poly.coefficients.data()), 1,
                                   reinterpret_cast<double*>(p_plus_q.coefficients.data()), 1,
                                   reinterpret_cast<double*>(p_plus_q.coefficients.data()), 1, min_count * 2);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<float>>)
        {
          acceleration::vDSP_vadd(reinterpret_cast<const float*>(smaller_poly.coefficients.data()), 1,
                                  reinterpret_cast<float*>(p_plus_q.coefficients.data()), 1,
                                  reinterpret_cast<float*>(p_plus_q.coefficients.data()), 1, min_count * 2);
          accelerated = true;
        }
#endif

#ifdef LAM_USE_BLAS
        if (!accelerated)
        {
           if constexpr (std::is_same_v<R, double>)
           {
             acceleration::cblas_daxpy(min_count, 1.0, smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1);
             accelerated = true;
           }
           else if constexpr (std::is_same_v<R, float>)
           {
             acceleration::cblas_saxpy(min_count, 1.0f, smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1);
             accelerated = true;
           }
           else if constexpr (std::is_same_v<R, std::complex<double>>)
           {
             double one[] = {1.0, 0.0};
             acceleration::cblas_zaxpy(min_count, one, smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1);
             accelerated = true;
           }
           else if constexpr (std::is_same_v<R, std::complex<float>>)
           {
             float one[] = {1.0f, 0.0f};
             acceleration::cblas_caxpy(min_count, one, smaller_poly.coefficients.data(), 1, p_plus_q.coefficients.data(), 1);
             accelerated = true;
           }
        }
#endif
      };

      if constexpr (M >= N)
        perform_acceleration(p, q);
      else
        perform_acceleration(q, p);
    }

    if (!accelerated)
    {
      constexpr std::size_t tbb_threshold = 16000;
      constexpr std::size_t jthread_threshold = 50000;
      std::size_t N_size = max_degree + 1;

      if constexpr (lam::polynomial::config::use_tbb)
      {
#ifdef LAM_USE_TBB
        if (N_size >= tbb_threshold)
        {
          tbb::parallel_for(tbb::blocked_range<std::size_t>(0, N_size), [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); ++i)
            {
              p_plus_q.coefficients[i] = p[i] + q[i];
            }
          });
          accelerated = true; // effectively handled
        }
#endif
      }

      if (!accelerated && !lam::polynomial::config::use_tbb && N_size >= jthread_threshold)
      {
        std::size_t num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
          num_threads = 2;

        std::size_t chunk_size = N_size / num_threads;
        std::vector<std::jthread> threads;
        threads.reserve(num_threads);

        auto worker = [&](std::size_t start, std::size_t end) {
          for (std::size_t i = start; i < end; ++i)
          {
            p_plus_q.coefficients[i] = p[i] + q[i];
          }
        };

        for (std::size_t t = 0; t < num_threads; ++t)
        {
          std::size_t start = t * chunk_size;
          std::size_t end = (t == num_threads - 1) ? N_size : start + chunk_size;
          threads.emplace_back(worker, start, end);
        }
        accelerated = true;
      }

      if (!accelerated)
      {
        stdr::transform(indexing_set(max_degree + 1), p_plus_q.coefficients.begin(),
                        [&](auto i) { return p[i] + q[i]; });
      }
    }
  }
  return p_plus_q;
}

// add a constant on the left
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator+(const R& r, const polynomial_nttp<R, N>& p) noexcept
{
  polynomial_nttp<R, N> r_plus_p = p;
  r_plus_p.coefficients[0] = r + r_plus_p.coefficients[0];
  return r_plus_p;
}

// add a constant on the right
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator+(const polynomial_nttp<R, N>& p, const R& r) noexcept
{
  polynomial_nttp<R, N> p_plus_r = p;
  p_plus_r.coefficients[0] = p_plus_r.coefficients[0] + r;
  return p_plus_r;
}

// enable syntax for subtracting polynomials in R[X]
export template<ring_element_c_weak R = double, std::size_t M, std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, M>& p, const polynomial_nttp<R, N>& q) noexcept
{
  constexpr auto max_degree = stdr::max(M, N);
  polynomial_nttp<R, max_degree> p_minus_q{};

  constexpr bool is_accelerated_type = std::is_same_v<R, double> || std::is_same_v<R, float> ||
                                       std::is_same_v<R, std::complex<double>> ||
                                       std::is_same_v<R, std::complex<float>>;
  // Benchmarking shows vDSP is faster even at N=16 (21ns vs 110ns Naive)
  constexpr std::size_t threshold = 16;

  if consteval
  {
    stdr::transform(indexing_set(max_degree + 1), p_minus_q.coefficients.begin(), [&](auto i) { return p[i] - q[i]; });
  }
  else
  {
    bool accelerated = false;
    if constexpr (is_accelerated_type && max_degree >= threshold)
    {
#ifdef LAM_USE_ACCELERATE
      // C = A - B
      // vDSP_vsubD(B, 1, A, 1, C, 1, N);  Note: vsub matches C = A - B in docs usually, but confirm param order
      // vDSP_vsubD(B, IB, A, IA, C, IC, N) -> C = A - B
      // Wait, vDSP_vsubD says "Vector Subtract Double".
      // Documentation says: C[i] = A[i] - B[i] or it's B, A...
      // Apple docs: vDSP_vsub(B, 1, A, 1, C, 1, N) -> C = A - B.
      // So FIRST arg is SUBTRAHEND (what overrides), SECOND is MINUEND.
      // We want p - q. So A=p, B=q.
      // Call: vDSP_vsub(q, 1, p, 1, res, 1, n).

      std::size_t min_count = std::min(M, N) + 1;
      std::size_t max_count = max_degree + 1;

      // Strategy:
      // 1. Copy p to result.
      // 2. Subtract q from result (prefix).
      // 3. If q is longer than p, we need to subtract the 'tail' of q from 0.
      //    This loop handles it? p[i] - q[i]. p[i] is 0. So -q[i].
      //    Ranges fallback handles this naturally.
      //    Acceleration needs care.

      // Simplified acceleration for M == N (most common case for huge poly ops usually? or not?)
      // Let's implement full logic.

      if (M >= N)
      {
        // Result size M. Copy p fully.
        stdr::copy(p.coefficients, p_minus_q.coefficients.begin());
        // Subtract q from first N+1 elements.
        // C = A - B -> res = res - q -> res = p - q.
        // A = res (p), B = q.
        // vDSP_vsub(q, ..., res, ..., res, ...)
        if constexpr (std::is_same_v<R, double>)
        {
          acceleration::vDSP_vsubD(q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1,
                                   p_minus_q.coefficients.data(), 1, min_count);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, float>)
        {
          acceleration::vDSP_vsub(q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1,
                                  p_minus_q.coefficients.data(), 1, min_count);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<double>>)
        {
          acceleration::vDSP_vsubD(reinterpret_cast<const double*>(q.coefficients.data()), 1,
                                   reinterpret_cast<double*>(p_minus_q.coefficients.data()), 1,
                                   reinterpret_cast<double*>(p_minus_q.coefficients.data()), 1, min_count * 2);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<float>>)
        {
          acceleration::vDSP_vsub(reinterpret_cast<const float*>(q.coefficients.data()), 1,
                                  reinterpret_cast<float*>(p_minus_q.coefficients.data()), 1,
                                  reinterpret_cast<float*>(p_minus_q.coefficients.data()), 1, min_count * 2);
          accelerated = true;
        }
      }
      else
      {
        // q is larger (degree N > M).
        // Result size N.
        // Copy -q to result?
        // Or Copy q then negate? vDSP_vneg?
        // Implementation:
        // res[0..M] = p[0..M] - q[0..M]
        // res[M+1..N] = 0 - q[M+1..N]

        // Step 1: Initialize res with -q.
        // Ranges negate copy?
        // Too complex for simple acceleration block?
        // Let's just fall back for N > M to ensure correctness or use ranges.
        // Performance gain is main goal, N > M might be rare for subtraction in some contexts, but let's handle M==N
        // efficiently.
        accelerated = false;
      }
#endif

// Support BLAS daxpy: Y = alpha*X + Y.
// p - q -> Y=p, X=q, alpha=-1.
#ifdef LAM_USE_BLAS
      if (!accelerated && M >= N)
      {
        stdr::copy(p.coefficients, p_minus_q.coefficients.begin());
        std::size_t min_count = std::min(M, N) + 1;

        if constexpr (std::is_same_v<R, double>)
        {
          acceleration::cblas_daxpy(min_count, -1.0, q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, float>)
        {
          acceleration::cblas_saxpy(min_count, -1.0f, q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<double>>)
        {
          double minus_one[] = {-1.0, 0.0};
          acceleration::cblas_zaxpy(min_count, minus_one, q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1);
          accelerated = true;
        }
        else if constexpr (std::is_same_v<R, std::complex<float>>)
        {
          float minus_one[] = {-1.0f, 0.0f};
          acceleration::cblas_caxpy(min_count, minus_one, q.coefficients.data(), 1, p_minus_q.coefficients.data(), 1);
          accelerated = true;
        }
      }
#endif
    }

    if (!accelerated)
    {
      constexpr std::size_t tbb_threshold = 16000;
      constexpr std::size_t jthread_threshold = 50000;
      std::size_t N_size = max_degree + 1;

      if constexpr (lam::polynomial::config::use_tbb)
      {
#ifdef LAM_USE_TBB
        if (N_size >= tbb_threshold)
        {
          tbb::parallel_for(tbb::blocked_range<std::size_t>(0, N_size), [&](const tbb::blocked_range<std::size_t>& r) {
            for (std::size_t i = r.begin(); i != r.end(); ++i)
            {
              p_minus_q.coefficients[i] = p[i] - q[i];
            }
          });
          accelerated = true; // effectively handled
        }
#endif
      }

      if (!accelerated && !lam::polynomial::config::use_tbb && N_size >= jthread_threshold)
      {
        std::size_t num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
          num_threads = 2;

        std::size_t chunk_size = N_size / num_threads;
        std::vector<std::jthread> threads;
        threads.reserve(num_threads);

        auto worker = [&](std::size_t start, std::size_t end) {
          for (std::size_t i = start; i < end; ++i)
          {
            p_minus_q.coefficients[i] = p[i] - q[i];
          }
        };

        for (std::size_t t = 0; t < num_threads; ++t)
        {
          std::size_t start = t * chunk_size;
          std::size_t end = (t == num_threads - 1) ? N_size : start + chunk_size;
          threads.emplace_back(worker, start, end);
        }
        accelerated = true;
      }

      if (!accelerated)
      {
        stdr::transform(indexing_set(max_degree + 1), p_minus_q.coefficients.begin(),
                        [&](auto i) { return p[i] - q[i]; });
      }
    }
  }
  return p_minus_q;
}

// subtract a constant on the left
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator-(const R& r, const polynomial_nttp<R, N>& p) noexcept
{
  polynomial_nttp<R, N> r_minus_p = -p;
  r_minus_p.coefficients[0] += r; // R needs +=
  return r_minus_p;
}

// subtract a constant on the right
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator-(const polynomial_nttp<R, N>& p, const R& r) noexcept
{
  polynomial_nttp<R, N> p_minus_r = p;
  p_minus_r.coefficients[0] = p_minus_r.coefficients[0] - r;
  return p_minus_r;
}

// enable syntax for multiplying polynomials in R[X]
export template<ring_element_c_weak R = double, std::size_t M, std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, M>& p, const polynomial_nttp<R, N>& q) noexcept
{
  // Optimization: Use FFT for large degrees (N >= 32)
  // Only for double or complex<double> where we have optimized FFT
  constexpr bool is_compatible_type = std::is_same_v<R, double> || std::is_same_v<R, std::complex<double>>;
  constexpr std::size_t fft_threshold = 32;

  if constexpr (is_compatible_type && (M + N) >= fft_threshold)
  {
    if consteval
    {
      // Compile-time FFT
      std::size_t n_fft = std::bit_ceil(M + N + 1);

      std::vector<std::complex<double>> lhs_c(n_fft, {0.0, 0.0});
      std::vector<std::complex<double>> rhs_c(n_fft, {0.0, 0.0});

      for (std::size_t i = 0; i <= M; ++i)
      {
        if constexpr (std::is_same_v<R, std::complex<double>>)
          lhs_c[i] = p.coefficients[i];
        else
          lhs_c[i] = {static_cast<double>(p.coefficients[i]), 0.0};
      }
      for (std::size_t i = 0; i <= N; ++i)
      {
        if constexpr (std::is_same_v<R, std::complex<double>>)
          rhs_c[i] = q.coefficients[i];
        else
          rhs_c[i] = {static_cast<double>(q.coefficients[i]), 0.0};
      }

      fft::fft_constexpr_impl(lhs_c, false);
      fft::fft_constexpr_impl(rhs_c, false);

      for (std::size_t i = 0; i < n_fft; ++i)
        lhs_c[i] *= rhs_c[i];

      fft::fft_constexpr_impl(lhs_c, true);

      polynomial_nttp<R, M + N> p_times_q{};
      for (std::size_t i = 0; i <= M + N; ++i)
      {
        p_times_q.coefficients[i] = static_cast<R>(lhs_c[i].real());
      }
      return p_times_q;
    }
    else
    { // Runtime: Use FFT
      // 1. Calculate optimal FFT size (Power of 2)
      std::size_t fft_size = std::bit_ceil(M + N + 1);
      // 2. Allocate padded vectors directly
      std::vector<std::complex<double>> p_vec(fft_size, {0.0, 0.0});
      std::vector<std::complex<double>> q_vec(fft_size, {0.0, 0.0});
      // 3. Fill data
      for (std::size_t i = 0; i <= M; ++i)
      {
        if constexpr (std::is_same_v<R, std::complex<double>>)
          p_vec[i] = p[i];
        else
          p_vec[i] = {static_cast<double>(p[i]), 0.0};
      }
      for (std::size_t i = 0; i <= N; ++i)
      {
        if constexpr (std::is_same_v<R, std::complex<double>>)
          q_vec[i] = q[i];
        else
          q_vec[i] = {static_cast<double>(q[i]), 0.0};
      }
      // 4. FFT In-Place (Span)
      fft::fft_transform(p_vec, false);
      fft::fft_transform(q_vec, false);
      // 5. Convolution
      for (std::size_t i = 0; i < fft_size; ++i)
        p_vec[i] *= q_vec[i];
      // 6. Inverse FFT
      fft::fft_transform(p_vec, true);
      // 7. Extract Result
      polynomial_nttp<R, M + N> result{};
      for (std::size_t i = 0; i <= M + N; ++i)
      {
        if constexpr (std::is_same_v<R, double>)
          result.coefficients[i] = p_vec[i].real();
        else
          result.coefficients[i] = p_vec[i];
      }
      return result;
    }
  }
  else
  {
    polynomial_nttp<R, M + N> p_times_q{};
    for (auto&& k : indexing_set(M + N + 1))
      for (auto&& i : indexing_set(k + 1))
        p_times_q.coefficients.at(k) += p[i] * q[k - i];
    return p_times_q;
  }
}

// multiply a constant on the left
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator*(const R& r, const polynomial_nttp<R, N>& p) noexcept
{
  polynomial_nttp<R, N> r_times_p{};
  stdr::transform(stdr::cbegin(p), stdr::cend(p), stdr::begin(r_times_p), [&](auto&& p_i) { return r * p_i; });
  return r_times_p;
}

// multiply a constant on the right
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator*(const polynomial_nttp<R, N>& p, const R& r) noexcept
{
  polynomial_nttp<R, N> p_times_r{};
  stdr::transform(stdr::cbegin(p), stdr::cend(p), stdr::begin(p_times_r), [&](auto&& p_i) { return p_i * r; });
  return p_times_r;
}

// divide a constant on the right
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto operator/(const polynomial_nttp<R, N>& p, const R& r) noexcept
{
  polynomial_nttp<R, N> p_quotient_r{};
  stdr::transform(stdr::cbegin(p), stdr::cend(p), stdr::begin(p_quotient_r), [&](auto&& p_i) { return p_i / r; });
  return p_quotient_r;
}

} // end namespace lam::polynomial::univariate

namespace lam::polynomial::univariate::algebra
{

const auto indexing_set = [](auto n) { return stdv::iota(static_cast<decltype(n)>(0), n); };

const auto indexing_set_from_to = [](auto m, auto n) { return stdv::iota(static_cast<decltype(n)>(m), n); };

/*
 *  the norm of a polynomial is its degree
 *  a `polynomial_nttp` can be of degree `N > 0` with all zero coefficients
 */
export template<ring_element_c_weak R = double, std::size_t N>
constexpr auto norm([[maybe_unused]] const polynomial_nttp<R, N>& p) noexcept
{
  return N;
}
/* returns the leading coefficient, even if it is zero! */
export template<ring_element_c_weak R = double, std::size_t N>
constexpr R leading(const polynomial_nttp<R, N>& p) noexcept
{
  return p.coefficients[N];
}
/* returns the leading nonzero coefficient, UNFINISHED */
/*
//export
template<ring_element_c_weak R = double,
         std::size_t N,
         polynomial_nttp<R, N> p_of_x>
constexpr R leading()
noexcept
{
  if constexpr (p_of_x.coefficients[N] != static_cast<R>(0)) // placeholder
  {
    return leading(p_of_x);
  } else
  {
    return p_of_x.coefficients[N];
  }
} */
/* returns a new monomial of degree `N` */
export template<ring_element_c_weak R = double, std::size_t N>
constexpr polynomial_nttp<R, N> make_monomial() noexcept
{
  polynomial_nttp<R, N> new_monomial{};
  new_monomial.coefficients[N] = static_cast<R>(1);
  return new_monomial;
}

/*
 *  division_prototype()
 *    takes several NTTPs:
 *      `std::size_t M`                 – the degree of the first polynomial
 *      `polynomial_nttp<R, M> a_of_x`  – polynomial a(x) of degree M in R[X]
 *      `std::size_t N`                 – the degree of the second polynomial
 *      `polynomial_nttp<R, N> b_of_x`  – polynomial b(x) of degree N in R[X]
 *
 *    returns a `std::pair` composed of the quotient `q` and the remainder `r`
 *      the exact return type is determined inside the body of the function
 *      in particular, it is determined within the immediately invoked lambda
 *      `sizes_and_oversized_arrays_q_and_r`, a `std::pair` of `std::pair`s
 *
 *    the quotient `q` is identically `0`, and
 *    the remainder `r` is identically `a_of_x` if
 *      - `b_of_x` is identically `0`, or
 *      - `N > M`, i.e., `norm(b_of_x) > norm(a_of_x)`
\*                                                                          */
export template<field_element_c_weak R = double, std::size_t M, polynomial_nttp<R, M> a_of_x, std::size_t N,
                polynomial_nttp<R, N> b_of_x>
constexpr auto division_prototype()
{
  using divisor_info = std::tuple<bool, std::size_t, R>;

  constexpr auto b_info = [&]() -> divisor_info {
    auto b_reversed = stdv::reverse(b_of_x);
    auto it = stdr::find_if(b_reversed, [&](auto&& coeff) { return !is_approx_equal(R(0), coeff); });

    if (it != b_reversed.end())
    {
      std::size_t dist = std::distance(b_reversed.begin(), it);
      return {true, N - dist, *it};
    }
    return {false, 0, R(0)};
  }();

  constexpr bool b_is_not_zero = std::get<0>(b_info);
  constexpr std::size_t degree_of_b = std::get<1>(b_info);
  constexpr R leading_coeff_of_b = std::get<2>(b_info);

  constexpr auto sizes_and_oversized_arrays_q_and_r = [&]() {
    std::array<R, M + 1> oversized_quotient{};
    std::array<R, M + 1> oversized_remainder{};
    std::vector<R> quotient{};
    std::vector<R> remainder{};
    remainder.reserve(a_of_x.coefficients.size()); // have not compared w/o
    stdr::copy(a_of_x.coefficients.begin(), a_of_x.coefficients.end(), std::back_inserter(remainder));
    std::size_t quotient_size = 0;
    if (N <= M and b_is_not_zero)
    {
      auto deg_b = degree_of_b;
      auto lead_coeff_b = leading_coeff_of_b;
      while (remainder.size() - 1 >= deg_b)
      {
        auto s = remainder.back() / lead_coeff_b;
        quotient.push_back(s);
        std::size_t i_offset = remainder.size() - deg_b - 1;
        auto remainder_span = std::span(remainder).subspan(i_offset, deg_b + 1);
        for (auto&& [r_i, b_i] : stdv::zip(remainder_span, b_of_x.coefficients))
          r_i = r_i - (s * b_i);
        /*for (std::size_t i = 0; i <= deg_b; ++i)
          remainder[i + i_offset] = remainder[i + i_offset]
                                    - (s * b_of_x.coefficients.at(i));*/
        remainder.pop_back();
      }
      if (is_approx_equal(R(0), remainder.back()))
        remainder.pop_back();
      if (remainder.empty())
        remainder.push_back(R(0));
      std::reverse(quotient.begin(), quotient.end());
      quotient_size = quotient.size() - 1;
    }
    else // will dereference a null without something like this
    {
      quotient.push_back(R(0));
    }
    for (std::size_t i = 0; i <= quotient_size; ++i)
      oversized_quotient.at(i) = quotient[i];
    std::size_t remainder_size = remainder.size() - 1;
    for (std::size_t i = 0; i <= remainder_size; ++i)
      oversized_remainder.at(i) = remainder[i];
    return std::make_pair(std::make_pair(quotient_size, remainder_size),
                          std::make_pair(oversized_quotient, oversized_remainder));
  }();

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.first> q{};
  stdr::copy(sizes_and_oversized_arrays_q_and_r.second.first.begin(),
             sizes_and_oversized_arrays_q_and_r.second.first.begin() + norm(q) + 1, q.coefficients.begin());

  polynomial_nttp<R, sizes_and_oversized_arrays_q_and_r.first.second> r{};
  stdr::copy(sizes_and_oversized_arrays_q_and_r.second.second.begin(),
             sizes_and_oversized_arrays_q_and_r.second.second.begin() + norm(r) + 1, r.coefficients.begin());
  return std::make_pair(q, r);
}

/*  R needs to model a field and static_cast<R>(i) needs to make sense
 *
 *  return type is `auto` because:
 *    the derivative of a constant is 0, and as a polynomial is of degree 0
 *    the derivative of a degree `N > 0` polynomial is of degree `N - 1`
 *    `auto` and `if constexpr` are how we choose to avoid `0 - 1` underflow
 */
export template<field_element_c_weak R = double, std::size_t N>
constexpr auto derivative(const polynomial_nttp<R, N>& p) noexcept
{

  if constexpr (N == 0)
    return polynomial_nttp<R, 0>{};
  else
  {
    polynomial_nttp<R, N - 1> ddxp{};
    auto ids = indexing_set(N);
    stdr::transform(stdr::cbegin(ids), stdr::cend(ids), stdr::begin(ddxp),
                    [&](auto&& index) { return static_cast<R>(index + 1) * p[index + 1]; });
    return ddxp;
  }
}

/*
 *  this function returns the antiderivative of p for which a_0 = 0
 *  R needs to model a field and static_cast<R>(i) needs to make sense
 */
export template<field_element_c_weak R = double, std::size_t N>
constexpr polynomial_nttp<R, N + 1> antiderivative(const polynomial_nttp<R, N>& p) noexcept
{
  return [&]() {
    polynomial_nttp<R, N + 1> antiderivative_of_p{};
    auto ids = indexing_set_from_to(1, N + 1 + 1);
    for (auto&& i : ids)
      antiderivative_of_p.coefficients.at(i) = 1 / static_cast<R>(i) * p[i - 1];
    return antiderivative_of_p;
  }();
}

} // end namespace lam::polynomial::univariate::algebra

namespace lam
{
export using polynomial::univariate::algebra::norm;
export using polynomial::univariate::algebra::leading;
export using polynomial::univariate::algebra::make_monomial;
export using polynomial::univariate::algebra::division_prototype;
export using polynomial::univariate::algebra::derivative;
export using polynomial::univariate::algebra::antiderivative;
} // end namespace lam
