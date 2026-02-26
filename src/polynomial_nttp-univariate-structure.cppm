/*
 *  polynomial_nttp-univariate-structure.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

module;


#ifdef LAM_USE_TBB
#include <tbb/parallel_for.h>
#endif

export module lam.polynomial_nttp:univariate.structure;

import std;
import lam.concepts;
import :config;
import :univariate.acceleration;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace lam::polynomial::univariate
{

template<typename R>
concept ring_element_c_weak = lam::concepts::experimental::ring_element_c_weak<R>;
template<typename R>
concept field_element_c_weak = lam::concepts::experimental::field_element_c_weak<R>;

// Use acceleration namespace for vDSP types/functions
#ifdef LAM_USE_ACCELERATE
using namespace acceleration;
#endif

// Primary template for finite field detection
// Users or interop modules specialize this for their field types
export template<typename K>
struct finite_field_traits
{
  static constexpr bool is_finite_field = false;
  static constexpr std::size_t field_order = 0;

  // Customization points for optimized arithmetic
  // By default, uses the global mul_mod
  static constexpr auto mul(const K& a, const K& b) { return a * b; }
  static constexpr auto add(const K& a, const K& b) { return a + b; }
  static constexpr auto sub(const K& a, const K& b) { return a - b; }
};

export template<typename T>
constexpr auto get_epsilon()
{
  if constexpr (std::is_floating_point_v<T>)
    return std::numeric_limits<T>::epsilon();
  else if constexpr (requires { typename T::value_type; })
    if constexpr (std::is_floating_point_v<typename T::value_type>)
      return std::numeric_limits<typename T::value_type>::epsilon();
  return 0.; // dummy
}

export template<typename T>
constexpr bool is_approx_equal(const T& a, const T& b, std::optional<double> abs_tol = std::nullopt)
{
  if constexpr (std::is_floating_point_v<T>)
  {
    auto diff = a - b;
    if (diff < static_cast<T>(0))
      diff = -diff;
    if (abs_tol)
      return diff < static_cast<T>(*abs_tol);
    return diff < (std::numeric_limits<T>::epsilon() * static_cast<T>(2));
  }
  else
  {
    if constexpr (requires { typename T::value_type; })
    {
      using value_t = typename T::value_type;
      if constexpr (std::is_floating_point_v<value_t> && requires(T r) {
                      { r.real() } -> std::convertible_to<value_t>;
                      { r.imag() } -> std::convertible_to<value_t>;
                    })
      {
        auto diff = a - b;
        auto real_diff = diff.real();
        auto imag_diff = diff.imag();
        if (real_diff < static_cast<value_t>(0))
          real_diff = -real_diff;
        if (imag_diff < static_cast<value_t>(0))
          imag_diff = -imag_diff;

        if (abs_tol)
        {
          auto tol = static_cast<value_t>(*abs_tol);
          return (real_diff < tol) && (imag_diff < tol);
        }

        auto eps = std::numeric_limits<value_t>::epsilon() * static_cast<value_t>(2);
        return (real_diff < eps) && (imag_diff < eps);
      }
    }
    return a == b;
  }
}

// Custom relative tolerance check to survive Runge's Phenomenon / scale variance.
// Error = |expected - actual| / max(|expected|, |actual|)
export 
template <std::floating_point T>
constexpr bool is_relatively_approx_equal(T expected, T actual, T tolerance = 1e-9) 
{
  if (expected == actual) 
    return true;
  T diff = expected > actual ? expected - actual : actual - expected;
  T abs_e = expected > 0 ? expected : -expected;
  T abs_a = actual > 0 ? actual : -actual;
  T max_val = abs_e > abs_a ? abs_e : abs_a;
  if (max_val == 0) 
    return diff < tolerance; // Fallback for 0 against near 0
  return (diff / max_val) < tolerance;
}

export 
template<typename T>
constexpr bool is_negligible(const T& val, std::optional<double> abs_tol = std::nullopt)
{ return is_approx_equal(val, T(0), abs_tol); }

export 
template<ring_element_c_weak R = double, std::size_t N = 0>
struct polynomial_nttp
{
  using coefficient_t = R; // std::type_identity_t<R>?
  using index_t = decltype(N);

  index_t degree = N;
  std::array<coefficient_t, N + 1> coefficients{};

  // possibly allegedly std::fill may not perform this sometimes?
  // could be it's only a thing with std::vector
  constexpr polynomial_nttp() noexcept
  { stdr::fill(stdr::begin(coefficients), stdr::end(coefficients), coefficient_t(0)); }

  constexpr polynomial_nttp(std::initializer_list<coefficient_t> li) noexcept
  { std::copy(li.begin(), li.end(), coefficients.begin()); }

  explicit constexpr polynomial_nttp(std::array<coefficient_t, N + 1> coeffs) noexcept : coefficients(coeffs) {}

  explicit constexpr polynomial_nttp(coefficient_t constant) noexcept
  {
    if constexpr (N > 0)
    {
      stdr::fill(stdr::begin(coefficients), stdr::end(coefficients), coefficient_t(0));
    }
    coefficients[0] = constant;
  }

  static constexpr polynomial_nttp zero() noexcept { return polynomial_nttp{}; }

  using begin_type = decltype(stdr::begin(coefficients));
  using end_type = decltype(stdr::end(coefficients));
  using cbegin_type = decltype(stdr::cbegin(coefficients));
  using cend_type = decltype(stdr::cend(coefficients));
  using rbegin_type = decltype(stdr::rbegin(coefficients));
  using rend_type = decltype(stdr::rend(coefficients));
  using crbegin_type = decltype(stdr::crbegin(coefficients));
  using crend_type = decltype(stdr::crend(coefficients));
  constexpr auto begin() { return std::forward<begin_type>(stdr::begin(coefficients)); }
  constexpr auto end() { return std::forward<end_type>(stdr::end(coefficients)); }
  constexpr auto begin() const { return std::forward<cbegin_type>(stdr::begin(coefficients)); }
  constexpr auto end() const { return std::forward<cend_type>(stdr::end(coefficients)); }
  constexpr auto cbegin() const { return std::forward<cbegin_type>(stdr::cbegin(coefficients)); }
  constexpr auto cend() const { return std::forward<cend_type>(stdr::cend(coefficients)); }
  constexpr auto rbegin() { return std::forward<rbegin_type>(stdr::rbegin(coefficients)); }
  constexpr auto rend() { return std::forward<rend_type>(stdr::rend(coefficients)); }
  constexpr auto crbegin() const { return std::forward<crbegin_type>(stdr::crbegin(coefficients)); }
  constexpr auto crend() const { return std::forward<crend_type>(stdr::crend(coefficients)); }

  constexpr coefficient_t operator[](index_t index) noexcept
  { return index <= N ? coefficients.at(index) : coefficient_t(0); }

  constexpr coefficient_t operator[](index_t index) const noexcept
  { return index <= N ? coefficients.at(index) : coefficient_t(0); }

  // Accelerated evaluation (Apple vDSP)
  // Uses vDSP_vsmaD for efficient Strided Multiply-Add
  constexpr coefficient_t evaluate_accelerate(coefficient_t x) const noexcept
  {
#ifdef LAM_USE_ACCELERATE
    if constexpr (std::is_same_v<coefficient_t, double>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<double, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      vDSP_Length active_size = static_cast<vDSP_Length>(N + 1);
      double pow_x = x;
      // vDSP_vsmaD: vector scalar multiply add double
      // D[i] = A[i] * S + B[i]
      // Implementation of Estrin's layer:
      // P_new[i] = P_old[2*i+1] * x + P_old[2*i]
      // A (input to multiply) : P_old[2*i+1] (odd elements) -> stride 2, offset 1
      // S (scalar)            : x
      // B (input to add)      : P_old[2*i]   (even elements) -> stride 2, offset 0
      // D (dest)              : workspace    (contiguous)    -> stride 1

      while (active_size > 1)
      {
        vDSP_Length half_size = active_size / 2;
        // Perform the layer reduction
        // A pointer = workspace + 1
        // B pointer = workspace
        // D pointer = workspace
        acceleration::vDSP_vsmaD(workspace.data() + 1, 2, &pow_x, workspace.data(), 2, workspace.data(), 1, half_size);
        // Handle odd element if size is odd (just copy it to the end of new sequence)
        if (active_size % 2 != 0)
        {
          workspace[half_size] = workspace[active_size - 1];
          active_size = half_size + 1;
        }
        else
          active_size = half_size;

        pow_x *= pow_x;
      }
      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, std::complex<double>>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<std::complex<double>, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      // Zero-Copy Stride Trick for Split Complex View
      // We point to the interleaved data with stride 2
      double* raw_ptr = reinterpret_cast<double*>(workspace.data());
      acceleration::DSP_double_split_complex split_view;
      split_view.realp = raw_ptr;     // Points to real part of 0th element
      split_view.imagp = raw_ptr + 1; // Points to imag part of 0th element

      vDSP_Length active_size = static_cast<vDSP_Length>(N + 1);
      std::complex<double> pow_x = x;

      // Scalar View: Point to the single scalar x
      // vDSP treats stride 0 as "stay on this element", effectively broadcasting it
      double* x_ptr = reinterpret_cast<double*>(&pow_x);
      acceleration::DSP_double_split_complex scalar_view;
      scalar_view.realp = x_ptr;
      scalar_view.imagp = x_ptr + 1;

      while (active_size > 1)
      {
        vDSP_Length half_size = active_size / 2;

        // vDSP_zmaD(A, IA, B, IB, C, IC, D, ID, N)
        // D = (A * B) + C
        // A: Odd elements.  Start at index 1 -> offset (1*2)=2 doubles. Stride 2.
        // B: Scalar x.      Start at x.     -> offset 0.            Stride 0.
        // C: Even elements. Start at index 0 -> offset 0.            Stride 2.
        // D: Even elements. Start at index 0 -> offset 0.            Stride 2.

        // Adjust pointers for A (Odd elements)
        acceleration::DSP_double_split_complex A_view;
        A_view.realp = split_view.realp + 2; // +1 complex element = +2 doubles
        A_view.imagp = split_view.imagp + 2;

        acceleration::vDSP_zvmaD(&A_view, 4, &scalar_view, 0, &split_view, 4, &split_view, 2, half_size);

        if (active_size % 2 != 0)
        {
          workspace[half_size] = workspace[active_size - 1];
          active_size = half_size + 1;
        }
        else
          active_size = half_size;

        pow_x *= pow_x;
      }
      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, float>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<float, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      vDSP_Length active_size = static_cast<vDSP_Length>(N + 1);
      float pow_x = x;

      while (active_size > 1)
      {
        vDSP_Length half_size = active_size / 2;
        acceleration::vDSP_vsma(workspace.data() + 1, 2, &pow_x, workspace.data(), 2, workspace.data(), 1, half_size);
        if (active_size % 2 != 0)
        {
          workspace[half_size] = workspace[active_size - 1];
          active_size = half_size + 1;
        }
        else
          active_size = half_size;

        pow_x *= pow_x;
      }
      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, std::complex<float>>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<std::complex<float>, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      float* raw_ptr = reinterpret_cast<float*>(workspace.data());
      acceleration::DSP_float_split_complex split_view;
      split_view.realp = raw_ptr;
      split_view.imagp = raw_ptr + 1;

      vDSP_Length active_size = static_cast<vDSP_Length>(N + 1);
      std::complex<float> pow_x = x;

      float* x_ptr = reinterpret_cast<float*>(&pow_x);
      acceleration::DSP_float_split_complex scalar_view;
      scalar_view.realp = x_ptr;
      scalar_view.imagp = x_ptr + 1;

      while (active_size > 1)
      {
        vDSP_Length half_size = active_size / 2;

        acceleration::DSP_float_split_complex A_view;
        A_view.realp = split_view.realp + 2;
        A_view.imagp = split_view.imagp + 2;

        acceleration::vDSP_zvma(&A_view, 4, &scalar_view, 0, &split_view, 4, &split_view, 2, half_size);

        if (active_size % 2 != 0)
        {
          workspace[half_size] = workspace[active_size - 1];
          active_size = half_size + 1;
        }
        else
          active_size = half_size;

        pow_x *= pow_x;
      }
      return workspace[0];
    }
#endif
    return coefficient_t(0); // Should be unreachable via checks
  }

  // Generic BLAS evaluation (cblas_daxpy)
  // Uses daxpy to update evens in-place: Y <- alpha*X + Y
  constexpr coefficient_t evaluate_blas(coefficient_t x) const noexcept
  {
#ifdef LAM_USE_BLAS
    if constexpr (std::is_same_v<coefficient_t, double>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<double, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      int active_size = static_cast<int>(N + 1);
      double pow_x = x;
      // Estrin's via daxpy:
      // P[2*i] = P[2*i] + x * P[2*i+1]
      // This is Y = Y + alpha * X
      // Y (dest/src) : Even elements (stride 2)
      // X (src)      : Odd elements (stride 2)
      // alpha        : x

      // Since cblas_daxpy doesn't straightforwardly support a "gather" stride for X and Y in the way
      // that reduces the array contiguousness, we maintain the "virtual" stride.
      // Layer 0: Stride 2. Pairs (0,1), (2,3), etc. Result in 0, 2, 4...
      // Layer 1: Stride 4. Pairs (0,2), (4,6)... Result in 0, 4, 8...
      int stride = 1;
      // Loop while we have at least one pair to process
      while (active_size > 1)
      {
        int half_size = active_size / 2;
        // cblas_daxpy(n, alpha, X, incX, Y, incY)
        // n: number of pairs (half_size)
        // alpha: pow_x
        // X: Odd elements. Start at workspace[stride]. Increment: 2*stride
        // Y: Even elements. Start at workspace[0].      Increment: 2*stride
        acceleration::cblas_daxpy(half_size, pow_x, workspace.data() + stride, stride * 2, workspace.data(),
                                  stride * 2);

        if (active_size % 2 != 0)
          active_size = half_size + 1;
        else
          active_size = half_size;
        stride *= 2;
        pow_x *= pow_x;
      }

      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, std::complex<double>>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      // Ensure alignment/layout is compatible (it is standard, but good to be careful)
      std::array<std::complex<double>, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      int active_size = static_cast<int>(N + 1);
      std::complex<double> pow_x = x;
      int stride = 1;

      while (active_size > 1)
      {
        int half_size = active_size / 2;
        // zaxpy takes alpha by pointer
        acceleration::cblas_zaxpy(half_size, &pow_x, workspace.data() + stride, stride * 2, workspace.data(),
                                  stride * 2);

        if (active_size % 2 != 0)
          active_size = half_size + 1;
        else
          active_size = half_size;
        stride *= 2;
        pow_x *= pow_x;
      }
      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, float>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<float, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      int active_size = static_cast<int>(N + 1);
      float pow_x = x;
      int stride = 1;

      while (active_size > 1)
      {
        int half_size = active_size / 2;
        acceleration::cblas_saxpy(half_size, pow_x, workspace.data() + stride, stride * 2, workspace.data(),
                                  stride * 2);

        if (active_size % 2 != 0)
          active_size = half_size + 1;
        else
          active_size = half_size;
        stride *= 2;
        pow_x *= pow_x;
      }

      return workspace[0];
    }
    else if constexpr (std::is_same_v<coefficient_t, std::complex<float>>)
    {
      if constexpr (N == 0)
        return coefficients[0];

      std::array<std::complex<float>, N + 1> workspace;
      std::copy(coefficients.begin(), coefficients.end(), workspace.begin());

      int active_size = static_cast<int>(N + 1);
      std::complex<float> pow_x = x;
      int stride = 1;

      while (active_size > 1)
      {
        int half_size = active_size / 2;
        acceleration::cblas_caxpy(half_size, &pow_x, workspace.data() + stride, stride * 2, workspace.data(),
                                  stride * 2);

        if (active_size % 2 != 0)
          active_size = half_size + 1;
        else
          active_size = half_size;
        stride *= 2;
        pow_x *= pow_x;
      }
      return workspace[0];
    }
#endif
    return coefficient_t(0);
  }

  // Parallel evaluation (Estrin's Scheme)
  constexpr coefficient_t evaluate_parallel(coefficient_t x) const noexcept
  {
    if constexpr (N == 0)
      return coefficients[0];
    std::array<coefficient_t, N + 1> workspace = coefficients;
    std::size_t active_size = N + 1;
    coefficient_t pow_x = x;

    while (active_size > 1)
    {
      const std::size_t half_size = active_size / 2;

      for (std::size_t i = 0; i < half_size; ++i)
        workspace[i] = workspace[2 * i] + workspace[2 * i + 1] * pow_x;
      if (active_size % 2 != 0)
      {
        workspace[half_size] = workspace[active_size - 1];
        active_size = half_size + 1;
      }
      else
        active_size = half_size;

      pow_x = pow_x * pow_x;
    }

    return workspace[0];
  }

  // Horner's method (Runtime version)
  constexpr coefficient_t evaluate_horner(coefficient_t x) const noexcept
  {
    if constexpr (N == 0)
      return coefficients[0];
    coefficient_t result = coefficients[N];
    for (index_t i = N; i > 0; --i)
      result = result * x + coefficients[i - 1];
    return result;
  }

  // Horner's method: p(x) = a₀ + x(a₁ + x(a₂ + x(...)))
  constexpr coefficient_t operator()(const coefficient_t x) const noexcept
  {
    if consteval
    {
      if constexpr (N == 0)
        return coefficients[0];
      coefficient_t result = coefficients[N];
      for (index_t i = N; i > 0; --i)
        result = result * x + coefficients[i - 1];
      return result;
    }
    else
    {
      if constexpr (lam::polynomial::config::use_accelerate)
      {
        if constexpr (std::is_same_v<coefficient_t, double> && N >= 80)
        {
#ifdef LAM_USE_ACCELERATE
          return evaluate_accelerate(x);
#endif
        }
        if constexpr (std::is_same_v<coefficient_t, std::complex<double>> && N >= 80)
        {
#ifdef LAM_USE_ACCELERATE
          return evaluate_accelerate(x);
#endif
        }
        if constexpr (std::is_same_v<coefficient_t, float> && N >= 80)
        {
#ifdef LAM_USE_ACCELERATE
          return evaluate_accelerate(x);
#endif
        }
        if constexpr (std::is_same_v<coefficient_t, std::complex<float>> && N >= 80)
        {
#ifdef LAM_USE_ACCELERATE
          return evaluate_accelerate(x);
#endif
        }
      }

      if constexpr (lam::polynomial::config::use_blas)
      {
        // Double precision support
        if constexpr (std::is_same_v<coefficient_t, double> && N >= 80)
        {
#ifdef LAM_USE_BLAS
          return evaluate_blas(x);
#endif
        }
        // Complex double support
        if constexpr (std::is_same_v<coefficient_t, std::complex<double>> && N >= 80)
        {
#ifdef LAM_USE_BLAS
          return evaluate_blas(x);
#endif
        }
        // Single precision support
        if constexpr (std::is_same_v<coefficient_t, float> && N >= 80)
        {
#ifdef LAM_USE_BLAS
          return evaluate_blas(x);
#endif
        }
        // Complex single support
        if constexpr (std::is_same_v<coefficient_t, std::complex<float>> && N >= 80)
        {
#ifdef LAM_USE_BLAS
          return evaluate_blas(x);
#endif
        }
      }

      return evaluate_horner(x);
    }
  }

  // Bulk Evaluation (One polynomial, many input points)
  constexpr void evaluate_bulk(std::span<const coefficient_t> inputs, std::span<coefficient_t> outputs) const
  {
    if (inputs.size() != outputs.size())
    {
      throw std::invalid_argument("evaluate_bulk: inputs and outputs size mismatch");
    }

    std::size_t size = inputs.size();

    // 1. TBB Path (Preferred for overhead)
    if constexpr (lam::polynomial::config::use_tbb)
    {
#ifdef LAM_USE_TBB
      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size), [&](const tbb::blocked_range<std::size_t>& r) {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
          outputs[i] = (*this)(inputs[i]);
      });
      return;
#endif
    }

    // 2. Std Threading Fallback (Runtime Only)
    // Only use if size is large enough to justify thread overhead (> 1000 items)
    if consteval
    {
      // Constexpr must be serial
      for (std::size_t i = 0; i < size; ++i)
        outputs[i] = (*this)(inputs[i]);
    }
    else
    {
      if (size > 1000)
      {
        std::size_t num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0)
          num_threads = 2; // Fallback

        std::size_t chunk_size = size / num_threads;
        std::vector<std::jthread> threads;
        threads.reserve(num_threads);

        for (std::size_t t = 0; t < num_threads; ++t)
        {
          std::size_t start = t * chunk_size;
          std::size_t end = (t == num_threads - 1) ? size : start + chunk_size;

          threads.emplace_back([this, start, end, inputs, outputs]() {
            for (std::size_t i = start; i < end; ++i)
            {
              outputs[i] = (*this)(inputs[i]);
            }
          });
        }
        // jthreads join automatically on destruction
      }
      else
      {
        // Serial Fallback for small sizes
        for (std::size_t i = 0; i < size; ++i)
          outputs[i] = (*this)(inputs[i]);
      }
    }
  }

  // raw loops may be better in constexpr context
  constexpr polynomial_nttp operator-() const noexcept
  {
    polynomial_nttp negated_polynomial{};
    stdr::transform(stdr::cbegin(coefficients), stdr::cend(coefficients), stdr::begin(negated_polynomial),
                    std::negate());
    return negated_polynomial;
  };
};

// ============================================================
// Roots Types (Moved from roots module to break circular dependency)
// ============================================================

// Root with algebraic multiplicity
export 
template<typename K>
struct root_with_multiplicity
{
  K value;
  std::size_t multiplicity{1};

  constexpr bool operator==(const root_with_multiplicity&) const = default;
};

// Fixed-capacity roots container for compile-time and runtime use
export 
template<typename K, std::size_t max_roots>
struct roots_result
{
  std::array<root_with_multiplicity<K>, max_roots> data{};
  std::size_t count{0};

  constexpr auto begin() const { return data.begin(); }
  constexpr auto end() const { return data.begin() + count; }
  constexpr auto begin() { return data.begin(); }
  constexpr auto end() { return data.begin() + count; }
  constexpr std::size_t size() const { return count; }
  constexpr bool empty() const { return count == 0; }
  constexpr const auto& operator[](std::size_t i) const { return data[i]; }
  constexpr auto& operator[](std::size_t i) { return data[i]; }

  constexpr void push(root_with_multiplicity<K> r)
  {
    if (count < max_roots)
      data[count++] = r;
  }
  constexpr void push(K value, std::size_t mult = 1) { push({value, mult}); }

  template<std::size_t M>
  constexpr void append(const roots_result<K, M>& other)
  {
    for (std::size_t i = 0; i < other.count && count < max_roots; ++i)
      data[count++] = other.data[i];
  }

  // Conversion to vector (for runtime convenience)
  operator std::vector<root_with_multiplicity<K>>() const
  { return std::vector<root_with_multiplicity<K>>{begin(), end()}; }

  constexpr std::vector<root_with_multiplicity<K>> to_vector() const
  { return std::vector<root_with_multiplicity<K>>{begin(), end()}; }

  // Extract just the root values
  constexpr std::array<K, max_roots> values() const
  {
    std::array<K, max_roots> vals{};
    for (std::size_t i = 0; i < count; ++i)
      vals[i] = data[i].value;
    return vals;
  }
};

} // end namespace lam::polynomial::univariate

// Export to lam::polynomial for convenient access
namespace lam::polynomial
{
export using univariate::polynomial_nttp;
export using univariate::is_approx_equal;
export using univariate::is_negligible;
export using univariate::get_epsilon;
} // end namespace lam::polynomial

// Export to lam for the simplest access
namespace lam
{
export using polynomial::univariate::polynomial_nttp;
export using polynomial::univariate::is_approx_equal;
export using polynomial::univariate::is_negligible;
export using polynomial::univariate::get_epsilon;
} // end namespace lam
