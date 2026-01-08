/*
 *  polynomial_nttp-univariate-roots.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  polynomial_nttp is a c++ module
 */

module;

export module lam.polynomial_nttp:univariate.roots;

import std;
import lam.concepts;
import :univariate.structure;
import :univariate.algebra;
import :univariate.berlekamp;

namespace stdr = std::ranges;
namespace stdv = std::views;

namespace lam::polynomial::univariate::roots
{
template<typename K>
concept field_element_c_weak = lam::concepts::experimental::field_element_c_weak<K>;

// Imported from structure module
using univariate::root_with_multiplicity;
using univariate::roots_result;

// Generic sqrt iterator step: x_{k+1} = 0.5 * (x_k + S / x_k)
template<field_element_c_weak K>
constexpr K sqrt_iterator(K of, K x_k)
{
  constexpr K one_half = K(1) / K(2);
  return one_half * (x_k + of / x_k);
}

// Generic cbrt iterator step: x_{k+1} = (2*x_k + S / x_k²) / 3
template<field_element_c_weak K>
constexpr K cbrt_iterator(K of, K x_k)
{
  return (K(2) * x_k + of / (x_k * x_k)) / K(3);
}

// Fixed point solver for sqrt using Newton-Raphson (Compile-time friendly)
template<field_element_c_weak K>
constexpr K sqrt_fixed_point_solver(K of)
{
  if (is_negligible(of))
    return K(0);

  // Handle Branch Cut for Complex Types (Negative Reals)
  if constexpr (requires(K k) {
                  k.real();
                  k.imag();
                  typename K::value_type;
                })
  {
    using value_t = typename K::value_type;
    auto r = of.real();
    auto i_mag = of.imag();
    // If imaginary part is negligible and real part is negative
    if (is_negligible(i_mag) && r < value_t(0))
    { // sqrt(-r) * i
      K pos_input = K(-r);
      K pos_root = sqrt_fixed_point_solver(pos_input);
      return K(0, pos_root.real());
    }
  }

  K running_root = K(1);
  K prev = K(0); // Different from running_root to ensure at least one iteration

  // Limit iterations for safety in constexpr context
  constexpr int max_iter = 100;
  for (int i = 0; i < max_iter; ++i)
  {
    if (is_approx_equal(running_root, prev))
      break;
    prev = running_root;
    running_root = sqrt_iterator(of, running_root);
  }
  return running_root;
}

// Fixed point solver for cube root using Newton-Raphson
template<field_element_c_weak K>
constexpr K cbrt_fixed_point_solver(K of)
{
  if (is_negligible(of))
    return K(0);

  // Handle negative reals for real types
  if constexpr (std::is_floating_point_v<K>)
  {
    if (of < K(0))
      return -cbrt_fixed_point_solver(-of);
  }

  K running_root = K(1);
  K prev = K(0);

  constexpr int max_iter = 100;
  for (int i = 0; i < max_iter; ++i)
  {
    if (is_approx_equal(running_root, prev))
      break;
    prev = running_root;
    running_root = cbrt_iterator(of, running_root);
  }
  return running_root;
}

// Helper concepts for feature detection
template<typename T>
concept has_optional_cbrt = requires(T x) {
  { cbrt(x) } -> std::same_as<std::optional<T>>;
};

// Generic cbrt dispatch (compile-time vs runtime)
template<field_element_c_weak K>
constexpr std::optional<K> generic_cbrt(K val)
{
  if constexpr (has_optional_cbrt<K>)
  {
    return cbrt(val);
  }

  if (std::is_constant_evaluated())
    return cbrt_fixed_point_solver(val);
  else
  {
    if constexpr (std::is_floating_point_v<K>)
      return std::cbrt(val);
    else
      return cbrt_fixed_point_solver(val);
  }
}

// Brute force root solver for small fields or specific ranges
// Useful for Characteristic 2 and 3 where analytic formulas fail or small fields
export template<typename K, std::size_t N, typename Range>
constexpr auto roots_brute_force(const polynomial_nttp<K, N>& p, Range&& range) -> roots_result<K, N>
{
  roots_result<K, N> result;
  for (const auto& val : range)
  {
    if (is_negligible(p(val)))
    {
      result.push(val, 1);
    }
  }
  return result;
}

export template<field_element_c_weak K>
constexpr auto roots_degree_1(const polynomial_nttp<K, 1>& p) -> roots_result<K, 1>
{
  roots_result<K, 1> result;
  if (is_negligible(p[1]))
    return result;
  result.push(-p[0] / p[1], 1);
  return result;
}

// Helper concepts for feature detection
template<typename T>
concept has_optional_sqrt = requires(T x) {
  { sqrt(x) } -> std::same_as<std::optional<T>>;
};

// Solver for ax^2 + bx + c = 0 with flexible sqrt dispatch
export template<field_element_c_weak K>
constexpr auto roots_degree_2(const polynomial_nttp<K, 2>& p) -> roots_result<K, 2>
{
  roots_result<K, 2> result;
  auto a = p[2];
  auto b = p[1];
  auto c = p[0];

  if (is_negligible(a))
  {
    if (!is_negligible(b))
      result.push(-c / b, 1);
    return result;
  }

  auto discriminant = b * b - K(4) * a * c;

  auto generic_sqrt = [](const auto& val) -> std::optional<K> {
    // Use native optional sqrt if available (e.g. finite fields)
    if constexpr (has_optional_sqrt<K>)
    {
      return sqrt(val);
    }

    if (std::is_constant_evaluated())
    {
      if constexpr (std::is_floating_point_v<K>)
        if (val < K(0))
          return std::nullopt;
      return sqrt_fixed_point_solver(val);
    }

    if constexpr (std::is_floating_point_v<K>)
    {
      if (val < K(0))
        return std::nullopt;
      return std::sqrt(val);
    }
    else
    {
      using std::sqrt;
      return sqrt(val);
    }
  };

  auto root_delta_opt = generic_sqrt(discriminant);

  if (!root_delta_opt)
    return result;

  K root_delta = *root_delta_opt;
  auto two_a = K(2) * a;

  if (is_negligible(discriminant))
  {
    // Double root
    result.push(-b / two_a, 2);
  }
  else
  {
    // Two distinct roots
    result.push((-b + root_delta) / two_a, 1);
    result.push((-b - root_delta) / two_a, 1);
  }

  return result;
}

// -----------------------------------------------------------------------------
// Constexpr Trigonometric Helpers (Minimax Approximations)
// -----------------------------------------------------------------------------
namespace math
{
template<typename Real>
constexpr Real PI = Real(3.14159265358979323846);
template<typename Real>
constexpr Real PI_2 = Real(1.57079632679489661923);

// Range reduction to [-PI, PI]
template<typename Real>
constexpr Real reduce_to_pi(Real x)
{
  constexpr Real pi = PI<Real>;
  constexpr Real two_pi = Real(2) * pi;
  while (x > pi)
    x -= two_pi;
  while (x < -pi)
    x += two_pi;
  return x;
}

template<typename Real>
constexpr Real cos_core(Real x)
{
  Real x2 = x * x;
  return Real(1.0) + x2 * (Real(-0.49999999999999994) +
                           x2 * (Real(0.04166666666666650) +
                                 x2 * (Real(-0.00138888888888667) +
                                       x2 * (Real(0.00002480158730105) +
                                             x2 * (Real(-0.00000027557319119) + x2 * Real(0.00000000208757209))))));
}

template<typename Real>
constexpr Real sin_core(Real x)
{
  Real x2 = x * x;
  return x *
         (Real(1.0) + x2 * (Real(-0.16666666666666666) +
                            x2 * (Real(0.00833333333333333) +
                                  x2 * (Real(-0.00019841269841270) +
                                        x2 * (Real(0.00000275573192240) +
                                              x2 * (Real(-0.00000002505210840) + x2 * Real(0.00000000016059044)))))));
}

template<typename Real>
constexpr Real cos(Real x)
{
  x = reduce_to_pi(x);
  if (x < Real(0))
    x = -x;
  constexpr Real pi = PI<Real>;

  if (x <= pi / Real(4))
    return cos_core(x);
  else if (x <= Real(3) * pi / Real(4))
    return -sin_core(x - PI_2<Real>);
  else
    return -cos_core(pi - x);
}

template<typename Real>
constexpr Real sin(Real x)
{
  return cos(x - PI_2<Real>);
}

template<typename Real>
constexpr Real acos(Real x)
{
  // Newton-Raphson: f(y) = cos(y) - x = 0
  if (x > Real(1))
    x = Real(1);
  if (x < Real(-1))
    x = Real(-1);

  Real y = PI_2<Real> - x; // Initial guess
  constexpr int max_iter = 25;

  for (int i = 0; i < max_iter; ++i)
  {
    Real cy = cos(y);
    Real sy = sin(y);

    // Use standard epsilon-based comparison
    if (is_approx_equal(cy, x))
      return y;

    // Avoid division by zero
    if (is_negligible(sy))
      break;

    y = y + (cy - x) / sy;
  }
  return y;
}
} // namespace math

// Solver for ax³ + bx² + cx + d = 0 using Cardano's formula
export template<field_element_c_weak K>
constexpr auto roots_degree_3(const polynomial_nttp<K, 3>& p) -> roots_result<K, 3>
{
  roots_result<K, 3> result;

  auto a = p[3];
  auto b = p[2];
  auto c = p[1];
  auto d = p[0];

  // Degenerate case: not actually cubic
  if (is_negligible(a))
  {
    polynomial_nttp<K, 2> reduced{{d, c, b}};
    auto deg2_result = roots_degree_2(reduced);
    result.append(deg2_result);
    return result;
  }

  // Convert to depressed cubic t³ + pt + q = 0 via x = t - b/(3a)
  auto shift = b / (K(3) * a);
  auto p_coef = (K(3) * a * c - b * b) / (K(3) * a * a);
  auto q_coef = (K(2) * b * b * b - K(9) * a * b * c + K(27) * a * a * d) / (K(27) * a * a * a);

  // Discriminant: Δ = -4p³ - 27q²
  auto discriminant = -K(4) * p_coef * p_coef * p_coef - K(27) * q_coef * q_coef;

  if (is_negligible(p_coef) && is_negligible(q_coef))
  { // Triple root at 0 (before shift)
    result.push(-shift, 3);
    return result;
  }

  if (is_negligible(discriminant))
  { // Repeated roots
    if (is_negligible(p_coef))
    {
      result.push(-shift, 3);
    }
    else
    { // One single root, one double root
      auto t1 = K(3) * q_coef / p_coef;
      auto t2 = -K(3) * q_coef / (K(2) * p_coef);
      result.push(t1 - shift, 1);
      result.push(t2 - shift, 2);
    }
    return result;
  }

  // For real coefficients with Δ > 0: three distinct real roots
  // Use trigonometric method
  bool use_cardano = true;

  if constexpr (std::is_floating_point_v<K>)
  {
    if (discriminant > K(0) && p_coef < K(0))
    {
      // Check if we can use Trig method
      // We can ALWAYS use it now thanks to math:: helpers for compile-time!

      K m;
      if (std::is_constant_evaluated())
        m = K(2) * sqrt_fixed_point_solver(-p_coef / K(3));
      else
        m = K(2) * std::sqrt(-p_coef / K(3));

      auto cos_arg = (K(3) * q_coef) / (p_coef * m);
      if (cos_arg > K(1))
        cos_arg = K(1);
      if (cos_arg < K(-1))
        cos_arg = K(-1);

      K theta;
      constexpr K pi = K(3.14159265358979323846);

      if (std::is_constant_evaluated())
      {
        theta = math::acos(cos_arg);
        result.push(m * math::cos(theta / K(3)) - shift, 1);
        result.push(m * math::cos((theta - K(2) * pi) / K(3)) - shift, 1);
        result.push(m * math::cos((theta + K(2) * pi) / K(3)) - shift, 1);
      }
      else
      {
        theta = std::acos(cos_arg);
        result.push(m * std::cos(theta / K(3)) - shift, 1);
        result.push(m * std::cos((theta - K(2) * pi) / K(3)) - shift, 1);
        result.push(m * std::cos((theta + K(2) * pi) / K(3)) - shift, 1);
      }
      use_cardano = false;
    }
  }

  if (use_cardano)
  {
    auto inner = q_coef * q_coef / K(4) + p_coef * p_coef * p_coef / K(27);
    K sqrt_inner;
    if constexpr (std::is_floating_point_v<K>)
    {
      if (inner < K(0))
        return result; // Should be handled by trig block
      sqrt_inner = std::is_constant_evaluated() ? sqrt_fixed_point_solver(inner) : std::sqrt(inner);
    }
    else
    {
      sqrt_inner = sqrt_fixed_point_solver(inner);
    }

    auto u_opt = generic_cbrt(-q_coef / K(2) + sqrt_inner);
    auto v_opt = generic_cbrt(-q_coef / K(2) - sqrt_inner);

    if (u_opt && v_opt)
    {
      auto t = *u_opt + *v_opt;
      result.push(t - shift, 1);
    }
  }

  return result;
}

// Solver for quartic using Ferrari's Method
export template<field_element_c_weak K>
constexpr auto roots_degree_4(const polynomial_nttp<K, 4>& p) -> roots_result<K, 4>
{
  roots_result<K, 4> result;

  // Handle a=0 case (cubic fallback)
  if (is_negligible(p[4]))
  {
    polynomial_nttp<K, 3> reduced{{p[0], p[1], p[2], p[3]}};
    auto deg3 = roots_degree_3(reduced);
    result.append(deg3);
    return result;
  }

  // Normalize: x^4 + ax^3 + bx^2 + cx + d = 0
  K lc = p[4];
  K a = p[3] / lc;
  K b = p[2] / lc;
  K c = p[1] / lc;
  K d = p[0] / lc;

  // Depress: x = y - a/4
  // y^4 + Py^2 + Qy + R = 0
  K shift = a / K(4);
  K P = b - K(3) * a * a / K(8);
  K Q = c + a * a * a / K(8) - a * b / K(2);
  K R = d - K(3) * a * a * a * a / K(256) + a * a * b / K(16) - a * c / K(4);

  // Helper for sqrt dispatch
  auto generic_sqrt = [](const auto& val) -> std::optional<K> {
    if constexpr (has_optional_sqrt<K>)
    {
      return sqrt(val);
    }
    if (std::is_constant_evaluated())
    {
      if constexpr (std::is_floating_point_v<K>)
        if (val < K(0))
          return std::nullopt;
      return sqrt_fixed_point_solver(val);
    }
    if constexpr (std::is_floating_point_v<K>)
    {
      if (val < K(0))
        return std::nullopt;
      return std::sqrt(val);
    }
    else
    {
      using std::sqrt;
      return sqrt(val);
    }
  };

  // Case 1: Biquadratic (Q = 0)
  if (is_negligible(Q))
  {
    // u^2 + Pu + R = 0 where u = y^2
    polynomial_nttp<K, 2> quad_u{{R, P, K(1)}};
    auto u_roots = roots_degree_2(quad_u);

    for (const auto& u_root : u_roots)
    {
      // y^2 = u (multiplicity logic handled by returning separate entries)
      if (is_negligible(u_root.value))
      {
        result.push(-shift, 2 * u_root.multiplicity);
      }
      else
      {
        auto sqrt_u_opt = generic_sqrt(u_root.value);
        if (sqrt_u_opt)
        {
          result.push(*sqrt_u_opt - shift, u_root.multiplicity);
          result.push(-(*sqrt_u_opt) - shift, u_root.multiplicity);
        }
      }
    }
    return result;
  }

  // Case 2: Ferrari's Resolvent Cubic
  // z^3 + 2Pz^2 + (P^2 - 4R)z - Q^2 = 0
  polynomial_nttp<K, 3> resolvent;
  resolvent.coefficients[3] = K(1);
  resolvent.coefficients[2] = K(2) * P;
  resolvent.coefficients[1] = P * P - K(4) * R;
  resolvent.coefficients[0] = -Q * Q;

  auto z_roots = roots_degree_3(resolvent);

  // Find a suitable non-zero root z0
  K z0 = K(0);
  bool found_z0 = false;

  for (const auto& r : z_roots)
  {
    if (!is_negligible(r.value))
    {
      // Prefer positive real for safety
      if constexpr (std::is_floating_point_v<K>)
      {
        if (r.value > K(0))
        {
          z0 = r.value;
          found_z0 = true;
          break;
        }
        if (!found_z0)
        {
          z0 = r.value;
          found_z0 = true;
        }
      }
      else
      {
        z0 = r.value;
        found_z0 = true;
        break;
      }
    }
  }

  if (!found_z0)
    return result;

  auto sqrt_z0_opt = generic_sqrt(z0);
  if (!sqrt_z0_opt)
    return result;
  K sqrt_z0 = *sqrt_z0_opt;

  // Two quadratics:
  // y^2 + sqrt(z0)y + (P + z0 - Q/sqrt(z0))/2 = 0
  // y^2 - sqrt(z0)y + (P + z0 + Q/sqrt(z0))/2 = 0

  K term1 = (P + z0 - Q / sqrt_z0) / K(2);
  K term2 = (P + z0 + Q / sqrt_z0) / K(2);

  polynomial_nttp<K, 2> q1{{term1, sqrt_z0, K(1)}};
  polynomial_nttp<K, 2> q2{{term2, -sqrt_z0, K(1)}};

  auto roots1 = roots_degree_2(q1);
  auto roots2 = roots_degree_2(q2);

  for (auto& r : roots1)
    result.push(r.value - shift, r.multiplicity);
  for (auto& r : roots2)
    result.push(r.value - shift, r.multiplicity);

  return result;
}

template<field_element_c_weak K, std::size_t N>
constexpr auto root_newton_raphson(const polynomial_nttp<K, N>& p, K initial_guess = K(0)) -> std::optional<K>
{
  if constexpr (N == 0)
    return std::nullopt;

  auto deriv = derivative(p);
  K x = initial_guess;
  constexpr int max_iters = 100;

  for (int i = 0; i < max_iters; ++i)
  {
    K y = p(x);
    if (is_negligible(y, 1e-6))
      return x;

    K dy = deriv(x);
    if (is_negligible(dy, 1e-6))
      return std::nullopt;
    x = x - y / dy;
  }
  if (is_negligible(p(x), 1e-6))
    return x;
  return std::nullopt;
}

// Try Newton-Raphson with multiple starting points
template<field_element_c_weak K, std::size_t N>
constexpr auto root_newton_raphson_multi(const polynomial_nttp<K, N>& p) -> std::optional<K>
{
  if constexpr (N == 0)
    return std::nullopt;

  // Try multiple starting points
  std::array<K, 5> starts;
  if constexpr (std::is_floating_point_v<K>)
    starts = {K(0), K(1), K(-1), K(0.5), K(-0.5)};
  else
    starts = {K(0), K(1), K(-1), K(2), K(-2)};

  for (const auto& start : starts)
  {
    auto result = root_newton_raphson(p, start);
    if (result)
      return result;
  }
  return std::nullopt;
}

// Check if a root has higher multiplicity by testing derivative
template<field_element_c_weak K, std::size_t N>
constexpr std::size_t compute_multiplicity(const polynomial_nttp<K, N>& p, K root)
{
  std::size_t mult = 1;

  if constexpr (N >= 1)
  {
    auto d1 = derivative(p);
    if (is_negligible(d1(root), 1e-6))
    {
      mult = 2;
      if constexpr (N >= 2)
      {
        auto d2 = derivative(d1);
        if (is_negligible(d2(root), 1e-6))
          mult = 3;
      }
    }
  }
  return mult;
}

// Synthetic division: divide p by (x - root)
template<field_element_c_weak K, std::size_t N>
constexpr polynomial_nttp<K, N - 1> deflate(const polynomial_nttp<K, N>& p, K root)
{
  polynomial_nttp<K, N - 1> quotient;
  quotient.coefficients[N - 1] = p[N];
  for (std::size_t i = N - 1; i > 0; --i)
    quotient.coefficients[i - 1] = p[i] + root * quotient.coefficients[i];
  return quotient;
}

// Forward declaration
export template<field_element_c_weak K, std::size_t N>
constexpr roots_result<K, N> roots(const polynomial_nttp<K, N>& p);

// Generic solver using Newton-Raphson and Deflation
template<field_element_c_weak K, std::size_t N>
constexpr roots_result<K, N> roots_via_newton(const polynomial_nttp<K, N>& p)
{
  roots_result<K, N> result;

  // Use Newton-Raphson to find one root
  auto r_opt = root_newton_raphson_multi(p);
  if (!r_opt)
    return result;

  K root = *r_opt;
  std::size_t mult = compute_multiplicity(p, root);

  result.push(root, mult);

  // Deflate the polynomial by (x - root)^mult
  if constexpr (N > 1)
  {
    auto quotient = deflate(p, root);

    // Deflate additional times for multiplicity > 1
    if constexpr (N > 2)
    {
      if (mult >= 2)
      {
        auto q2 = deflate(quotient, root);
        if constexpr (N > 3)
        {
          if (mult >= 3)
          {
            auto q3 = deflate(q2, root);
            auto rec_roots = roots(q3);
            result.append(rec_roots);
            return result;
          }
        }
        auto rec_roots = roots(q2);
        result.append(rec_roots);
        return result;
      }
    }

    auto rec_roots = roots(quotient);
    result.append(rec_roots);
  }

  return result;
}

// Generic Dispatcher with Deflation and Multiplicity Tracking
export template<field_element_c_weak K, std::size_t N>
constexpr roots_result<K, N> roots(const polynomial_nttp<K, N>& p)
{
  roots_result<K, N> result;

  // ============================================================
  // Finite Field Dispatch (Prime Characteristic)
  // ============================================================
  if constexpr (univariate::finite_field_traits<K>::is_finite_field)
  {
    constexpr auto P = univariate::finite_field_traits<K>::modulus;
    return lam::polynomial::univariate::berlekamp::roots_berlekamp<K, P, N>(
      p, K(0), K(1));
  }

  if constexpr (N == 0)
    return result;
  else if constexpr (N == 1)
  {
    auto deg1 = roots_degree_1(p);
    result.append(deg1);
    return result;
  }
  else if constexpr (N == 2)
  {
    auto deg2 = roots_degree_2(p);
    result.append(deg2);
    return result;
  }
  else if constexpr (N == 3)
  {
    // N=3: Analytic solution (Cardano/Trigonometric)
    // Works at runtime and compile-time (via constexpr math:: helpers)
    auto deg3 = roots_degree_3(p);
    result.append(deg3);
    return result;
  }
  else if constexpr (N == 4)
  {
    auto deg4 = roots_degree_4(p);
    result.append(deg4);
    return result;
  }
  else
  {
    // For degree >= 5, use Newton-Raphson with deflation
    return roots_via_newton(p);
  }
}

} // end namespace lam::polynomial::univariate::roots

namespace lam::polynomial
{
  export using univariate::roots::root_with_multiplicity;
  export using univariate::roots::roots_result;
  export using univariate::roots::roots_degree_1;
  export using univariate::roots::roots_degree_2;
  export using univariate::roots::roots_degree_3;
  export using univariate::roots::roots_degree_4;
  export using univariate::roots::roots_brute_force;
  export using univariate::roots::roots;
}

// Export to lam for the simplest access
namespace lam
{
  export using polynomial::univariate::roots::roots;
  export using polynomial::univariate::roots::roots_result;
  export using polynomial::univariate::roots::root_with_multiplicity;
}
