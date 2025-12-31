/*
 *  concepts.cppm – written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *  experimental.concepts is a c++ module
 *    as the name suggests, the contents found here are strictly experimental,
 *    and while the project will compile and run just fine, these concepts are
 *    still very much in development
 *
 *  these concepts are part of another library I am working on, and are
 *    subject to reappear within that library under a different license,
 *    whenever said library is released
 */

module;

import std;

export module lam.experimental.concepts;

/*
 *  all concepts will be declared with _c to distinguish them from types, etc
 *    at this time I am also appending _weak in many situations as concepts
 *    are limited in what they can check. With time it would be nice to
 *    work toward removing the _weak appendage.
 */
namespace lam::experimental::concepts::internals
{
  template<typename T>
  concept has_additive_identity_c_weak = 
  requires {
    { T(0) } -> std::same_as<T>;
  } || requires {
    { T::zero() } -> std::same_as<T>;
  };
  template<typename T>
  concept has_multiplicative_identity_c_weak = 
  requires {
    { T(1) } -> std::same_as<T>;
  } || requires {
    { T::one() } -> std::same_as<T>;
  };
  /*
   *  operator overloaded weak syntactical 'requirements' modeling group
   *  elements (does not check values)
   */
  template<typename G>
  concept additive_group_element_c_weak = has_additive_identity_c_weak<G> and 
  requires(G g, G h, G k)
  {
    { -g } -> std::same_as<decltype(g)>;
    { g + h } -> std::same_as<decltype(h + g)>;
    { g - h } -> std::same_as<decltype(h - g)>;
    { (g + h) + k } -> std::same_as<decltype(g + (h + k))>;
  };

  template<typename G>
  concept multiplicative_group_element_c_weak = has_multiplicative_identity_c_weak<G> and 
  requires(G g, G h, G k)
  {
    { g * h } -> std::same_as<decltype(h * g)>;
    { (g * h) * k } -> std::same_as<decltype(g * (h * k))>;
    { g / h } -> std::same_as<decltype(g)>;
    requires requires { { 1 / g } -> std::same_as<decltype(g)>; } 
          or requires { { G::one() / g } -> std::same_as<decltype(g)>; }
          or requires { { inv(g) } -> std::same_as<decltype(g)>; }; 
  };

  template<typename G>
  concept group_element_c_weak = additive_group_element_c_weak<G> or multiplicative_group_element_c_weak<G>;

  /*
   *  operator overloaded weak syntactical 'requirements' modeling group
   *  elements (does not check values)
   */
  template<typename R>
  concept ring_element_c_weak = additive_group_element_c_weak<R> 
                             and has_multiplicative_identity_c_weak<R>
                             and requires(R r, R s, R t)
  {
    { r * s } -> std::same_as<decltype(s * r)>;
    { (r * s) * t } -> std::same_as<decltype(r * (s * t))>;
    { r * (s + t) } -> std::same_as<decltype(r * s + r * t)>;
    { (r + s) * t } -> std::same_as<decltype(r * t + s * t)>;
  };

  /*
   *  operator overloaded weak syntactical 'requirements' modeling field
   *  elements (does not check values)
   */
  template<typename K>
  concept field_element_c_weak = ring_element_c_weak<K> 
                              and multiplicative_group_element_c_weak<K>;
} // end namespace lam::experimental::concepts::internals
namespace lam::experimental::concepts
{
  export template<typename G>
  concept additive_group_element_c_weak = internals::additive_group_element_c_weak<G>;
  export template<typename G>
  concept multiplicative_group_element_c_weak = internals::multiplicative_group_element_c_weak<G>;
  export template<typename G>
  concept group_element_c_weak = internals::group_element_c_weak<G>;
  export template<typename R>
  concept ring_element_c_weak = internals::ring_element_c_weak<R>;
  export template<typename K>
  concept field_element_c_weak = internals::field_element_c_weak<K>;
} // end namespace lam::experimental::concepts
