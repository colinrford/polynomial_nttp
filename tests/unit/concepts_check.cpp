/*
 *  concepts_check.cpp - written by Colin Ford
 *    see github.com/colinrford/polynomial_nttp for AGPL-3.0 License, and
 *                                              for more info
 *
 *  Compile-time verification of concepts.
 */
import std;
import lam.concepts;

using namespace lam::concepts::experimental;

// 1. Additive Group Checks
static_assert(additive_group_element_c_weak<int>, "int should be additive group");
static_assert(additive_group_element_c_weak<double>, "double should be additive group");
// 2. Multiplicative Group Checks
static_assert(multiplicative_group_element_c_weak<double>, "double should be multiplicative group");
// 3. Ring Checks
static_assert(ring_element_c_weak<int>, "int should be a ring");
static_assert(ring_element_c_weak<double>, "double should be a ring");
// 4. Field Checks
static_assert(field_element_c_weak<double>, "double should be a field");
static_assert(field_element_c_weak<int>, "int is weakly a field (due to integer division)");

int main() 
{
  return 0;
}
