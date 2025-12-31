/*
 *  concepts_check.cpp
 *  Compile-time verification of concepts.
 */
import std;
import lam.experimental.concepts;

using namespace lam::experimental::concepts;

// 1. Additive Group Checks
static_assert(additive_group_element_c_weak<int>, "int should be additive group");
static_assert(additive_group_element_c_weak<double>, "double should be additive group");

// 2. Multiplicative Group Checks
static_assert(multiplicative_group_element_c_weak<double>, "double should be multiplicative group");
// int has no inverse (1/2 is 0), so strictly it usually fails field concepts, 
// but does it fail this weak concept? 
// '1/g' is valid syntax for int, returning int. 
// So int MIGHT pass this "weak" concept even though it shouldn't mathematically.
// Let's verify what happens.
// If valid syntax is all that matters, int will pass.

// 3. Ring Checks
static_assert(ring_element_c_weak<int>, "int should be a ring");
static_assert(ring_element_c_weak<double>, "double should be a ring");

// 4. Field Checks
static_assert(field_element_c_weak<double>, "double should be a field");
// int fails strict field check? 
// It passes ring. It arguably passes multiplicative_group (1/a syntax).
// So 'int' might be a weak field.
static_assert(field_element_c_weak<int>, "int is weakly a field (due to integer division)");

int main() 
{
  return 0;
}
