# Helper function is now in cmake/PolynomialHelpers.cmake and included by root CMakeLists.txt

# ============================================================
# Structure Tests
# ============================================================

add_polynomial_test(constructor LIBS lam::concepts)
add_polynomial_test(initializer_list_constructor SOURCE init_list_constructor.cpp)
add_polynomial_test(constant_cast_constructor)
add_polynomial_test(array_cast_constructor)

# ============================================================
# Accessor Tests
# ============================================================

add_polynomial_test(iterators)
add_polynomial_test(subscript_operator TEST_NAME "Polynomial_Nttp_Subscript_operator[]")
add_polynomial_test(parentheses_operator TEST_NAME "Polynomial_Nttp_Parentheses_operator()")

# ============================================================
# Unary Operator Tests
# ============================================================

add_polynomial_test(unary_negate_operator TEST_NAME "Polynomial_Nttp_Unary_Negate_operator-")

# ============================================================
# Binary Operator Tests: Addition
# ============================================================

add_polynomial_test(binary_add_polynomials
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Add_Polynomials_operator+")

add_polynomial_test(binary_add_polynomial_with_constant
  TEST_NAME "Polynomial_Nttp_Binary_Add_Polynomial_With_Constant_operator+")

# ============================================================
# Binary Operator Tests: Subtraction
# ============================================================

add_polynomial_test(binary_subtract_polynomials
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Subtract_Polynomials_operator-")

add_polynomial_test(binary_subtract_polynomial_with_constant
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Subtract_Polynomial_With_Constant_operator-")

# ============================================================
# Binary Operator Tests: Multiplication
# ============================================================

add_polynomial_test(binary_multiply_polynomials
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Multiply_Polynomials_operator*")

add_polynomial_test(binary_multiply_polynomial_with_constant
  TEST_NAME "Polynomial_Nttp_Binary_Multiply_Polynomial_With_Constant_operator*")

# ============================================================
# Algebra Tests
# ============================================================

add_polynomial_test(norm)
add_polynomial_test(leading_coefficient)
add_polynomial_test(make_monomial)

# ============================================================
# Division Tests
# ============================================================

add_polynomial_test(division)
add_polynomial_test(division_complex)

# ============================================================
# Calculus Tests
# ============================================================

add_polynomial_test(derivative)
add_polynomial_test(antiderivative TEST_NAME "Polynomial_Nttp_Anti_Derivative")

# ============================================================
# Advanced Tests
# ============================================================

add_polynomial_test(concepts_check
  LIBS lam::concepts
  TEST_NAME "Concepts_Check")

add_polynomial_test(roots)

# ctbignum-based tests
find_package(lam_ctbignum QUIET)
if(TARGET lam::ctbignum)
  set_property(TARGET lam::ctbignum APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_23)
  add_polynomial_test(berlekamp LIBS lam::ctbignum)
  
  # roots_finite_field requires lam_interop which needs ctbignum AND bitvector
  if(LAM_INTEROP_AVAILABLE)
    add_polynomial_test(roots_finite_field LIBS lam::interop)
  else()
    message(STATUS "Skipping roots_finite_field test (requires lam_interop)")
  endif()
endif()
