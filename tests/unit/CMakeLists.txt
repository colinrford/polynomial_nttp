# Helper function is now in cmake/PolynomialHelpers.cmake and included by root CMakeLists.txt

# ============================================================
# Structure Tests
# ============================================================

# ============================================================
# Structure Tests
# ============================================================

add_polynomial_test(constructor SOURCE structure/constructor.cpp LIBS lam::concepts)
add_polynomial_test(initializer_list_constructor SOURCE structure/init_list_constructor.cpp)
add_polynomial_test(constant_cast_constructor SOURCE structure/constant_cast_constructor.cpp)
add_polynomial_test(array_cast_constructor SOURCE structure/array_cast_constructor.cpp)

# ============================================================
# Accessor Tests
# ============================================================

add_polynomial_test(iterators SOURCE structure/iterators.cpp)
add_polynomial_test(subscript_operator SOURCE structure/subscript_operator.cpp TEST_NAME "Polynomial_Nttp_Subscript_operator[]")
add_polynomial_test(parentheses_operator SOURCE structure/parentheses_operator.cpp TEST_NAME "Polynomial_Nttp_Parentheses_operator()")

# ============================================================
# Unary Operator Tests
# ============================================================

add_polynomial_test(unary_negate_operator SOURCE structure/unary_negate_operator.cpp TEST_NAME "Polynomial_Nttp_Unary_Negate_operator-")

# ============================================================
# Binary Operator Tests: Addition
# ============================================================

add_polynomial_test(binary_add_polynomials
  SOURCE algebra/binary_add_polynomials.cpp
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Add_Polynomials_operator+")

add_polynomial_test(binary_add_polynomial_with_constant
  SOURCE algebra/binary_add_polynomial_with_constant.cpp
  TEST_NAME "Polynomial_Nttp_Binary_Add_Polynomial_With_Constant_operator+")

# ============================================================
# Binary Operator Tests: Subtraction
# ============================================================

add_polynomial_test(binary_subtract_polynomials
  SOURCE algebra/binary_subtract_polynomials.cpp
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Subtract_Polynomials_operator-")

add_polynomial_test(binary_subtract_polynomial_with_constant
  SOURCE algebra/binary_subtract_polynomial_with_constant.cpp
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Subtract_Polynomial_With_Constant_operator-")

# ============================================================
# Binary Operator Tests: Multiplication
# ============================================================

add_polynomial_test(binary_multiply_polynomials
  SOURCE algebra/binary_multiply_polynomials.cpp
  LIBS lam::concepts
  TEST_NAME "Polynomial_Nttp_Binary_Multiply_Polynomials_operator*")

add_polynomial_test(binary_multiply_polynomial_with_constant
  SOURCE algebra/binary_multiply_polynomial_with_constant.cpp
  TEST_NAME "Polynomial_Nttp_Binary_Multiply_Polynomial_With_Constant_operator*")

add_polynomial_test(binary_multiply_fft SOURCE fft/binary_multiply_fft.cpp)
add_polynomial_test(ct_fft_benchmark SOURCE fft/ct_fft_benchmark.cpp)
add_polynomial_test(precision_benchmark SOURCE fft/precision_benchmark.cpp)
add_polynomial_test(real_world_signal_fft SOURCE fft/real_world_signal_fft.cpp)

# FFT tests require increased constexpr step limits for compile-time evaluation
# 16384-point FFT requires ~1 billion steps
target_compile_options(real_world_signal_fft PRIVATE -fconstexpr-steps=1000000000)
# 1024-point FFT round-trip requires ~100 million steps
target_compile_options(precision_benchmark PRIVATE -fconstexpr-steps=100000000)
# Binary multiply with FFT uses fixed-size polynomials
target_compile_options(binary_multiply_fft PRIVATE -fconstexpr-steps=100000000)

# ============================================================
# Algebra Tests
# ============================================================

add_polynomial_test(norm SOURCE algebra/norm.cpp)
add_polynomial_test(leading_coefficient SOURCE algebra/leading_coefficient.cpp)
add_polynomial_test(make_monomial SOURCE algebra/make_monomial.cpp)

# ============================================================
# Division Tests
# ============================================================

add_polynomial_test(division SOURCE algebra/division.cpp)
add_polynomial_test(division_complex SOURCE algebra/division_complex.cpp)

# ============================================================
# Calculus Tests
# ============================================================

add_polynomial_test(derivative SOURCE algebra/derivative.cpp)
add_polynomial_test(antiderivative SOURCE algebra/antiderivative.cpp TEST_NAME "Polynomial_Nttp_Anti_Derivative")

# ============================================================
# Advanced Tests
# ============================================================

add_polynomial_test(concepts_check
  SOURCE concepts/concepts_check.cpp
  LIBS lam::concepts
  TEST_NAME "Concepts_Check")

add_polynomial_test(roots SOURCE roots/roots.cpp)

# ============================================================
# Coverage Gap Tests
# ============================================================

add_polynomial_test(coverage_gap_fillers SOURCE coverage_gap_fillers.cpp)

# ctbignum-based tests
find_package(lam_ctbignum QUIET)
if(TARGET lam::ctbignum)
  get_target_property(_aliased lam::ctbignum ALIASED_TARGET)
  if(_aliased)
    set(_ctbignum_target ${_aliased})
  else()
    set(_ctbignum_target lam::ctbignum)
  endif()
  set_property(TARGET ${_ctbignum_target} APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_23)
  add_polynomial_test(berlekamp SOURCE berlekamp/berlekamp.cpp LIBS lam::ctbignum)
  
endif()
