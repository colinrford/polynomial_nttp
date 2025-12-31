# constructor
add_executable(constructor constructor.cpp)

target_include_directories(constructor PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(constructor polynomial_nttp lam::experimental_concepts)
target_compile_features(constructor PRIVATE cxx_std_23)

add_test(NAME Polynomial_Nttp_Constructor
         COMMAND constructor)

# initializer_list_constructor
add_executable(initializer_list_constructor init_list_constructor.cpp)

target_include_directories(initializer_list_constructor
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(initializer_list_constructor polynomial_nttp)
target_compile_features(initializer_list_constructor PRIVATE cxx_std_23)

add_test(NAME Polynomial_Nttp_Initializer_List_Constructor
         COMMAND initializer_list_constructor)

# constant_cast_constructor
add_executable(constant_cast_constructor constant_cast_constructor.cpp)

target_include_directories(constant_cast_constructor
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(constant_cast_constructor polynomial_nttp)
target_compile_features(constant_cast_constructor PRIVATE cxx_std_23)

add_test(NAME Polynomial_Nttp_Constant_Cast_Constructor
         COMMAND constant_cast_constructor)

# array_cast_constructor
add_executable(array_cast_constructor array_cast_constructor.cpp)

target_include_directories(array_cast_constructor
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(array_cast_constructor polynomial_nttp)
target_compile_features(array_cast_constructor PRIVATE cxx_std_23)

add_test(NAME Polynomial_Nttp_Array_Cast_Constructor
         COMMAND array_cast_constructor)

# iterators
add_executable(iterators iterators.cpp)
target_include_directories(iterators
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(iterators polynomial_nttp)
target_compile_features(iterators PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Iterators
         COMMAND iterators)

# subscript_operator
add_executable(subscript_operator subscript_operator.cpp)
target_include_directories(subscript_operator
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(subscript_operator polynomial_nttp)
target_compile_features(subscript_operator PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Subscript_operator[]"
         COMMAND subscript_operator)

# parentheses_operator
add_executable(parentheses_operator parentheses_operator.cpp)
target_include_directories(parentheses_operator
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(parentheses_operator polynomial_nttp)
target_compile_features(parentheses_operator PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Parentheses_operator()"
         COMMAND parentheses_operator)

# unary_negate_operator
add_executable(unary_negate_operator unary_negate_operator.cpp)
target_include_directories(unary_negate_operator
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(unary_negate_operator polynomial_nttp)
target_compile_features(unary_negate_operator PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Unary_Negate_operator-"
         COMMAND unary_negate_operator
)
# binary_add_polynomials
add_executable(binary_add_polynomials binary_add_polynomials.cpp)
target_include_directories(binary_add_polynomials
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_add_polynomials
  polynomial_nttp
  lam::experimental_concepts
)
target_compile_features(binary_add_polynomials PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Binary_Add_Polynomials_operator+"
         COMMAND binary_add_polynomials
)
# binary_add_polynomial_with_constant
add_executable(binary_add_polynomial_with_constant
        binary_add_polynomial_with_constant.cpp
)
target_include_directories(binary_add_polynomial_with_constant
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_add_polynomial_with_constant polynomial_nttp)
target_compile_features(binary_add_polynomial_with_constant PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Binary_Add_Polynomial_With_Constant_operator+"
         COMMAND binary_add_polynomial_with_constant)

# binary_subtract_polynomials
add_executable(binary_subtract_polynomials binary_subtract_polynomials.cpp)
target_include_directories(binary_subtract_polynomials
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_subtract_polynomials
  polynomial_nttp
  lam::experimental_concepts
)
target_compile_features(binary_subtract_polynomials PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Binary_Subtract_Polynomials_operator-"
         COMMAND binary_subtract_polynomials
)
# binary_subtract_polynomial_with_constant
add_executable(binary_subtract_polynomial_with_constant
        binary_subtract_polynomial_with_constant.cpp
)
target_include_directories(binary_subtract_polynomial_with_constant
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_subtract_polynomial_with_constant
  polynomial_nttp
  lam::experimental_concepts
)
target_compile_features(binary_subtract_polynomial_with_constant
  PRIVATE cxx_std_23
)
add_test(NAME
        "Polynomial_Nttp_Binary_Subtract_Polynomial_With_Constant_operator-"
         COMMAND binary_subtract_polynomial_with_constant)

# binary_multiply_polynomials
add_executable(binary_multiply_polynomials binary_multiply_polynomials.cpp)
target_include_directories(binary_multiply_polynomials
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_multiply_polynomials
  polynomial_nttp
  lam::experimental_concepts
)
target_compile_features(binary_multiply_polynomials PRIVATE cxx_std_23)
add_test(NAME "Polynomial_Nttp_Binary_Multiply_Polynomials_operator*"
         COMMAND binary_multiply_polynomials
)
# binary_multiply_polynomial_with_constant
add_executable(binary_multiply_polynomial_with_constant
        binary_multiply_polynomial_with_constant.cpp
)
target_include_directories(binary_multiply_polynomial_with_constant
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(binary_multiply_polynomial_with_constant polynomial_nttp)
target_compile_features(binary_multiply_polynomial_with_constant
  PRIVATE cxx_std_23
)
add_test(NAME
        "Polynomial_Nttp_Binary_Multiply_Polynomial_With_Constant_operator*"
         COMMAND binary_multiply_polynomial_with_constant)

# norm
add_executable(norm norm.cpp)
target_include_directories(norm PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(norm polynomial_nttp)
target_compile_features(norm PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Norm COMMAND norm)

# leading_coefficient
add_executable(leading_coefficient leading_coefficient.cpp)
target_include_directories(leading_coefficient PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(leading_coefficient polynomial_nttp)
target_compile_features(leading_coefficient PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Leading_Coefficient
         COMMAND leading_coefficient)

# make_monomial
add_executable(make_monomial make_monomial.cpp)
target_include_directories(make_monomial
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(make_monomial polynomial_nttp)
target_compile_features(make_monomial PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Make_Monomial
         COMMAND make_monomial)
# division
add_executable(division division.cpp)
target_include_directories(division
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(division polynomial_nttp)
target_compile_features(division PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Division
         COMMAND division)

add_executable(division_complex division_complex.cpp)
target_include_directories(division_complex
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(division_complex polynomial_nttp)
target_compile_features(division_complex PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Division_Complex
         COMMAND division_complex)
# derivative
add_executable(derivative derivative.cpp)
target_include_directories(derivative
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(derivative polynomial_nttp)
target_compile_features(derivative PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Derivative
         COMMAND derivative)

# antiderivative
add_executable(antiderivative antiderivative.cpp)
target_include_directories(antiderivative
  PRIVATE ${CMAKE_BINARY_DIR}
)
target_link_libraries(antiderivative polynomial_nttp)
target_compile_features(antiderivative PRIVATE cxx_std_23)
add_test(NAME Polynomial_Nttp_Anti_Derivative
         COMMAND antiderivative)

# concepts_check
add_executable(concepts_check concepts_check.cpp)
target_include_directories(concepts_check PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(concepts_check lam::experimental_concepts)
target_compile_features(concepts_check PRIVATE cxx_std_23)
add_test(NAME Concepts_Check COMMAND concepts_check)
