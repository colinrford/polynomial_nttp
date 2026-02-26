add_polynomial_executable(wilkinson_stress SOURCE wilkinson.cpp LIBS lam::concepts)
add_polynomial_executable(mignotte_stress SOURCE mignotte.cpp LIBS lam::concepts)

add_polynomial_test(cyclotomic_stress SOURCE cyclotomic.cpp TEST_NAME Cyclotomic_Stress_Test LIBS lam::concepts)
add_polynomial_test(lagrange_stress SOURCE lagrange.cpp TEST_NAME Lagrange_Stress_Test LIBS lam::concepts)

if(LAM_BUILD_THOUSAND_DIVISIONS)
  add_polynomial_executable(thousand_divisions SOURCE thousand_divisions.cpp LIBS lam::concepts)
  target_compile_options(thousand_divisions PRIVATE -fconstexpr-steps=200000000)
endif()
