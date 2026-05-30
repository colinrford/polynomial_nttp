add_polynomial_executable(wilkinson_stress SOURCE wilkinson.cpp LIBS lam_concepts::concepts)
add_polynomial_executable(mignotte_stress SOURCE mignotte.cpp LIBS lam_concepts::concepts)

add_polynomial_test(cyclotomic_stress SOURCE cyclotomic.cpp TEST_NAME Cyclotomic_Stress_Test LIBS lam_concepts::concepts)
add_polynomial_test(lagrange_stress SOURCE lagrange.cpp TEST_NAME Lagrange_Stress_Test LIBS lam_concepts::concepts)

if(LAM_BUILD_THOUSAND_DIVISIONS)
  include(LamConstexprLimit)
  add_polynomial_executable(thousand_divisions SOURCE thousand_divisions.cpp LIBS lam_concepts::concepts)
  lam_target_constexpr_limit(thousand_divisions 200000000)
endif()
