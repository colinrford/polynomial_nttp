
add_executable(wilkinson_stress wilkinson.cpp)
target_link_libraries(wilkinson_stress polynomial_nttp lam::concepts)
target_compile_features(wilkinson_stress PRIVATE cxx_std_23)
set_target_properties(wilkinson_stress PROPERTIES CXX_SCAN_FOR_MODULES ON)


add_executable(mignotte_stress mignotte.cpp)
target_link_libraries(mignotte_stress polynomial_nttp lam::concepts)
target_compile_features(mignotte_stress PRIVATE cxx_std_23)
set_target_properties(mignotte_stress PROPERTIES CXX_SCAN_FOR_MODULES ON)


add_executable(cyclotomic_stress cyclotomic.cpp)
target_link_libraries(cyclotomic_stress polynomial_nttp lam::concepts)
target_compile_features(cyclotomic_stress PRIVATE cxx_std_23)
set_target_properties(cyclotomic_stress PROPERTIES CXX_SCAN_FOR_MODULES ON)

add_test(NAME Cyclotomic_Stress_Test COMMAND cyclotomic_stress)

if(LAM_BUILD_THOUSAND_DIVISIONS)
  add_executable(thousand_divisions thousand_divisions.cpp)
  target_link_libraries(thousand_divisions polynomial_nttp lam::concepts)
  target_compile_features(thousand_divisions PRIVATE cxx_std_23)
  target_compile_options(thousand_divisions PRIVATE -fconstexpr-steps=200000000)
  set_target_properties(thousand_divisions PROPERTIES CXX_SCAN_FOR_MODULES ON)
endif()
