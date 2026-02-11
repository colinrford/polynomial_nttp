include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(test_polynomial_nttp test_polynomial_nttp.cpp)

target_link_libraries(test_polynomial_nttp polynomial_nttp lam::concepts)
target_compile_features(test_polynomial_nttp PRIVATE cxx_std_23)
set_target_properties(test_polynomial_nttp PROPERTIES CXX_SCAN_FOR_MODULES ON)

add_executable(test_polynomial_acceleration test_polynomial_acceleration.cpp)
target_link_libraries(test_polynomial_acceleration polynomial_nttp lam::concepts)
target_compile_features(test_polynomial_acceleration PRIVATE cxx_std_23)
set_target_properties(test_polynomial_acceleration PROPERTIES CXX_SCAN_FOR_MODULES ON)


add_subdirectory(unit)
add_subdirectory(stress)

if(TARGET lam::ctbignum)
  add_executable(check_ntt_logic manual/check_ntt_logic.cpp)
  target_link_libraries(check_ntt_logic polynomial_nttp lam::concepts lam::ctbignum)
  target_compile_features(check_ntt_logic PRIVATE cxx_std_23)
  set_target_properties(check_ntt_logic PROPERTIES CXX_SCAN_FOR_MODULES ON)
endif()
