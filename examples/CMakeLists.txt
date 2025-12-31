include_directories(${PROJECT_SOURCE_DIR}/src)
#find_package(lam_ctbignum QUIET)
add_executable(type_erased type_erased.cpp)

target_include_directories(type_erased PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(type_erased lam::polynomial_nttp)
target_compile_features(type_erased PRIVATE cxx_std_23)

# ctbignum_poly_example
#if(ctbignum_FOUND)
#  add_executable(ctbignum_poly ctbignum_poly.cpp)
#  target_include_directories(ctbignum_poly PRIVATE ${CMAKE_BINARY_DIR})
#  target_link_libraries(ctbignum_poly lam::polynomial_nttp lam::ctbignum)
#  target_compile_features(ctbignum_poly PRIVATE cxx_std_23)
#endif()
