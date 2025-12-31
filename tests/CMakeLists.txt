include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(test_polynomial_nttp test_polynomial_nttp.cpp)

target_include_directories(test_polynomial_nttp PRIVATE ${CMAKE_BINARY_DIR})
target_link_libraries(test_polynomial_nttp polynomial_nttp lam::experimental_concepts)
target_compile_features(test_polynomial_nttp PRIVATE cxx_std_23)
set_target_properties(test_polynomial_nttp PROPERTIES CXX_SCAN_FOR_MODULES ON)
# this is for thousand_divisions
#target_compile_options(thousand_divisions PRIVATE
#    # Example for GCC/Clang
#    $<$<CXX_COMPILER_ID:GNU>:-fconstexpr-steps=1000000>
#    $<$<CXX_COMPILER_ID:Clang>:-fconstexpr-steps=1000000>
#    # Example for MSVC
#    $<$<CXX_COMPILER_ID:MSVC>:/constexpr:steps1000000>
#)

add_subdirectory(unit)
