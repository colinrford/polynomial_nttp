# Benchmarks Directory CMakeLists

# Re-enable standard benchmarks
add_polynomial_executable(trig_benchmark LIBS polynomial_nttp)
add_polynomial_executable(eval_benchmark)
add_polynomial_executable(bulk_benchmark)
add_polynomial_executable(benchmark_complex_speed)
add_polynomial_executable(benchmark_multiplication_crossover)
add_polynomial_executable(benchmark_addition_crossover)
add_polynomial_executable(benchmark_threading_crossover)

# Find ctbignum for benchmarks
find_package(lam_ctbignum QUIET)
if(TARGET lam::ctbignum)
  get_target_property(_aliased lam::ctbignum ALIASED_TARGET)
  if(_aliased)
    set(_ctbignum_target ${_aliased})
  else()
    set(_ctbignum_target lam::ctbignum)
  endif()
  set_property(TARGET ${_ctbignum_target} APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_23)

  add_polynomial_executable(benchmark_operator_speed)
  add_polynomial_executable(benchmark_ntt_vs_naive LIBS lam::ctbignum)
  add_polynomial_executable(benchmark_ntt_pure LIBS lam::ctbignum)

  if(TARGET lam::interop)
    add_polynomial_executable(benchmark_parallel_ntt LIBS lam::ctbignum lam::interop)
  endif()

  # NTL Benchmark (if available)
  find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h PATHS /usr/local/include /opt/homebrew/include)
  find_library(NTL_LIBRARY NAMES ntl PATHS /usr/local/lib /opt/homebrew/lib)

  if(NTL_INCLUDE_DIR AND NTL_LIBRARY)
    message(STATUS "NTL found: enabling NTL benchmark")
    add_polynomial_executable(benchmark_ntt_vs_ntl
      SOURCE benchmark_ntt_vs_ntl.cpp ntl_benchmark_impl.cpp
      LIBS lam::ctbignum
    )
    target_include_directories(benchmark_ntt_vs_ntl PRIVATE ${NTL_INCLUDE_DIR})
    target_link_libraries(benchmark_ntt_vs_ntl PRIVATE ${NTL_LIBRARY})
  else()
    message(STATUS "NTL not found")
  endif()

  add_polynomial_executable(benchmark_ntt_vector_vs_array
    SOURCE benchmark_ntt_vector_vs_array.cpp
    LIBS lam::ctbignum
  )
endif()

# Compile-Time Benchmarks

# 1. Naive Multiplication (N=256)
add_executable(benchmark_ct_naive benchmark_ct_naive.cpp)
target_link_libraries(benchmark_ct_naive PRIVATE polynomial_nttp)
target_compile_options(benchmark_ct_naive PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_naive PRIVATE cxx_std_23)
set_target_properties(benchmark_ct_naive PROPERTIES CXX_SCAN_FOR_MODULES ON)

# 2. FFT Multiplication (N=256)
add_executable(benchmark_ct_fft benchmark_ct_fft.cpp)
target_link_libraries(benchmark_ct_fft PRIVATE polynomial_nttp)
target_compile_options(benchmark_ct_fft PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_fft PRIVATE cxx_std_23)
set_target_properties(benchmark_ct_fft PROPERTIES CXX_SCAN_FOR_MODULES ON)

# 3. Sunspots (N=16384)
add_executable(benchmark_ct_sunspots benchmark_ct_sunspots.cpp)
target_link_libraries(benchmark_ct_sunspots PRIVATE polynomial_nttp)
target_compile_options(benchmark_ct_sunspots PRIVATE -fconstexpr-steps=200000000 -fbracket-depth=1024)
target_compile_features(benchmark_ct_sunspots PRIVATE cxx_std_23)
set_target_properties(benchmark_ct_sunspots PROPERTIES CXX_SCAN_FOR_MODULES ON)



# 5. Division (N=1024)
add_executable(benchmark_ct_div benchmark_ct_div.cpp)
target_link_libraries(benchmark_ct_div PRIVATE polynomial_nttp)
target_compile_options(benchmark_ct_div PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_div PRIVATE cxx_std_23)
set_target_properties(benchmark_ct_div PROPERTIES CXX_SCAN_FOR_MODULES ON)

# 6. GCD (N=256)
add_executable(benchmark_ct_gcd benchmark_ct_gcd.cpp)
target_link_libraries(benchmark_ct_gcd PRIVATE polynomial_nttp)
target_compile_options(benchmark_ct_gcd PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_gcd PRIVATE cxx_std_23)
set_target_properties(benchmark_ct_gcd PROPERTIES CXX_SCAN_FOR_MODULES ON)

# Enable module scanning for standard benchmarks using modules
foreach(t 
    benchmark_multiplication_crossover
)
  if(TARGET ${t})
    set_target_properties(${t} PROPERTIES CXX_SCAN_FOR_MODULES ON)
  endif()
endforeach()

# Boost Integration for Benchmarks
if(Boost_FOUND)
  target_link_libraries(eval_benchmark PRIVATE Boost::headers)
  target_compile_definitions(eval_benchmark PRIVATE HAS_BOOST_MATH)
  
  target_link_libraries(bulk_benchmark PRIVATE Boost::headers)
  target_compile_definitions(bulk_benchmark PRIVATE HAS_BOOST_MATH)
endif()

# TBB Integration for Benchmarks
if(TBB_FOUND)
  target_link_libraries(bulk_benchmark PRIVATE TBB::tbb)
  target_compile_definitions(bulk_benchmark PRIVATE LAM_USE_TBB)
  
  target_link_libraries(benchmark_threading_crossover PRIVATE TBB::tbb)
  target_compile_definitions(benchmark_threading_crossover PRIVATE LAM_USE_TBB)
endif()
