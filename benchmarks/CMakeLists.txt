# Benchmarks Directory CMakeLists

# Re-enable standard benchmarks
add_polynomial_executable(trig_benchmark LIBS polynomial_nttp)
add_polynomial_executable(eval_benchmark)
add_polynomial_executable(bulk_benchmark)
add_polynomial_executable(benchmark_complex_speed)
add_polynomial_executable(benchmark_multiplication_crossover)
add_polynomial_executable(benchmark_addition_crossover)
add_polynomial_executable(benchmark_threading_crossover)
add_polynomial_executable(benchmark_operator_speed)



# Compile-Time Benchmarks

# 1. Naive Multiplication (N=256)
add_polynomial_executable(benchmark_ct_naive)
target_compile_options(benchmark_ct_naive PRIVATE -fconstexpr-steps=100000000)

# 2. FFT Multiplication (N=256)
add_polynomial_executable(benchmark_ct_fft)
target_compile_options(benchmark_ct_fft PRIVATE -fconstexpr-steps=100000000)

# 3. Sunspots (N=16384)
add_polynomial_executable(benchmark_ct_sunspots)
target_compile_options(benchmark_ct_sunspots PRIVATE -fconstexpr-steps=200000000 -fbracket-depth=1024)

# 5. Division (N=1024)
add_polynomial_executable(benchmark_ct_div)
target_compile_options(benchmark_ct_div PRIVATE -fconstexpr-steps=100000000)

# 6. GCD (N=256)
add_polynomial_executable(benchmark_ct_gcd)
target_compile_options(benchmark_ct_gcd PRIVATE -fconstexpr-steps=100000000)

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
