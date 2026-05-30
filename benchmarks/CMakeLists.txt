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
include(LamConstexprLimit)

# 1. Naive Multiplication (N=256)
add_polynomial_executable(benchmark_ct_naive)
lam_target_constexpr_limit(benchmark_ct_naive 100000000)

# 2. FFT Multiplication (N=256)
add_polynomial_executable(benchmark_ct_fft)
lam_target_constexpr_limit(benchmark_ct_fft 100000000)

# 3. Sunspots (N=16384)
add_polynomial_executable(benchmark_ct_sunspots)
lam_target_constexpr_limit(benchmark_ct_sunspots 200000000)
target_compile_options(benchmark_ct_sunspots PRIVATE "$<$<CXX_COMPILER_ID:Clang,AppleClang>:-fbracket-depth=1024>")

# 5. Division (N=1024)
add_polynomial_executable(benchmark_ct_div)
lam_target_constexpr_limit(benchmark_ct_div 100000000)

# 6. GCD (N=256)
add_polynomial_executable(benchmark_ct_gcd)
lam_target_constexpr_limit(benchmark_ct_gcd 100000000)

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
