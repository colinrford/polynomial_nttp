# Benchmarks Directory CMakeLists

# Re-enable standard benchmarks
add_polynomial_executable(trig_benchmark LIBS polynomial_nttp)
add_polynomial_executable(eval_benchmark)
add_polynomial_executable(bulk_benchmark)
add_polynomial_executable(benchmark_complex_speed)
add_polynomial_executable(benchmark_multiplication_crossover)
add_polynomial_executable(benchmark_operator_speed)

# Compile-Time Benchmarks
add_executable(benchmark_ct_naive benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_naive PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_naive PRIVATE BENCH_MODE=0 BENCH_N=256)
target_compile_options(benchmark_ct_naive PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_naive PRIVATE cxx_std_23)

add_executable(benchmark_ct_fft benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_fft PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_fft PRIVATE BENCH_MODE=1 BENCH_N=256)
target_compile_options(benchmark_ct_fft PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_fft PRIVATE cxx_std_23)

add_executable(benchmark_ct_sunspots benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_sunspots PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_sunspots PRIVATE BENCH_MODE=2 BENCH_N=16384)
# Increase steps significantly for 16k FFT: O(N log N) is fast, but consteval overhead is high.
target_compile_options(benchmark_ct_sunspots PRIVATE -fconstexpr-steps=200000000 -fbracket-depth=1024)
target_compile_features(benchmark_ct_sunspots PRIVATE cxx_std_23)

add_executable(benchmark_ct_32k benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_32k PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_32k PRIVATE BENCH_MODE=3 BENCH_N=32768)
# Double the steps for 32k: 400M
target_compile_options(benchmark_ct_32k PRIVATE -fconstexpr-steps=400000000 -fbracket-depth=1024)
target_compile_features(benchmark_ct_32k PRIVATE cxx_std_23)

# Mode 4: Division (NTTP) - N=1024
add_executable(benchmark_ct_div benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_div PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_div PRIVATE BENCH_MODE=4 BENCH_N=1024)
target_compile_options(benchmark_ct_div PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_div PRIVATE cxx_std_23)

# Mode 5: GCD (Euclidean) - N=256
add_executable(benchmark_ct_gcd benchmark_compile_time.cpp)
target_link_libraries(benchmark_ct_gcd PRIVATE polynomial_nttp)
target_compile_definitions(benchmark_ct_gcd PRIVATE BENCH_MODE=5 BENCH_N=256)
target_compile_options(benchmark_ct_gcd PRIVATE -fconstexpr-steps=100000000)
target_compile_features(benchmark_ct_gcd PRIVATE cxx_std_23)

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
endif()
