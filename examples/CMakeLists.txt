# Helper function is now in cmake/PolynomialHelpers.cmake and included by root CMakeLists.txt

# ============================================================
# Basic Examples
# ============================================================

add_polynomial_executable(type_erased)
add_polynomial_executable(audio_filtering_demo LIBS polynomial_nttp)

# Embed Wav Demo needs dependency on audio filtering to get clean wav
add_polynomial_executable(embed_wav_demo LIBS polynomial_nttp)
add_dependencies(embed_wav_demo audio_filtering_demo)

# ============================================================
# Orthogonal Polynomial Examples
# ============================================================

add_polynomial_executable(laguerre SOURCE orthogonal/laguerre.cpp)
add_polynomial_executable(hermite  SOURCE orthogonal/hermite.cpp)
add_polynomial_executable(chebyshev SOURCE orthogonal/chebyshev.cpp)
add_polynomial_executable(legendre SOURCE orthogonal/legendre.cpp)
add_polynomial_executable(gegenbauer SOURCE orthogonal/gegenbauer.cpp)
add_polynomial_executable(jacobi SOURCE orthogonal/jacobi.cpp)

# ============================================================
# Optional Boost Integration
# ============================================================

find_package(Boost QUIET)
if(Boost_FOUND)
  message(STATUS "Boost found: enabling boost::math comparison")
  
  set(BOOST_EXAMPLES laguerre hermite chebyshev legendre gegenbauer jacobi)
  foreach(example ${BOOST_EXAMPLES})
    target_link_libraries(${example} PRIVATE Boost::headers)
    target_compile_definitions(${example} PRIVATE HAS_BOOST_MATH)
  endforeach()
else()
  message(STATUS "Boost NOT found: using local reference implementation")
endif()

find_package(TBB QUIET)
if(TBB_FOUND)
  message(STATUS "TBB found: enabling parallel evaluation")
  # TBB mainly used for bulk benchmarks which moved to benchmarks dir
endif()


