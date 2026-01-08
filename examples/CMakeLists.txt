# Helper function is now in cmake/PolynomialHelpers.cmake and included by root CMakeLists.txt

# ============================================================
# Basic Examples
# ============================================================

add_polynomial_executable(type_erased)

# ============================================================
# Benchmark Examples
# ============================================================

add_polynomial_executable(trig_benchmark LIBS polynomial_nttp)
add_polynomial_executable(eval_benchmark)

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

find_package(lam_ctbignum QUIET)
if(TARGET lam::ctbignum)
    set_property(TARGET lam::ctbignum APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_23)
endif()

# Use helper for compatibility examples
add_polynomial_executable(ctbignum_compatibility 
  SOURCE compatibility/ctbignum_compat.cpp 
  LIBS lam::ctbignum
)

add_polynomial_executable(roots_brute_force 
  SOURCE compatibility/roots_brute_force.cpp 
  LIBS lam::ctbignum
)
