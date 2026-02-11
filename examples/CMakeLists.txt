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

find_package(lam_ctbignum QUIET)
if(TARGET lam::ctbignum)
  get_target_property(_aliased lam::ctbignum ALIASED_TARGET)
  if(_aliased)
    set(_ctbignum_target ${_aliased})
  else()
    set(_ctbignum_target lam::ctbignum)
  endif()
  set_property(TARGET ${_ctbignum_target} APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_23)

  # Use helper for compatibility examples
  add_polynomial_executable(ctbignum_compatibility 
    SOURCE compatibility/ctbignum_compat.cpp 
    LIBS lam::ctbignum
  )

  target_compile_options(ctbignum_compatibility PRIVATE -fconstexpr-steps=150000000 -fbracket-depth=1024)

  add_polynomial_executable(roots_brute_force 
    SOURCE compatibility/roots_brute_force.cpp 
    LIBS lam::ctbignum
  )

  add_polynomial_executable(root_finding_demo 
    SOURCE root_finding_demo.cpp 
    LIBS lam::ctbignum
  )
endif()
