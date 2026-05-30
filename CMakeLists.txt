
cmake_minimum_required(VERSION 3.31.6 FATAL_ERROR)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.3.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "451f2fe2-a8a2-47c3-bc32-94786d8fc91b")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.0)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")
elseif(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31.6)
  set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
else()
  message(FATAL_ERROR "Version lower than 3.31.6 not supported (by me)")
endif()

set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Use compiler extensions" FORCE)
set(CMAKE_CXX_STANDARD_REQUIRED True) #is this needed?
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wno-experimental-header-units -Wno-module-import-in-extern-c)
endif()

# Single source of truth for the version: the top-level VERSION file.
# The Conan recipe reads the same file via set_version(), so the two never drift.
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _lam_polynomial_nttp_version)
string(STRIP "${_lam_polynomial_nttp_version}" _lam_polynomial_nttp_version)

project(lam_polynomial_nttp VERSION ${_lam_polynomial_nttp_version} LANGUAGES CXX)

# Re-run CMake configure whenever VERSION changes.
set_property(DIRECTORY APPEND PROPERTY
  CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")

# Support centralized cmake modules for in-tree builds
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(PolynomialHelpers)
include(LamModule)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

if(ENABLE_COVERAGE)
    message(STATUS "Coverage enabled")
    add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
    add_link_options(-fprofile-instr-generate)
endif()

include(Sanitizers)
include(BlasSetup)

option(LAM_USE_TBB "Enable TBB for bulk operations" OFF)
option(LAM_BUILD_THOUSAND_DIVISIONS "Enable the heavy compile-time stress test thousand_divisions" OFF)
if(LAM_USE_TBB)
  find_package(TBB)
  if(TBB_FOUND)
    message(STATUS "Found TBB: ${TBB_DIR} (Version ${TBB_VERSION})")
    set(LAM_USE_TBB_BOOL "true")
    add_compile_definitions(LAM_USE_TBB)
  else()
    message(STATUS "TBB not found, disabling TBB support")
    set(LAM_USE_TBB_BOOL "false")
    set(LAM_USE_TBB OFF)
  endif()
else()
  set(LAM_USE_TBB_BOOL "false")
endif()

macro(map_bool_to_string VAR RESULT)
  if(${VAR})
    set(${RESULT} "true")
  else()
    set(${RESULT} "false")
  endif()
endmacro()

map_bool_to_string(LAM_USE_ACCELERATE LAM_USE_ACCELERATE_BOOL)
map_bool_to_string(LAM_ENABLE_TSAN IS_TSAN_BUILD)
if(NOT DEFINED LAM_USE_BLAS_BOOL)
    map_bool_to_string(LAM_USE_BLAS LAM_USE_BLAS_BOOL)
endif()

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  int main() {
    unsigned __int128 a = 1;
    unsigned __int128 b = 2;
    unsigned __int128 c = a * b;
    return 0;
  }
" LAM_HAS_INT128)
map_bool_to_string(LAM_HAS_INT128 LAM_HAS_INT128_BOOL)

if(NOT TARGET lam_concepts::concepts)
  find_package(lam_concepts REQUIRED)
endif()
# BMI synthesis property dance — set on imported lam_concepts::concepts so the
# consumer's synthesizer builds the concepts modules with std + strict C++23.
# Without these, the synthesized BMI compile of installed concepts.cppm fails
# with "module 'std' not found" or a -std=gnu++23/c++23 mismatch.
if(TARGET lam_concepts::concepts)
  set_target_properties(lam_concepts::concepts PROPERTIES
    CXX_MODULE_STD ON
    CXX_EXTENSIONS OFF)
endif()

add_library(polynomial_nttp STATIC)
add_library(lam_polynomial_nttp::polynomial_nttp ALIAS polynomial_nttp)
add_library(lam::polynomial_nttp ALIAS polynomial_nttp)
target_sources(polynomial_nttp
  PUBLIC
  FILE_SET CXX_MODULES
  TYPE CXX_MODULES
  BASE_DIRS ${PROJECT_SOURCE_DIR}/src
  FILES src/polynomial_nttp.cppm
        src/polynomial_nttp-univariate-structure.cppm
        src/polynomial_nttp-univariate-algebra.cppm
        src/polynomial_nttp-univariate-finite-field.cppm
        src/polynomial_nttp-univariate-roots.cppm
        src/polynomial_nttp-univariate-berlekamp.cppm
        src/polynomial_nttp-univariate-interpolation.cppm
        src/polynomial_nttp-univariate-print.cppm
        src/polynomial_nttp-univariate-compatibility.cppm
        src/polynomial_nttp-univariate.cppm
        src/polynomial_nttp-univariate-math.cppm
        src/polynomial_nttp-univariate-fft.cppm
        src/polynomial_nttp-univariate-ntt.cppm
        src/polynomial_nttp-univariate-acceleration.cppm
)
target_sources(polynomial_nttp PRIVATE src/polynomial_nttp-univariate-fft.cpp)
target_link_libraries(polynomial_nttp PUBLIC lam_concepts::concepts)

lam_link_blas(polynomial_nttp)

if(LAM_USE_TBB)
  target_link_libraries(polynomial_nttp PRIVATE TBB::tbb)
endif()

target_compile_features(polynomial_nttp PUBLIC cxx_std_23)
set_target_properties(polynomial_nttp PROPERTIES CXX_SCAN_FOR_MODULES ON)
set_target_properties(polynomial_nttp PROPERTIES OUTPUT_NAME "lam_polynomial_nttp")

# Generate + attach the lam.polynomial_nttp:config partition. EXTRA_DEFINES feed
# the template's BLAS / Accelerate / TBB / TSAN / int128 fields (computed above).
lam_add_config_partition(polynomial_nttp NAME polynomial_nttp
  EXTRA_DEFINES
    -DLAM_USE_ACCELERATE_BOOL=${LAM_USE_ACCELERATE_BOOL}
    -DLAM_ACCELERATE_BACKEND=${LAM_ACCELERATE_BACKEND}
    -DLAM_USE_BLAS_BOOL=${LAM_USE_BLAS_BOOL}
    -DLAM_BLAS_BACKEND=${LAM_BLAS_BACKEND}
    -DLAM_USE_TBB_BOOL=${LAM_USE_TBB_BOOL}
    -DIS_TSAN_BUILD=${IS_TSAN_BUILD}
    -DLAM_HAS_INT128_BOOL=${LAM_HAS_INT128_BOOL})

add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(benchmarks)

# On macOS with -DCMAKE_INSTALL_PREFIX=/usr/local:
#   Libraries     -> /usr/local/lib/liblam_polynomial_nttp.a
#   CMake config  -> /usr/local/lib/cmake/lam_polynomial_nttp/
#   C++ modules   -> /usr/local/include/lam/polynomial_nttp/
#   CPS metadata  -> /usr/local/lib/cps/
lam_install_module_package(polynomial_nttp
  NAME        polynomial_nttp
  DESCRIPTION "lam.polynomial_nttp — number-theoretic transform polynomial arithmetic (C++23 modules)"
  LICENSE     "AGPL-3.0-or-later"
  REQUIRES    lam_concepts
)
