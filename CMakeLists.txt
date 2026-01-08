
cmake_minimum_required(VERSION 3.31.6 FATAL_ERROR)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 4.0.3)
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

project(lam_polynomial_nttp VERSION 0.2.260107 LANGUAGES CXX)

# Make helpers available to subdirectories
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(PolynomialHelpers)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
else()
  set(GIT_HASH "unknown")
endif()

enable_testing()

#find_program(CLANG_TIDY_EXE "clang-tidy")
#if(CLANG_TIDY_EXE)
#  set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE} "-checks=-*,bugprone-*")
#endif()

if(CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
  set(IS_DEBUG "true")
else()
  set(IS_DEBUG "false")
endif()

option(ENABLE_COVERAGE "Enable coverage build" OFF)
if(ENABLE_COVERAGE)
    message(STATUS "Coverage enabled")
    add_compile_options(-fprofile-instr-generate -fcoverage-mapping)
    add_link_options(-fprofile-instr-generate)
endif()

configure_file(
  config.cppm.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)

# Use installed lam_concepts library
find_package(lam_concepts REQUIRED)

add_library(polynomial_nttp STATIC)
add_library(lam::polynomial_nttp ALIAS polynomial_nttp)
target_sources(polynomial_nttp
  PUBLIC
  FILE_SET CXX_MODULES
  TYPE CXX_MODULES
  FILES src/polynomial_nttp.cppm
        src/polynomial_nttp-univariate-structure.cppm
        src/polynomial_nttp-univariate-algebra.cppm
        src/polynomial_nttp-univariate-roots.cppm
        src/polynomial_nttp-univariate-berlekamp.cppm
        src/polynomial_nttp-univariate-print.cppm
        src/polynomial_nttp-univariate-compatibility.cppm
        src/polynomial_nttp-univariate.cppm
        ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)
target_link_libraries(polynomial_nttp lam::concepts)
target_compile_features(polynomial_nttp PUBLIC cxx_std_23)
set_target_properties(polynomial_nttp PROPERTIES CXX_SCAN_FOR_MODULES ON)
set_target_properties(polynomial_nttp PROPERTIES OUTPUT_NAME "lam_polynomial_nttp")

# lam_interop library (bridges polynomial_nttp with ctbignum and bitvector)
# Only available when installed as a package
find_package(lam_interop QUIET)

if(lam_interop_FOUND)
  message(STATUS "Found lam_interop package")
  set(LAM_INTEROP_AVAILABLE TRUE CACHE BOOL "lam_interop library is available" FORCE)
else()
  message(STATUS "lam_interop not found - roots_finite_field test will be skipped")
  set(LAM_INTEROP_AVAILABLE FALSE CACHE BOOL "lam_interop library is available" FORCE)
endif()

add_subdirectory(tests)
add_subdirectory(examples)

# Installation Configuration
include(GNUInstallDirs)

# Library Installation
install(TARGETS polynomial_nttp
  EXPORT polynomial_nttpTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_LIBDIR}/c++/v1
)

install(EXPORT polynomial_nttpTargets
  FILE polynomial_nttpTargets.cmake
  NAMESPACE lam::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_polynomial_nttp
)

include(CMakePackageConfigHelpers)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/lam_polynomial_nttpConfig.cmake" [[
include(CMakeFindDependencyMacro)
find_dependency(lam_concepts)
include("${CMAKE_CURRENT_LIST_DIR}/polynomial_nttpTargets.cmake")
]])

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/lam_polynomial_nttpConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/lam_polynomial_nttpConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/lam_polynomial_nttpConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_polynomial_nttp
)
