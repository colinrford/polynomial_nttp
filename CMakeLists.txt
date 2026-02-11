
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
add_compile_options(-Wno-experimental-header-units -Wno-module-import-in-extern-c)

project(lam_polynomial_nttp VERSION 0.2.260210 LANGUAGES CXX)

# Support centralized cmake modules for in-tree builds
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
endif()
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

if(CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE MATCHES "^[Dd][Ee][Bb][Uu][Gg]$")
  set(IS_DEBUG "true")
else()
  set(IS_DEBUG "false")
endif()

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
  find_package(TBB REQUIRED)
  message(STATUS "Found TBB: ${TBB_DIR} (Version ${TBB_VERSION})")
  set(LAM_USE_TBB_BOOL "true")
  add_compile_definitions(LAM_USE_TBB)
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

configure_file(
  polynomial_nttp_config.cppm.in
  ${CMAKE_CURRENT_BINARY_DIR}/polynomial_nttp_config.cppm
)

if(NOT TARGET lam::concepts)
  find_package(lam_concepts REQUIRED)
endif()

add_library(polynomial_nttp STATIC)
add_library(lam_polynomial_nttp::polynomial_nttp ALIAS polynomial_nttp)
add_library(lam::polynomial_nttp ALIAS polynomial_nttp)
target_sources(polynomial_nttp
  PUBLIC
  FILE_SET CXX_MODULES
  TYPE CXX_MODULES
  BASE_DIRS ${PROJECT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}
  FILES src/polynomial_nttp.cppm
        src/polynomial_nttp-univariate-structure.cppm
        src/polynomial_nttp-univariate-algebra.cppm
        src/polynomial_nttp-univariate-roots.cppm
        src/polynomial_nttp-univariate-berlekamp.cppm
        src/polynomial_nttp-univariate-print.cppm
        src/polynomial_nttp-univariate-compatibility.cppm
        src/polynomial_nttp-univariate.cppm
        src/polynomial_nttp-univariate-math.cppm 
        src/polynomial_nttp-univariate-fft.cppm
        src/polynomial_nttp-univariate-ntt.cppm
        src/polynomial_nttp-univariate-acceleration.cppm
        ${CMAKE_CURRENT_BINARY_DIR}/polynomial_nttp_config.cppm
)
target_sources(polynomial_nttp PRIVATE src/polynomial_nttp-univariate-fft.cpp)
target_link_libraries(polynomial_nttp PUBLIC lam::concepts)

lam_link_blas(polynomial_nttp)

if(LAM_USE_TBB)
  target_link_libraries(polynomial_nttp PRIVATE TBB::tbb)
endif()

target_compile_features(polynomial_nttp PUBLIC cxx_std_23)
set_target_properties(polynomial_nttp PROPERTIES CXX_SCAN_FOR_MODULES ON)
set_target_properties(polynomial_nttp PROPERTIES OUTPUT_NAME "lam_polynomial_nttp")

add_subdirectory(tests)
add_subdirectory(examples)
add_subdirectory(benchmarks)

include(GNUInstallDirs)

install(TARGETS polynomial_nttp
  EXPORT polynomial_nttpTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET CXX_MODULES DESTINATION ${CMAKE_INSTALL_LIBDIR}/c++/v1/lam/polynomial_nttp
)

install(EXPORT polynomial_nttpTargets
  FILE polynomial_nttpTargets.cmake
  NAMESPACE lam_polynomial_nttp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lam_polynomial_nttp
)

include(CMakePackageConfigHelpers)

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/lam_polynomial_nttpConfig.cmake" [[
include(CMakeFindDependencyMacro)
find_dependency(lam_concepts)
include("${CMAKE_CURRENT_LIST_DIR}/polynomial_nttpTargets.cmake")
# Backward compatibility alias
if(TARGET lam_polynomial_nttp::polynomial_nttp AND NOT TARGET lam::polynomial_nttp)
  get_target_property(_aliased lam_polynomial_nttp::polynomial_nttp ALIASED_TARGET)
  if(_aliased)
    add_library(lam::polynomial_nttp ALIAS ${_aliased})
  else()
    add_library(lam::polynomial_nttp ALIAS lam_polynomial_nttp::polynomial_nttp)
  endif()
  unset(_aliased)
endif()
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
