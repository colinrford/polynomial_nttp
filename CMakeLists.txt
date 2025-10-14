
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

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED True) #is this needed?
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_SCAN_FOR_MODULES ON)

project(math_nttp_polynomial LANGUAGES CXX)

enable_testing()

configure_file(
  config.cppm.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.cppm
)

add_library(experimental_concepts)
target_sources(experimental_concepts PUBLIC
  FILE_SET CXX_MODULES
  TYPE CXX_MODULES
  FILES src/experimental/concepts.cppm
)
target_compile_features(experimental_concepts PUBLIC cxx_std_23)

add_library(polynomial_nttp)

target_sources(polynomial_nttp
  PUBLIC
  FILE_SET CXX_MODULES
  TYPE CXX_MODULES
  FILES src/polynomial_nttp.cppm
        src/polynomial_nttp-univariate-structure.cppm
        src/polynomial_nttp-univariate-algebra.cppm
        src/polynomial_nttp-univariate.cppm
)

target_link_libraries(polynomial_nttp experimental_concepts)
target_compile_features(polynomial_nttp PUBLIC cxx_std_23)

add_subdirectory(tests)
#add_subdirectory(examples)
