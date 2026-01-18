# Helper functions for polynomial_nttp project

# Define a function to add a polynomial_nttp executable
function(add_polynomial_executable NAME)
  cmake_parse_arguments(ARG "" "SOURCE" "LIBS;MODULE_FILES" ${ARGN})
  
  # Use NAME.cpp as default source if not specified
  if(NOT ARG_SOURCE)
    set(ARG_SOURCE "${NAME}.cpp")
  endif()
  
  add_executable(${NAME} ${ARG_SOURCE})
  
  # Add module files if specified
  if(ARG_MODULE_FILES)
    target_sources(${NAME} PUBLIC FILE_SET CXX_MODULES FILES ${ARG_MODULE_FILES})
  endif()
  
  target_include_directories(${NAME} PRIVATE ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}/src)
  target_link_libraries(${NAME} PRIVATE lam_polynomial_nttp::polynomial_nttp ${ARG_LIBS})
  target_compile_features(${NAME} PRIVATE cxx_std_23)
  set_target_properties(${NAME} PROPERTIES CXX_SCAN_FOR_MODULES ON)
endfunction()

# Define a function to add a polynomial_nttp test
function(add_polynomial_test NAME)
  cmake_parse_arguments(ARG "" "SOURCE;TEST_NAME" "LIBS;MODULE_FILES" ${ARGN})
  
  # Use NAME as test name if not specified
  if(NOT ARG_TEST_NAME)
    set(ARG_TEST_NAME "Polynomial_Nttp_${NAME}")
  endif()
  
  # Forward arguments to add_polynomial_executable
  add_polynomial_executable(${NAME} 
    SOURCE ${ARG_SOURCE} 
    LIBS ${ARG_LIBS}
    MODULE_FILES ${ARG_MODULE_FILES}
  )
  
  add_test(NAME ${ARG_TEST_NAME} COMMAND ${NAME})
endfunction()
