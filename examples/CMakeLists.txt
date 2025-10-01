include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(type_erased type_erased.cpp)

target_include_directories(type_erased PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(type_erased polynomial_nttp)
target_compile_features(type_erased PRIVATE cxx_std_23)
