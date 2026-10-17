cmake_minimum_required(VERSION 3.16)
project(atomsim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)

add_library(physics src/hydrogen.cpp src/nucleus.cpp src/electrons.cpp src/fission.cpp)
target_include_directories(physics PUBLIC src)

add_executable(atomsim src/main.cpp)
target_link_libraries(atomsim PRIVATE physics)

enable_testing()
add_executable(test_physics tests/test_physics.cpp)
target_link_libraries(test_physics PRIVATE physics)
add_test(NAME physics COMMAND test_physics)
