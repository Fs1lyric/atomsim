cmake_minimum_required(VERSION 3.16)
project(atomsim C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
# Warnings apply to our own targets only, not to fetched dependencies.
option(ATOMSIM_WERROR "Treat warnings as errors" OFF)
set(ATOMSIM_WARNINGS -Wall -Wextra -Wpedantic $<$<BOOL:${ATOMSIM_WERROR}>:-Werror>)

find_package(Threads REQUIRED)

add_library(physics src/hydrogen.cpp src/nucleus.cpp src/electrons.cpp src/fission.cpp)
target_include_directories(physics PUBLIC src)
target_link_libraries(physics PUBLIC Threads::Threads)
target_compile_options(physics PRIVATE ${ATOMSIM_WARNINGS})
target_compile_definitions(physics PRIVATE ATOMSIM_DATA_DIR="${CMAKE_SOURCE_DIR}/data")

add_executable(atomsim src/main.cpp)
target_link_libraries(atomsim PRIVATE physics)
target_compile_options(atomsim PRIVATE ${ATOMSIM_WARNINGS})

option(ATOMSIM_TESTS "Build the GoogleTest suite" ON)
if(ATOMSIM_TESTS)
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  add_executable(test_physics tests/test_physics.cpp)
  target_link_libraries(test_physics PRIVATE physics GTest::gtest_main)
  target_compile_options(test_physics PRIVATE ${ATOMSIM_WARNINGS})
  include(GoogleTest)
  gtest_discover_tests(test_physics)
endif()
