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

set(ATOMSIM_DATA_DIR "${CMAKE_SOURCE_DIR}/data" CACHE STRING "Where the nuclear data files live at run time")

add_library(physics src/hydrogen.cpp src/nucleus.cpp src/electrons.cpp src/fission.cpp
            src/chain.cpp)
target_include_directories(physics PUBLIC src)
if(NOT EMSCRIPTEN)  # the browser build runs single-threaded
  find_package(Threads REQUIRED)
  target_link_libraries(physics PUBLIC Threads::Threads)
endif()
target_compile_options(physics PRIVATE ${ATOMSIM_WARNINGS})
target_compile_definitions(physics PRIVATE ATOMSIM_DATA_DIR="${ATOMSIM_DATA_DIR}")

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

option(ATOMSIM_VIZ "Build the raylib visualisation" OFF)
if(ATOMSIM_VIZ)
  include(FetchContent)
  set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(raylib
    URL https://github.com/raysan5/raylib/archive/refs/tags/5.5.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SYSTEM)
  FetchContent_MakeAvailable(raylib)
  add_executable(atomsim-viz viz/main.cpp)
  target_link_libraries(atomsim-viz PRIVATE physics raylib)
  target_compile_options(atomsim-viz PRIVATE ${ATOMSIM_WARNINGS})
  if(EMSCRIPTEN)
    set_target_properties(atomsim-viz PROPERTIES OUTPUT_NAME index SUFFIX ".html")
    target_link_options(atomsim-viz PRIVATE
      -sUSE_GLFW=3 -sALLOW_MEMORY_GROWTH=1 -sGL_ENABLE_GET_PROC_ADDRESS=1
      --shell-file ${CMAKE_SOURCE_DIR}/viz/shell.html
      "--preload-file=${CMAKE_SOURCE_DIR}/data@/data")
  endif()
endif()
