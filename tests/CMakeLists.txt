# Copyright (c) David C. Anastasiu
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

cmake_minimum_required(VERSION 3.17)

enable_testing()

# file(GLOB SNNLIB_TEST_SRC CONFIGURE_DEPENDS
#   ${CMAKE_CURRENT_LIST_DIR}/test_*.cpp
# )
SET(SNNLIB_TEST_SRC
  test_bit_vector.cpp
  test_kv_vector.cpp
  test_kvv_view.cpp
  test_pqueue.cpp
  test_sort.cpp
  test_csr.cpp
  test_queue.cpp
  test_cosknnindex.cpp
)
message(STATUS "SNNLib tests: ${SNNLIB_TEST_SRC}")

add_executable(snnlib_test ${SNNLIB_TEST_SRC})
set_target_properties(snnlib_test PROPERTIES
  CXX_STANDARD 20
)

include(FetchContent)
FetchContent_Declare(googletest
  GIT_REPOSITORY    https://github.com/google/googletest.git
  GIT_TAG           389cb68b87193358358ae87cc56d257fd0d80189
  SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-src"
  BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ""
  INSTALL_COMMAND   ""
  TEST_COMMAND      ""
)
set(BUILD_GMOCK CACHE BOOL OFF)
set(INSTALL_GTEST CACHE BOOL OFF)
FetchContent_MakeAvailable(googletest)

find_package(OpenMP REQUIRED)

target_link_libraries(snnlib_test PRIVATE
  OpenMP::OpenMP_CXX
  snnlib
  gtest_main
)

# Defines `gtest_discover_tests()`.
include(GoogleTest)
gtest_discover_tests(snnlib_test)
