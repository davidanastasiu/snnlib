cmake_minimum_required(VERSION 3.16)

# Add the library
add_library(snnlib INTERFACE)

# explicitly set version and soversion
set_target_properties(
    snnlib PROPERTIES
    SOVERSION ${SOVERSION}
    VERSION ${VERSION}
)

# set compile standard for C++
target_compile_features(snnlib INTERFACE cxx_std_20)

# Handle `#include <snnlib/foo.h>`.
target_include_directories(snnlib SYSTEM INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>)

find_package(OpenMP REQUIRED)

message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
target_link_libraries(snnlib INTERFACE OpenMP::OpenMP_CXX)

# installation instructions
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/snnlib-config-version.cmake
    VERSION ${VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/snnlib-config-version.cmake
    DESTINATION "lib/cmake/snnlib"
)