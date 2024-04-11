cmake_minimum_required(VERSION 3.17)

# Get the version of the library
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(snnlib)
snnlib_getversion(VERSION SOVERSION)

project(snnlib
    VERSION ${VERSION}
    DESCRIPTION "Efficient nearest neighor search and graph construction for sparse data."
    HOMEPAGE_URL "https://github.com/davidanastasiu/snnlib"
    LANGUAGES CXX)

# Include guards
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(
        FATAL_ERROR
        "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
    )
endif()

#--- Project options ---
option(SNNLIB_ENABLE_DEMOS "Build demos." OFF)
option(SNNLIB_ENABLE_PYTHON "Build Python extension." OFF)
option(SNNLIB_BUILD_TESTING "Build tests." OFF)
set(SNNLIB_BUILD_TESTING ON)
set(SNNLIB_ENABLE_DEMOS ON)
set(SNNLIB_ENABLE_PYTHON ON)

# Add targets
message(STATUS "Adding target: snnlib")
find_package(OpenMP REQUIRED)
add_subdirectory(snnlib)

if(SNNLIB_ENABLE_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)
    add_subdirectory(python_bindings)
endif()

if(SNNLIB_ENABLE_DEMOS)
    message(STATUS "Adding target: demos")
    add_subdirectory(demos)
endif()

# CTest must be included in the top level to enable the `make test` target.
if(SNNLIB_BUILD_TESTING)
    include(CTest)
    message(STATUS "Adding target: tests")
    add_subdirectory(tests)
endif()