cmake_minimum_required(VERSION 3.17)
project(pysnnlib)

find_package(pybind11 REQUIRED)
pybind11_add_module(pysnnlib bindings.cpp)

target_link_libraries(pysnnlib PUBLIC)

# install(TARGETS pysnnlib
#   COMPONENT python
#   LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
# )