# Copyright (c) David C. Anastasiu & Gheorghi Guzun
# All rights reserved.
#
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree.

# file(GLOB SNNLIB_DEMOS CONFIGURE_DEPENDS
#   ${CMAKE_CURRENT_LIST_DIR}/*.cpp
# )

SET(SNNLIB_DEMOS
  cosknng.cpp
)

foreach(program ${SNNLIB_DEMOS})
  get_filename_component(program_name ${program} NAME_WE )
  message(STATUS "  Setting up demo: ${program_name}")
  add_executable(${program_name} ${program})
  set_target_properties(${program_name} PROPERTIES
    CXX_STANDARD 20
  )
  target_link_libraries(${program_name} PRIVATE
    OpenMP::OpenMP_CXX
    snnlib
  )
endforeach()