cmake_minimum_required(VERSION 3.13)
project(jupyter-kernel CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory(xeus-cling-cmake-setup)

# Fiddle out potential library dirs. This will all vanish with our new
# target-based build system - hopefully!
set(libdirs)
set(incdirs)
foreach(mod ${ALL_DEPENDENCIES})
  if(EXISTS ${${mod}_DIR})
    set(libdirs ${libdirs} ${${mod}_DIR}/lib)
  endif()
  if(EXISTS ${${mod}_INCLUDE_DIRS})
    set(incdirs ${incdirs} ${${mod}_INCLUDE_DIRS})
  endif()
endforeach()

xeus_cling_setup(
  INCLUDE_DIRECTORIES ${incdirs} ${CMAKE_SOURCE_DIR}
  LINK_LIBRARIES ${DUNE_LIBS}
  LIBRARY_DIRECTORIES ${libdirs}
  SETUP_HEADERS
    ${CMAKE_BINARY_DIR}/config.h
  DOXYGEN_URLS https://dune-project.org/doxygen/2.7.1
  DOXYGEN_TAGFILES "Dune Core Modules.tag"
  KERNEL_NAME "C++17 (DUNE)"
  REQUIRED
)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
