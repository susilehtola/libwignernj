# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# Downstream-consumer test for the libwignernj CMake package config.
#
# This is a *separate* CMake project, not a subdirectory of the main build.
# It is configured by tests/cmake_downstream/run.sh (and by the CI) after
# libwignernj has been built and installed to a temporary prefix, and it
# exercises three things:
#
#   1. find_package(wignernj REQUIRED) imports wignernj::wignernj.
#   2. Linking a C    consumer against wignernj::wignernj works.
#   3. Linking a C++  consumer against wignernj::wignernj works
#      (the wigner.hpp header is a header-only wrapper around the C lib).
#   4. With -DWITH_FORTRAN=ON, find_package(wignernj REQUIRED COMPONENTS
#      Fortran) imports wignernj::wignernj_fortran, propagates the Fortran
#      .mod-file include path, and a Fortran consumer compiles and runs.

cmake_minimum_required(VERSION 3.16)

option(WITH_FORTRAN "Also build the Fortran consumer test" OFF)

if(WITH_FORTRAN)
    project(wignernj_downstream_test LANGUAGES C CXX Fortran)
    find_package(wignernj REQUIRED COMPONENTS Fortran)
    if(NOT wignernj_Fortran_FOUND)
        message(FATAL_ERROR
            "wignernj_Fortran_FOUND is FALSE despite REQUIRED COMPONENTS Fortran")
    endif()
    if(NOT TARGET wignernj::wignernj_fortran)
        message(FATAL_ERROR
            "wignernj::wignernj_fortran target was not imported")
    endif()
else()
    project(wignernj_downstream_test LANGUAGES C CXX)
    find_package(wignernj REQUIRED)
endif()

# C consumer
add_executable(test_c test_c.c)
target_link_libraries(test_c PRIVATE wignernj::wignernj)

# C++ consumer (header-only wrapper, still links the C library)
add_executable(test_cpp test_cpp.cpp)
set_target_properties(test_cpp PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)
target_link_libraries(test_cpp PRIVATE wignernj::wignernj)

# Fortran consumer (only when requested)
if(WITH_FORTRAN)
    add_executable(test_fort test_fort.f90)
    target_link_libraries(test_fort PRIVATE wignernj::wignernj_fortran)
endif()

# Run them through CTest so the CI surfaces failures as test failures.
enable_testing()
add_test(NAME downstream_c   COMMAND test_c)
add_test(NAME downstream_cpp COMMAND test_cpp)
if(WITH_FORTRAN)
    add_test(NAME downstream_fortran COMMAND test_fort)
endif()
