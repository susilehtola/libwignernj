# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
cmake_minimum_required(VERSION 3.16)
project(wignernj VERSION 0.1.0 LANGUAGES C)

option(BUILD_SHARED_LIBS  "Build shared library"         ON)
option(BUILD_FORTRAN      "Build Fortran interface"      ON)
option(BUILD_PYTHON       "Build Python extension"       OFF)
option(BUILD_MPFR         "Build MPFR arbitrary-precision interface" OFF)
option(BUILD_TESTS        "Build test suite"             ON)
option(BUILD_CXX_TESTS    "Build C++ header tests"      ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)

# Core C library
set(LIB_SOURCES
    src/primes.c
    src/bigint.c
    src/pfrac.c
    src/wigner_exact.c
    src/wigner3j.c
    src/wigner6j.c
    src/wigner9j.c
    src/clebsch.c
    src/racah.c
    src/gaunt.c
)

add_library(wignernj ${LIB_SOURCES})
set_target_properties(wignernj PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
)
target_include_directories(wignernj
    PUBLIC  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
)
if(NOT MSVC)
    target_compile_options(wignernj PRIVATE -Wall -Wextra -O2)
else()
    target_compile_options(wignernj PRIVATE /W3 /O2)
endif()
target_link_libraries(wignernj PRIVATE m)

# MPFR arbitrary-precision interface (optional)
if(BUILD_MPFR)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
        target_link_libraries(wignernj PUBLIC PkgConfig::MPFR)
    else()
        find_library(MPFR_LIB NAMES mpfr REQUIRED)
        find_path   (MPFR_INC mpfr.h     REQUIRED)
        target_include_directories(wignernj PUBLIC ${MPFR_INC})
        target_link_libraries     (wignernj PUBLIC ${MPFR_LIB})
    endif()
    target_compile_definitions(wignernj PUBLIC WIGNER_HAVE_MPFR)
    install(FILES include/wigner_mpfr.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# Fortran interface (optional, only if compiler found)
if(BUILD_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        add_library(wignernj_fortran src/fortran/wigner_f90.F90)
        target_link_libraries(wignernj_fortran PUBLIC wignernj)
        set_target_properties(wignernj_fortran PROPERTIES
            Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
        )
        install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules/
                DESTINATION include/wignernj/fortran)
    else()
        message(STATUS "No Fortran compiler found -- skipping Fortran interface")
    endif()
endif()

# Python extension (optional)
if(BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        Python3_add_library(_wigner MODULE
            src/python/wignermodule.c
            ${LIB_SOURCES}
        )
        target_include_directories(_wigner PRIVATE src include)
        set_target_properties(_wigner PROPERTIES
            C_STANDARD 99
            C_STANDARD_REQUIRED ON
            OUTPUT_NAME "_wigner"
            LIBRARY_OUTPUT_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/wigner
        )
        configure_file(wigner/__init__.py
            ${CMAKE_CURRENT_BINARY_DIR}/wigner/__init__.py COPYONLY)
    else()
        message(STATUS "Python3 not found -- skipping Python extension")
    endif()
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Installation
install(TARGETS wignernj
    EXPORT  wignernjTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(FILES include/wigner.h   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES include/wigner.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(libwignernj.pc.in libwignernj.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libwignernj.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

install(EXPORT wignernjTargets
    FILE      wignernjTargets.cmake
    NAMESPACE wignernj::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wignernj)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    wignernjConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/wignernjConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wignernj
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/wignernjConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/wignernjConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/wignernjConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/wignernj
)
