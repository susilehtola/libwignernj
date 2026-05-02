# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
cmake_minimum_required(VERSION 3.16)
project(wignernj VERSION 0.2.0 LANGUAGES C)

option(BUILD_SHARED_LIBS  "Build shared library"         ON)
option(BUILD_FORTRAN      "Build Fortran interface"      ON)
option(BUILD_PYTHON       "Build Python extension"       OFF)
option(BUILD_MPFR         "Build MPFR arbitrary-precision interface" OFF)
option(BUILD_TESTS        "Build test suite"             ON)
option(BUILD_CXX_TESTS    "Build C++ header tests"      ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# On Windows, the dynamic linker only searches the executable's own
# directory and PATH for required DLLs.  Without intervention, CMake
# would put wignernj.dll in <build>/Release and the test executables
# in <build>/tests/Release, and ctest would fail with status 0xc0000135
# (STATUS_DLL_NOT_FOUND).  Co-locate every runtime output in a single
# directory so the in-tree tests resolve the library naturally.  The
# Linux/macOS layout is left untouched: rpath/runpath already takes
# care of in-tree resolution there, and downstream tooling such as
# benchmarks/Makefile expects libwignernj.so directly under build/.
if(WIN32)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

include(GNUInstallDirs)

# Core C library
set(LIB_SOURCES
    src/xalloc.c
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
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    # Auto-generate the MSVC import library by exporting every function
    # symbol.  On ELF/Mach-O the public-by-default visibility makes this
    # a no-op; on Windows it is what saves us from sprinkling
    # __declspec(dllexport) over every public function (or supplying a
    # hand-written .def file) to make wignernj.lib appear next to
    # wignernj.dll.  Note that data symbols (the prime tables) are *not*
    # auto-exported and are decorated explicitly via the WIGNERNJ_DATA
    # macro in src/primes.h, gated on WIGNERNJ_BUILDING_DLL below.
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)
if(BUILD_SHARED_LIBS)
    # Drive the WIGNERNJ_DATA macro alternation in src/primes.h:
    # WIGNERNJ_BUILDING_DLL only on the library's own translation units
    # (PRIVATE), WIGNERNJ_DLL also on every consumer (PUBLIC) so that
    # the data symbols carry __declspec(dllimport) where consumed.  In
    # a static build neither macro is defined and the macro stays
    # inert, which is what we want for static linking.
    target_compile_definitions(wignernj
        PRIVATE WIGNERNJ_BUILDING_DLL
        PUBLIC  WIGNERNJ_DLL)
endif()
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
if(NOT WIN32)
    # libm is needed for sqrt/ldexp/etc. on every Unix-like target;
    # MSVC bundles the math functions into the default C runtime.
    target_link_libraries(wignernj PRIVATE m)
endif()

# MPFR arbitrary-precision interface (optional)
#
# mpfr.h does `#include <gmp.h>`, but on multi-prefix packagings (notably
# Homebrew on macOS) GMP lives in its own prefix and mpfr.pc does not
# propagate the GMP include path.  We therefore pull GMP in alongside MPFR.
if(BUILD_MPFR)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
        pkg_check_modules(GMP            IMPORTED_TARGET gmp)
        target_link_libraries(wignernj PUBLIC PkgConfig::MPFR)
        if(GMP_FOUND)
            target_link_libraries(wignernj PUBLIC PkgConfig::GMP)
        endif()
    else()
        find_library(MPFR_LIB NAMES mpfr REQUIRED)
        find_path   (MPFR_INC mpfr.h     REQUIRED)
        find_path   (GMP_INC  gmp.h)
        target_include_directories(wignernj PUBLIC ${MPFR_INC})
        if(GMP_INC)
            target_include_directories(wignernj PUBLIC ${GMP_INC})
        endif()
        target_link_libraries     (wignernj PUBLIC ${MPFR_LIB})
    endif()
    target_compile_definitions(wignernj PUBLIC WIGNER_HAVE_MPFR)
    install(FILES include/wigner_mpfr.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()

# Fortran interface (optional, only if compiler found)
set(WIGNERNJ_HAS_FORTRAN OFF)
if(BUILD_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
        add_library(wignernj_f03 src/fortran/wigner_f90.F90)
        target_link_libraries(wignernj_f03 PUBLIC wignernj)
        set_target_properties(wignernj_f03 PROPERTIES
            Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
            VERSION   ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION_MAJOR}
            WINDOWS_EXPORT_ALL_SYMBOLS ON
            # Bake $ORIGIN into the installed library so that, when a
            # downstream binary's DT_RUNPATH points at some non-system
            # prefix where libwignernj_f03.so lives, the dynamic linker
            # still finds its co-located libwignernj.so.0 transitive
            # dependency.  Linkers configured with
            # --no-copy-dt-needed-entries (the Ubuntu default) do not
            # propagate libwignernj.so.0 into the consuming binary's
            # NEEDED list, and DT_RUNPATH (unlike DT_RPATH) is *not*
            # searched for transitive deps -- so without this the
            # downstream Fortran consumer fails on those distros.
            INSTALL_RPATH "$ORIGIN"
        )
        # Propagate the .mod-file directory to consumers, both at build
        # time and after install, so downstream Fortran compilations
        # find `use wigner` automatically.
        target_include_directories(wignernj_f03 INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/fortran_modules>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/wignernj/fortran>
        )
        install(TARGETS wignernj_f03
            EXPORT  wignernjTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
        install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules/
                DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wignernj/fortran)
        set(WIGNERNJ_HAS_FORTRAN ON)
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
