# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
cmake_minimum_required(VERSION 3.16)
project(wignernj VERSION 0.2.0 LANGUAGES C)

option(BUILD_SHARED_LIBS  "Build shared library"         ON)
option(BUILD_FORTRAN      "Build Fortran interface"      ON)
option(BUILD_PYTHON       "Build Python extension"       OFF)
option(BUILD_MPFR         "Build MPFR arbitrary-precision interface" OFF)
option(BUILD_QUADMATH     "Build libquadmath / __float128 interface"  OFF)
option(BUILD_FLINT        "Use FLINT/GMP/MPFR for the bigint backend"  OFF)
option(BUILD_TESTS        "Build test suite"             ON)
option(BUILD_CXX_TESTS    "Build C++ header tests"      ON)
option(BUILD_LTO          "Enable link-time optimisation (LTO/IPO)"  OFF)

# Probe the toolchain for LTO/IPO support so consumers get a clear error
# rather than a confusing link failure when the option is requested on a
# compiler/linker pair that doesn't provide it.
if(BUILD_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _wignernj_ipo_supported OUTPUT _wignernj_ipo_error)
    if(NOT _wignernj_ipo_supported)
        message(FATAL_ERROR
            "BUILD_LTO=ON but link-time optimisation is not supported by "
            "this toolchain: ${_wignernj_ipo_error}")
    endif()
endif()

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

# Core C library.  The bigint backend is selected at configure time:
# the schoolbook src/bigint.c is the default; the FLINT backend
# (src/bigint_flint.c) replaces it when BUILD_FLINT=ON.
if(BUILD_FLINT)
    set(BIGINT_SOURCE src/bigint_flint.c)
else()
    set(BIGINT_SOURCE src/bigint.c)
endif()
set(LIB_SOURCES
    src/xalloc.c
    src/primes.c
    ${BIGINT_SOURCE}
    src/pfrac.c
    src/scratch.c
    src/wigner_exact.c
    src/wigner3j.c
    src/wigner6j.c
    src/wigner9j.c
    src/clebsch.c
    src/racah.c
    src/fano_x.c
    src/gaunt.c
)

add_library(wignernj ${LIB_SOURCES})
set_target_properties(wignernj PROPERTIES
    C_STANDARD 99
    C_STANDARD_REQUIRED ON
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    INTERPROCEDURAL_OPTIMIZATION ${BUILD_LTO}
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

# libquadmath / __float128 interface (optional)
#
# __float128 is a GCC/Clang/ICC extension on Linux/macOS; not provided
# by MSVC.  The runtime functions (sqrtq, ldexpq, ...) and the header
# <quadmath.h> live in the libquadmath package shipped with GCC.
# Clang supports __float128 as a builtin and links libquadmath at link
# time, but its driver does *not* add GCC's internal include directory
# to the search path, so <quadmath.h> is unreachable by default on
# Clang/Linux.  We probe in two stages:
#
#   1. Try a direct compile-and-link with <quadmath.h> + libquadmath.
#   2. If that fails, ask `gcc -print-file-name=include` for the
#      directory containing quadmath.h and retry with that on the
#      include path.  If the second probe succeeds, propagate the
#      include directory to consumers via target_include_directories.
#
# Apple Clang and MSVC have no __float128 at all and fall through to
# the FATAL_ERROR.
if(BUILD_QUADMATH)
    include(CheckCSourceCompiles)
    set(_qm_test
        "#include <quadmath.h>\nint main(void){__float128 x=sqrtq(2.0Q);return (int)x;}\n")
    set(_qm_saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES quadmath m)
    check_c_source_compiles("${_qm_test}" WIGNERNJ_QUADMATH_DIRECT)
    set(CMAKE_REQUIRED_LIBRARIES "${_qm_saved_libs}")

    set(WIGNERNJ_QUADMATH_INCLUDE "")
    if(NOT WIGNERNJ_QUADMATH_DIRECT)
        # Fall back to gcc's include directory (where libquadmath-devel /
        # libgcc-N-dev installs quadmath.h).
        find_program(WIGNERNJ_GCC_FOR_QUADMATH
            NAMES gcc gcc-15 gcc-14 gcc-13 gcc-12 gcc-11)
        if(WIGNERNJ_GCC_FOR_QUADMATH)
            execute_process(
                COMMAND ${WIGNERNJ_GCC_FOR_QUADMATH} -print-file-name=include
                OUTPUT_VARIABLE _qm_gcc_inc
                OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(_qm_gcc_inc AND EXISTS "${_qm_gcc_inc}/quadmath.h")
                set(_qm_saved_inc "${CMAKE_REQUIRED_INCLUDES}")
                set(_qm_saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
                set(CMAKE_REQUIRED_INCLUDES  ${_qm_gcc_inc})
                set(CMAKE_REQUIRED_LIBRARIES quadmath m)
                check_c_source_compiles("${_qm_test}" WIGNERNJ_QUADMATH_GCC_INC)
                set(CMAKE_REQUIRED_INCLUDES  "${_qm_saved_inc}")
                set(CMAKE_REQUIRED_LIBRARIES "${_qm_saved_libs}")
                if(WIGNERNJ_QUADMATH_GCC_INC)
                    set(WIGNERNJ_QUADMATH_INCLUDE "${_qm_gcc_inc}")
                endif()
            endif()
        endif()
    endif()

    if(NOT WIGNERNJ_QUADMATH_DIRECT AND NOT WIGNERNJ_QUADMATH_GCC_INC)
        message(FATAL_ERROR
            "BUILD_QUADMATH=ON, but neither <quadmath.h> nor the gcc-shipped "
            "fallback path was usable.  __float128 / libquadmath are "
            "available only on GCC, Clang (with libquadmath-devel / "
            "libgcc-N-dev installed), and Intel compilers on Linux/macOS.  "
            "Disable BUILD_QUADMATH on unsupported toolchains.")
    endif()

    target_compile_definitions(wignernj PUBLIC WIGNER_HAVE_QUADMATH)
    target_link_libraries     (wignernj PUBLIC quadmath)
    if(WIGNERNJ_QUADMATH_INCLUDE)
        target_include_directories(wignernj PUBLIC ${WIGNERNJ_QUADMATH_INCLUDE})
    endif()
    install(FILES include/wigner_quadmath.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
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

# FLINT bigint backend (optional)
#
# FLINT pulls in GMP and MPFR transitively, so when BUILD_FLINT=ON we
# unconditionally link all three.  Floating-point conversions in the
# FLINT path go through MPFR (correct round-to-nearest-even at every
# IEEE 754 binary precision), and the binary128 conversion extracts
# the 113-bit mantissa via mpfr_get_z_2exp (so no MPFR --enable-float128
# build option is required).  Detection uses find_library/find_path
# directly:  FLINT does not ship a pkg-config file on every distribution
# (notably Ubuntu's libflint-dev), so the explicit search is more robust
# than pkg_check_modules.
if(BUILD_FLINT)
    find_library(FLINT_LIB NAMES flint REQUIRED)
    find_path   (FLINT_INC flint/fmpz.h REQUIRED)
    find_library(MPFR_LIB  NAMES mpfr REQUIRED)
    find_path   (MPFR_INC  mpfr.h REQUIRED)
    find_library(GMP_LIB   NAMES gmp REQUIRED)
    find_path   (GMP_INC   gmp.h REQUIRED)
    target_include_directories(wignernj PRIVATE
        ${FLINT_INC} ${MPFR_INC} ${GMP_INC})
    target_link_libraries(wignernj PRIVATE
        ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
    target_compile_definitions(wignernj PRIVATE WIGNERNJ_USE_FLINT)
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
        # Forward the quadmath gate to the Fortran preprocessor so the
        # corresponding interface block and convenience wrappers compile
        # in.  iso_fortran_env::real128 is also gated on this define
        # because, on toolchains where __float128 is unavailable, real128
        # may be -1 (no kind), which would render the module unparseable.
        if(BUILD_QUADMATH)
            target_compile_definitions(wignernj_f03 PRIVATE WIGNERNJ_HAVE_QUADMATH)
        endif()
        set_target_properties(wignernj_f03 PROPERTIES
            Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules
            VERSION   ${PROJECT_VERSION}
            SOVERSION ${PROJECT_VERSION_MAJOR}
            INTERPROCEDURAL_OPTIMIZATION ${BUILD_LTO}
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
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )
        install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fortran_modules/
                DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/wignernj/fortran)

        # pkg-config file for the Fortran wrapper -- exposes the
        # .mod-file directory via -I and depends transitively on the C
        # library's .pc file via Requires: libwignernj.
        configure_file(libwignernj_f03.pc.in libwignernj_f03.pc @ONLY)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libwignernj_f03.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

        set(WIGNERNJ_HAS_FORTRAN ON)
    else()
        message(STATUS "No Fortran compiler found -- skipping Fortran interface")
    endif()
endif()

# Python extension (optional)
#
# The Python extension recompiles the C sources directly into the
# extension module rather than linking against libwignernj, so the
# resulting `_wigner.so` is a self-contained shared object.  The
# pip-installable wheel built from setup.py is the canonical Python
# distribution path and uses the in-tree schoolbook bigint with no
# system-library dependencies.
#
# The CMake-driven Python build (BUILD_PYTHON=ON) honours the
# BUILD_FLINT switch and propagates the FLINT bigint backend
# (together with its FLINT/MPFR/GMP runtime dependencies) into
# `_wigner.so`.  This is purely an internal speedup: the Python API
# is unchanged, but the underlying bigint multiplication uses
# FLINT's sub-quadratic algorithms, which is observable as a faster
# wall-clock at very large angular momenta.  BUILD_QUADMATH and
# BUILD_MPFR are not propagated because the Python wrapper does not
# expose `__float128` or `mpfr_t` return types -- exposing those
# would require additional CPython glue (e.g. a string return or a
# gmpy2 bridge) and is left for a future release.
if(BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development)
    if(Python3_FOUND)
        Python3_add_library(_wigner MODULE
            src/python/wignermodule.c
            ${LIB_SOURCES}
        )
        target_include_directories(_wigner PRIVATE src include)
        if(BUILD_FLINT)
            target_compile_definitions(_wigner PRIVATE WIGNERNJ_USE_FLINT)
            target_include_directories(_wigner PRIVATE
                ${FLINT_INC} ${MPFR_INC} ${GMP_INC})
            target_link_libraries(_wigner PRIVATE
                ${FLINT_LIB} ${MPFR_LIB} ${GMP_LIB})
        endif()
        set_target_properties(_wigner PROPERTIES
            C_STANDARD 99
            C_STANDARD_REQUIRED ON
            INTERPROCEDURAL_OPTIMIZATION ${BUILD_LTO}
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
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES include/wigner.h   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES include/wigner.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# pkg-config file: libm only on Unix-like targets (MSVC bundles math
# into the default C runtime, so an explicit -lm there would resolve
# to a non-existent m.lib).
if(WIN32)
    set(WIGNERNJ_PC_LIBM)
else()
    set(WIGNERNJ_PC_LIBM "-lm")
endif()
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
