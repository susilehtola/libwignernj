# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola

# libm is needed for sqrt/ldexp/etc. on every Unix-like target;
# MSVC bundles the math functions into the default C runtime, so an
# explicit -lm there would resolve to a non-existent m.lib.  Use a
# helper variable that is empty on Windows and "m" everywhere else.
if(WIN32)
    set(LIBM_LINK)
else()
    set(LIBM_LINK m)
endif()

# ── C unit tests ──────────────────────────────────────────────────────────────

foreach(test_name
    test_primes
    test_bigint
    test_pfrac
    test_3j
    test_6j
    test_9j
    test_derived
    test_gaunt_real
    test_precisions
    test_symmetry
)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE wignernj ${LIBM_LINK})
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# ── Allocation-failure injection harness ──────────────────────────────────────
# fork()-based; only built on POSIX systems.
if(NOT WIN32)
    add_executable(test_oom test_oom.c)
    target_link_libraries(test_oom PRIVATE wignernj ${LIBM_LINK})
    target_include_directories(test_oom PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME test_oom COMMAND test_oom)
endif()

# ── C++ header tests ──────────────────────────────────────────────────────────

if(BUILD_CXX_TESTS)
    include(CheckLanguage)
    check_language(CXX)
    if(CMAKE_CXX_COMPILER)
        enable_language(CXX)
        add_executable(test_cpp test_cpp.cpp)
        target_link_libraries(test_cpp PRIVATE wignernj ${LIBM_LINK})
        target_include_directories(test_cpp PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_SOURCE_DIR}/include
        )
        set_property(TARGET test_cpp PROPERTY CXX_STANDARD 11)
        add_test(NAME test_cpp COMMAND test_cpp)
    else()
        message(STATUS "No C++ compiler found -- skipping C++ tests")
    endif()
endif()

# ── MPFR tests ────────────────────────────────────────────────────────────────

if(BUILD_MPFR)
    add_executable(test_mpfr test_mpfr.c)
    target_link_libraries(test_mpfr PRIVATE wignernj)
    target_include_directories(test_mpfr PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME test_mpfr COMMAND test_mpfr)
endif()

# ── Fortran interface tests ───────────────────────────────────────────────────

if(BUILD_FORTRAN AND TARGET wignernj_f03)
    add_executable(test_wigner_fortran
        fortran/test_wigner_fortran.f90
    )
    target_link_libraries(test_wigner_fortran PRIVATE wignernj_f03 wignernj ${LIBM_LINK})
    target_include_directories(test_wigner_fortran PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    set_target_properties(test_wigner_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    add_test(NAME test_wigner_fortran COMMAND test_wigner_fortran)
endif()
