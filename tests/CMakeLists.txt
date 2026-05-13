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

# tests/test_bigint.c probes schoolbook-internal layout
# (uint64_t-word arrays, round/sticky bit extraction); it is skipped
# when the FLINT backend is active because its bigint_t has a
# different shape and the rounding path goes through MPFR rather
# than through the bit extractor.  Every other test exercises only
# the public C API and runs with both backends.
set(C_UNIT_TESTS
    test_primes
    test_pfrac
    test_3j
    test_6j
    test_9j
    test_derived
    test_gaunt_real
    test_precisions
    test_real_ylm_in_complex_ylm
    test_symmetry
)
if(NOT BUILD_FLINT)
    list(APPEND C_UNIT_TESTS test_bigint)
endif()
foreach(test_name ${C_UNIT_TESTS})
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

# ── Warmup / zero-allocation verification ─────────────────────────────────────
# Wraps glibc malloc/calloc/realloc via dlsym(RTLD_NEXT, ...) and asserts
# that wignernj_warmup_to() makes every subsequent symbol evaluation
# allocation-free.  Built only on Linux (glibc + dlsym semantics required).
#
# Skipped under:
#   - BUILD_FLINT: FLINT's internal fmpz/mpz arithmetic allocates from
#     glibc on its own; those allocations are outside the libwignernj
#     scratch cache and would cause spurious failures.
#   - sanitizer builds: AddressSanitizer intercepts malloc itself, and
#     our dlsym(RTLD_NEXT, "malloc") overlay conflicts with that
#     interception (the test aborts in __asan_init).
string(REGEX MATCH "fsanitize" _wigner_sanitizer_flags
       "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux"
   AND NOT BUILD_FLINT
   AND NOT _wigner_sanitizer_flags)
    add_executable(test_warmup test_warmup.c)
    target_link_libraries(test_warmup PRIVATE wignernj ${LIBM_LINK} dl)
    target_include_directories(test_warmup PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME test_warmup COMMAND test_warmup)
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

# ── libquadmath tests ─────────────────────────────────────────────────────────

if(BUILD_QUADMATH)
    add_executable(test_quadmath test_quadmath.c)
    target_link_libraries(test_quadmath PRIVATE wignernj quadmath ${LIBM_LINK})
    target_include_directories(test_quadmath PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME test_quadmath COMMAND test_quadmath)
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
    add_executable(test_wignernj_fortran
        fortran/test_wignernj_fortran.f90
    )
    target_link_libraries(test_wignernj_fortran PRIVATE wignernj_f03 wignernj ${LIBM_LINK})
    target_include_directories(test_wignernj_fortran PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    set_target_properties(test_wignernj_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    add_test(NAME test_wignernj_fortran COMMAND test_wignernj_fortran)

    if(BUILD_QUADMATH)
        add_executable(test_wignernj_quadmath
            fortran/test_wignernj_quadmath.F90
        )
        target_link_libraries(test_wignernj_quadmath
            PRIVATE wignernj_f03 wignernj quadmath ${LIBM_LINK})
        target_compile_definitions(test_wignernj_quadmath
            PRIVATE WIGNERNJ_HAVE_QUADMATH)
        target_include_directories(test_wignernj_quadmath PRIVATE
            ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
        )
        set_target_properties(test_wignernj_quadmath PROPERTIES
            Fortran_MODULE_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
        )
        add_test(NAME test_wignernj_quadmath COMMAND test_wignernj_quadmath)
    endif()
endif()
