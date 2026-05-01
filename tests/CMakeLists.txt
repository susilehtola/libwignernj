# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola

# ── C unit tests ──────────────────────────────────────────────────────────────

foreach(test_name
    test_primes
    test_bigint
    test_pfrac
    test_3j
    test_6j
    test_9j
    test_derived
    test_precisions
)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE wignernj m)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

# ── C++ header tests ──────────────────────────────────────────────────────────

if(BUILD_CXX_TESTS)
    include(CheckLanguage)
    check_language(CXX)
    if(CMAKE_CXX_COMPILER)
        enable_language(CXX)
        add_executable(test_cpp test_cpp.cpp)
        target_link_libraries(test_cpp PRIVATE wignernj m)
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

# ── Fortran interface tests ───────────────────────────────────────────────────

if(BUILD_FORTRAN AND TARGET wignernj_fortran)
    add_executable(test_wigner_fortran
        fortran/test_wigner_fortran.f90
    )
    target_link_libraries(test_wigner_fortran PRIVATE wignernj_fortran wignernj m)
    target_include_directories(test_wigner_fortran PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    set_target_properties(test_wigner_fortran PROPERTIES
        Fortran_MODULE_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/../fortran_modules
    )
    add_test(NAME test_wigner_fortran COMMAND test_wigner_fortran)
endif()
