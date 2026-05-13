# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026 Susi Lehtola
#
# Build language-binding examples and register them as ctest tests.
# Each `all_symbols.{c,cpp,f90,py}` exercises every public symbol family
# in libwignernj (wigner3j/6j/9j, clebsch_gordan, racah_w, fano_x, gaunt,
# gaunt_real, real_ylm_in_complex_ylm) and exits non-zero on any tolerance
# violation.  The `real_basis_lz.*` examples are focused demos of the
# real-basis l_z construction via the real/complex Y_lm change of basis.

if(WIN32)
    set(LIBM_LINK)
else()
    set(LIBM_LINK m)
endif()

# ── C examples ───────────────────────────────────────────────────────────
add_executable(example_c c/all_symbols.c)
target_link_libraries(example_c PRIVATE wignernj ${LIBM_LINK})
target_include_directories(example_c PRIVATE
    ${CMAKE_SOURCE_DIR}/include)
add_test(NAME example_c COMMAND example_c)

# Real-basis l_z built from the real <-> complex change of basis.
add_executable(example_real_basis_lz c/real_basis_lz.c)
target_link_libraries(example_real_basis_lz PRIVATE wignernj ${LIBM_LINK})
target_include_directories(example_real_basis_lz PRIVATE
    ${CMAKE_SOURCE_DIR}/include)
add_test(NAME example_real_basis_lz COMMAND example_real_basis_lz)

# ── C++ example ──────────────────────────────────────────────────────────
# Enable CXX if it isn't already.  tests/CMakeLists.txt enables it under
# BUILD_CXX_TESTS only when an actual CXX compiler is present, so we
# repeat the same probe here -- BUILD_EXAMPLES=ON without BUILD_TESTS is
# a valid configuration on machines that only have a C compiler.
if(BUILD_CXX_TESTS)
    include(CheckLanguage)
    check_language(CXX)
    if(CMAKE_CXX_COMPILER)
        # enable_language is idempotent; tests/CMakeLists.txt may have
        # called it already (when BUILD_TESTS is also on), in which case
        # this is a no-op.  Calling it unconditionally is safer than
        # querying ENABLED_LANGUAGES, which has been observed to lag
        # the actual state when subdirectories run in some orders.
        enable_language(CXX)
        add_executable(example_cpp cpp/all_symbols.cpp)
        set_target_properties(example_cpp PROPERTIES
            CXX_STANDARD          11
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS        OFF)
        target_link_libraries(example_cpp PRIVATE wignernj ${LIBM_LINK})
        target_include_directories(example_cpp PRIVATE
            ${CMAKE_SOURCE_DIR}/include)
        add_test(NAME example_cpp COMMAND example_cpp)

        # Real-basis l_z C++ parallel.
        add_executable(example_real_basis_lz_cpp cpp/real_basis_lz.cpp)
        set_target_properties(example_real_basis_lz_cpp PROPERTIES
            CXX_STANDARD          11
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS        OFF)
        target_link_libraries(example_real_basis_lz_cpp
            PRIVATE wignernj ${LIBM_LINK})
        target_include_directories(example_real_basis_lz_cpp PRIVATE
            ${CMAKE_SOURCE_DIR}/include)
        add_test(NAME example_real_basis_lz_cpp
            COMMAND example_real_basis_lz_cpp)
    else()
        message(STATUS "C++ compiler not found -- skipping example_cpp")
    endif()
endif()

# ── Fortran example ──────────────────────────────────────────────────────
if(BUILD_FORTRAN)
    add_executable(example_fortran fortran/all_symbols.f90)
    target_link_libraries(example_fortran
        PRIVATE wignernj_f03 wignernj ${LIBM_LINK})
    add_test(NAME example_fortran COMMAND example_fortran)

    # Real-basis l_z Fortran parallel.
    add_executable(example_real_basis_lz_fortran fortran/real_basis_lz.f90)
    target_link_libraries(example_real_basis_lz_fortran
        PRIVATE wignernj_f03 wignernj ${LIBM_LINK})
    add_test(NAME example_real_basis_lz_fortran
        COMMAND example_real_basis_lz_fortran)
endif()

# ── Python example ───────────────────────────────────────────────────────
# Run only when the Python extension is built in-tree (BUILD_PYTHON=ON);
# point PYTHONPATH at the build directory so `import wignernj` finds
# the freshly compiled _wignernj module without needing an editable
# install.  When BUILD_PYTHON is off the Python example still ships in
# the source tree but is not exercised by ctest.
if(BUILD_PYTHON AND TARGET _wignernj)
    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME example_python
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/python/all_symbols.py)
        # The compiled _wignernj extension lands in the build root; the
        # Python source tree (`wignernj/__init__.py`) lives at the project
        # root.  Both paths must be on PYTHONPATH for `import wignernj`
        # to resolve to the just-built extension.
        set_tests_properties(example_python PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}")

        # Real-basis l_z Python parallel.
        add_test(NAME example_real_basis_lz_python
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/python/real_basis_lz.py)
        set_tests_properties(example_real_basis_lz_python PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:${CMAKE_SOURCE_DIR}")
    endif()
endif()
