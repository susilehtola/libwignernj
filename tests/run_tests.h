/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 * Minimal C99 test harness. */
#ifndef RUN_TESTS_H
#define RUN_TESTS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

static int g_tests_run    = 0;
static int g_tests_passed = 0;
static int g_tests_failed = 0;

#define TEST_ASSERT(cond) do { \
    g_tests_run++; \
    if (!(cond)) { \
        fprintf(stderr, "FAIL  %s:%d  %s\n", __FILE__, __LINE__, #cond); \
        g_tests_failed++; \
    } else { \
        g_tests_passed++; \
    } \
} while(0)

#define TEST_NEAR(a, b, tol) do { \
    g_tests_run++; \
    double _a = (double)(a), _b = (double)(b), _t = (double)(tol); \
    double _d = fabs(_a - _b); \
    double _scale = fabs(_b) > 1e-300 ? fabs(_b) : 1.0; \
    if (_d > _t * _scale) { \
        fprintf(stderr, "FAIL  %s:%d  |%s - %s| = %.3e (tol=%.3e*|val|)\n" \
                        "       got=%.17g  expected=%.17g\n", \
                __FILE__, __LINE__, #a, #b, _d/_scale, _t, _a, _b); \
        g_tests_failed++; \
    } else { \
        g_tests_passed++; \
    } \
} while(0)

/* Relative tolerance for double: 1e-13 */
#define TEST_DBL(a, b)  TEST_NEAR(a, b, 1e-13)
/* Absolute tolerance when expected value may be 0 */
#define TEST_ABS(a, b, tol) do { \
    g_tests_run++; \
    double _d = fabs((double)(a) - (double)(b)); \
    if (_d > (double)(tol)) { \
        fprintf(stderr, "FAIL  %s:%d  |%s - %s| = %.3e > %.3e\n" \
                        "       got=%.17g  expected=%.17g\n", \
                __FILE__, __LINE__, #a, #b, _d, (double)(tol), \
                (double)(a), (double)(b)); \
        g_tests_failed++; \
    } else { g_tests_passed++; } \
} while(0)

#define SUMMARY() do { \
    printf("%s: %d/%d passed", __FILE__, g_tests_passed, g_tests_run); \
    if (g_tests_failed) printf("  (%d FAILED)", g_tests_failed); \
    printf("\n"); \
    return (g_tests_failed > 0) ? 1 : 0; \
} while(0)

#endif /* RUN_TESTS_H */
