/* SPDX-License-Identifier: BSD-3-Clause
 * Copyright (c) 2026 Susi Lehtola
 *
 * Allocation-failure injection harness.
 *
 * For each public symbol path (3j, 6j, 9j, gaunt) we force the (N+1)-th
 * allocation to fail and verify that the library aborts cleanly with
 * SIGABRT — i.e. that the diagnostic-and-abort path in xalloc.c is reached
 * rather than a NULL deref / SIGSEGV.
 *
 * Each forced failure runs in a forked child so the parent process keeps a
 * clean state and can sweep N over a range of values.  The child's stderr
 * is silenced because xmalloc/xcalloc/xrealloc print a diagnostic before
 * abort() and we don't want it cluttering ctest output.
 *
 * Skipped on Windows (no fork()) — print a SKIP line and return success.
 */

#include "run_tests.h"
#include "../include/wigner.h"

#if defined(_WIN32)

int main(void)
{
    printf("test_oom: SKIP (no fork() on Windows)\n");
    return 0;
}

#else

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

/* Test-only entry point in src/xalloc.c — declared here to avoid exposing
 * the private header via the public install layout. */
extern void xalloc_set_test_failure_countdown(long n);

/* When the test is built with --coverage (gcc) or
 * -fprofile-instr-generate (clang) the .gcda / .profraw file is
 * written by an at-exit handler that _exit() bypasses, so children
 * that complete the workload without tripping the injector would
 * otherwise contribute zero coverage data.  Call the appropriate
 * flush primitive before _exit to capture them.
 *
 * The flush functions are gated on WIGNERNJ_COVERAGE, which CMake's
 * BUILD_COVERAGE option defines.  We deliberately avoid plain
 * `__attribute__((weak))` here: ELF leaves an unresolved weak
 * reference null at link time, but Mach-O (macOS) demands the
 * symbol be defined and would fail to link non-coverage builds. */
#if defined(WIGNERNJ_COVERAGE) && (defined(__GNUC__) || defined(__clang__))
extern void __gcov_dump(void)              __attribute__((weak));   /* gcc 9+ */
extern int  __llvm_profile_write_file(void) __attribute__((weak));  /* clang */
static void coverage_flush(void)
{
    if (__gcov_dump)               __gcov_dump();
    if (__llvm_profile_write_file) (void)__llvm_profile_write_file();
}
#else
static void coverage_flush(void) { /* no instrumentation */ }
#endif

typedef void (*work_fn)(void);

static void work_3j  (void) { (void)wigner3j (10,10,10,  2, 4,-6); }
static void work_6j  (void) { (void)wigner6j (10,10,10, 10,10,10); }
static void work_9j  (void) { (void)wigner9j ( 4, 4, 4,  4, 4, 4,  4, 4, 4); }
/* gaunt() argument order is interleaved: (tl1, tm1, tl2, tm2, tl3, tm3).
 * l1+l2+l3 must be even or the parity rule short-circuits before any
 * Racah-sum allocations.  Use (l,m) = (2,0)(2,0)(4,0): sum=8, valid. */
static void work_gaunt(void){ (void)gaunt    (4, 0, 4, 0, 8, 0); }

/*
 * Run `work` in a child process with the (N+1)-th allocation forced to fail
 * (or, when N is huge, with the injector effectively disabled).
 *
 * Returns:
 *   +1  child exited normally (work completed before any forced failure)
 *    0  child died with SIGABRT (the expected abort path)
 *   -1  child died with any other signal or non-zero exit (test failure)
 */
static int run_with_failure_at(work_fn work, long n)
{
    pid_t pid = fork();
    if (pid < 0) {
        perror("fork");
        return -1;
    }
    if (pid == 0) {
        /* Child: silence the diagnostic so ctest output stays readable. */
        FILE *devnull = freopen("/dev/null", "w", stderr);
        (void)devnull;
        xalloc_set_test_failure_countdown(n);
        work();
        /* If we get here the work completed without tripping the injector.
         * Flush coverage data before _exit (the at-exit handlers that
         * normally write .gcda are bypassed by _exit). */
        coverage_flush();
        _exit(0);
    }

    int status = 0;
    if (waitpid(pid, &status, 0) < 0) {
        perror("waitpid");
        return -1;
    }
    if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return +1;
    if (WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT) return 0;
    fprintf(stderr,
            "  unexpected child status (n=%ld): exited=%d code=%d signaled=%d sig=%d\n",
            n,
            WIFEXITED(status), WIFEXITED(status)  ? WEXITSTATUS(status) : -1,
            WIFSIGNALED(status), WIFSIGNALED(status) ? WTERMSIG(status)  : -1);
    return -1;
}

/*
 * Sweep the failure countdown from 0..max_n-1.  Every trial must end in
 * either SIGABRT (forced failure produced a clean abort) or a clean exit
 * (countdown was higher than the function's allocation count for these
 * arguments).  No trial may segfault.  We require at least one SIGABRT to
 * confirm the injector is actually being exercised.
 */
static void sweep(const char *name, work_fn work, long max_n)
{
    long aborts = 0, completes = 0;
    for (long n = 0; n < max_n; n++) {
        int r = run_with_failure_at(work, n);
        TEST_ASSERT(r >= 0);
        if (r == 0) aborts++;
        else if (r == +1) completes++;
    }
    /* Sanity: the injector must have triggered the abort path at least once. */
    if (aborts == 0)
        fprintf(stderr, "  %s: no aborts in sweep[0..%ld) — injector unreachable?\n",
                name, max_n);
    TEST_ASSERT(aborts > 0);

    /* And: with the injector disabled, the function must complete normally. */
    int r = run_with_failure_at(work, /*disabled*/ -1);
    TEST_ASSERT(r == +1);

    printf("  %-8s sweep[0..%ld): %ld aborted, %ld completed\n",
           name, max_n, aborts, completes);
}

int main(void)
{
    /* The exact allocation counts depend on internals, so we sweep a range
     * wide enough to cover early allocations (pfrac_init etc.) but small
     * enough to keep total wall time reasonable.  Each trial is one fork. */
    sweep("3j",    work_3j,    24);
    sweep("6j",    work_6j,    24);
    sweep("9j",    work_9j,    32);
    sweep("gaunt", work_gaunt, 24);

    SUMMARY();
}

#endif /* !_WIN32 */
