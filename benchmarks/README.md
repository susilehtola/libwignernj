# Comparative benchmarks

This directory contains the benchmark code used in
[doc/paper.tex](../docs/paper.tex), Section "Benchmarks", to compare
`libwignernj` against the GNU Scientific Library (GSL) and the
prime-factorization reference implementation
[WIGXJPF](https://fy.chalmers.se/subatom/wigxjpf/).

## What is measured

`bench_compare.c` evaluates each of three symbols on a sweep of
**all-equal-*j*** inputs over a wide range of *j*:

| Symbol | Input |
|---|---|
| 3*j* | $(j,j,j;\,m,-m,0)$  with $m$ swept over $-j,\ldots,+j$ |
| 6*j* | $\{j,j,j;\,j,j,j\}$ (single fixed input per *j*) |
| 9*j* | $\{j,j,j;\,j,j,j;\,j,j,j\}$ (single fixed input per *j*) |

Reported per evaluation: wall-clock time in nanoseconds, plus the
floating-point sum of the returned values across the inner loop as
a sanity check that the libraries agree on the numerical result.

## Reproducibility

For a fair comparison, all three libraries should be built with the
**same compile flags**.  The paper's published numbers used the
default Fedora optimization flags, which on a Fedora-derived system
can be queried directly with

```sh
rpm -E %optflags
```

The Makefile in this directory takes its default `CFLAGS` from
`rpm -E %optflags` when `rpm` is available, and falls back to a
verbatim copy of the same flag string otherwise.  Override on the
command line if needed (e.g.\ `make CFLAGS="-O2 -march=native"`).

The same `CFLAGS` should be used to rebuild WIGXJPF
(`make -C $WIGXJPF_DIR CFLAGS="$(rpm -E %optflags) -fPIC -I inc/ -I cfg/ -I gen/"`)
and `libwignernj` (`cmake -DCMAKE_C_FLAGS="$(rpm -E %optflags)"`).
GSL is typically installed as a system package; on Fedora it is
already built with these flags.

## Build & run

```sh
# 1. Build libwignernj normally (it picks up CFLAGS through CMake)
cmake -B ../build -DBUILD_TESTS=OFF
cmake --build ../build

# 2. Build WIGXJPF (https://fy.chalmers.se/subatom/wigxjpf/)
tar xfz wigxjpf-1.13.tar.gz
make -C wigxjpf-1.13

# 3. Build and run the benchmark
make WIGXJPF_DIR=/path/to/wigxjpf-1.13 \
     WIGNERNJ_DIR=/path/to/lib3j/build \
     run
```

For best repeatability, run on a quiet machine with the CPU
frequency governor pinned to `performance`.  The benchmark loop
does not flush caches or randomize the input order; it is intended
for an order-of-magnitude comparison, not a microbenchmark down to
the cycle level.

## License

BSD 3-Clause, like the rest of `libwignernj`.
