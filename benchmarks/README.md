# Benchmarks

Library-versioning benchmarks for `libwignernj`.  These harnesses
measure libwignernj against prior versions of itself, supporting the
paired-alternating-runs methodology described in
[../docs/optimization_notes.md](../docs/optimization_notes.md).
None of them require external libraries; building requires only
`libwignernj` itself.

## What's here

| File | Measures |
|---|---|
| `bench_term_cache.c` | per-symbol microbench (3j, 6j, 9j; small-to-medium *j*; min-of-trials) |
| `bench_sweep.c`      | logarithmic *j* sweep for one symbol family at a time (3j, 6j, 9j, gaunt) |
| `bench_mul.c`        | `bigint_mul_ws` primitive at sweeping *m* (limb count) |
| `bench_div128.c`     | `bigint_div128` primitive across all three dispatch arms |
| `bench_div128_gm.c`  | Möller-Granlund vs hardware `divq`, swept over divisor classes and *m* |
| `profile_3j_4000.c`  | single-driver loop for `perf record` at the j=4000 hot point |
| `profile_6j_9j.c`    | same shape for 6j/9j |

## Building and running

Each bench compiles standalone against the libwignernj static
library.  From the repository root:

```sh
cmake -B build -DBUILD_SHARED_LIBS=OFF
cmake --build build --parallel

# example: 3j sweep up to j=4000, 0.5s budget per j, 60s ceiling per call
gcc -O3 -DNDEBUG -Iinclude benchmarks/bench_sweep.c \
    -o /tmp/bs build/libwignernj.a -lm
/tmp/bs 3j 4000 0.5 60
```

For paired-alternating bench against a prior commit, build twice (once
on each git ref), keep both binaries, and alternate runs to control
for thermal/frequency drift on consumer hardware:

```sh
for run in 1 2 3; do
  /tmp/bs_HEAD       3j 4000 0.5 60 > /tmp/sw_head_r$run.dat
  /tmp/bs_PARENT     3j 4000 0.5 60 > /tmp/sw_parent_r$run.dat
done
```

Take the minimum-across-runs at each j when comparing.

## Comparative benchmarks against external libraries

Head-to-head comparisons of `libwignernj` against
[WIGXJPF](https://fy.chalmers.se/subatom/wigxjpf/) and the GNU
Scientific Library (GSL) -- the data behind the published timing and
accuracy tables -- are not included in this directory.  That harness
links against WIGXJPF and GSL, which are outside `libwignernj`'s own
dependency surface; it is published as supplementary material with
the corresponding paper.

## License

BSD 3-Clause, like the rest of `libwignernj`.
