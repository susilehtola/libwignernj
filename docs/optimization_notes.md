# Optimization notes

This document records the optimization decisions taken in libwignernj —
what's been tried, what worked, and what didn't.  Future contributors
considering "obvious" speedups should check here first; several of the
items below look like good ideas on paper and have been implemented,
benchmarked, and reverted as net-negative.

The optimisation methodology is bench-driven: every perf-relevant
change is benched against the parent commit using paired alternating
runs (the variant under test and the parent are run in alternation
to control for thermal/frequency drift on consumer laptops), and the
headline numbers go in the commit message.  Benches live under
`benchmarks/` and are part of the source tree (see CLAUDE.md for
the policy distinguishing library-versioning benches from
paper-scope comparisons).

## Optimisations that landed

These commits stack on `granlund-mueller-3j` (or wherever the work
ultimately merges); inspect `git log --grep` for the precise ranges.

| Commit (or branch tag) | What                                              | Headline gain                              |
|------------------------|---------------------------------------------------|--------------------------------------------|
| (previously merged)    | Two-pass Racah sum (Pass 1 builds pfracs / LCM, Pass 2 accumulates LCM-scaled scaled bigint) | foundational                                |
| (previously merged)    | Per-thread cached `wignernj_scratch_t`              | allocation-free steady state                |
| (previously merged)    | Per-thread factorial-decomposition cache          | replaces per-call `legendre_valuation` loop |
| (previously merged)    | uint64-batched prime-power accumulator (`pfrac_mul_pow_into_acc`) | ~62% time saved in `pfrac_lcm_scaled_product` at j=4000 |
| (previously merged)    | Hardware divq via inline asm on x86-64            | ~3× faster `bigint_div128` than libgcc generic |
| (previously merged)    | Algorithm D for arbitrary 64-bit divisors         | corrects + speeds up the no-divq backend    |
| (previously merged)    | Small-integer ratio recurrence in 3j/6j/9j Pass 2 | replaces O(π(j)) prime-power expansion per term with one batched mul/div |
| (previously merged)    | Hensel exact division (`bigint_div_u64_exact`)    | 6j +1.13–1.20×; 9j ~flat; 3j untouched (Hensel setup didn't pay off at small bigints) |
| (previously merged)    | `BUILD_LTO=ON` default                            | small but free perf win                     |
| (previously merged)    | Drop explicit `-O2` override; default Release    | unblocks `-O3` for users who don't pass `-DCMAKE_BUILD_TYPE` |
| (previously merged)    | `restrict` qualifier on pfrac vector adds         | unblocks SSE2 vpaddd auto-vectorisation     |
| (recent)               | Möller-Granlund division for non-divq backends    | no-divq AlgoD: +16–21% on 3j at j ≥ 200     |
| (recent)               | Gaunt Pass-2 small-integer ratio recurrence       | 3.5–4.2× headline at l ≥ 50 across all backends |
| (recent)               | Gaunt Pass-1 cache skip (peeled rolling buffer)   | small bench win + multi-MiB real-world cache savings |

## Optimisations that didn't land (and why)

### Toom-3 multiplication for the in-tree bigint backend

Implemented as Bodrato 5-point evaluation (0, 1, −1, 2, ∞) on top of
the existing Karatsuba routine for the sub-mults, with signed-magnitude
handling of `vm`, exact division by 3 in the interpolation, and a
`TOOM3_THRESHOLD` cutoff below which it falls through to Karatsuba.
~432 LOC added to `src/bigint.c`.  All existing tests passed.

**Microbench (m × m bigint multiplication, isolated):** Toom-3 was
*slower* than Karatsuba almost everywhere below m = 512.  At m = 192
(close to the j ~ 3000 hot range), Toom-3 was 41% slower; at m = 256
(j ~ 4000), 18% slower.  Only at m = 512 did Toom-3 pull ahead, by 14%.

**End-to-end 3j sweep on x86-64 divq, paired against the parent
commit:** geo-mean speedup 0.91× at j ≥ 100, 0.96× at j ≥ 2000, 1.01×
(within noise) at j = 3995.  At j = 30 (where Toom-3 isn't even
invoked because the bigint is too small) the sweep was 29% slower —
a tell that the additional 432 lines of code change icache layout for
the rest of the binary, and the dispatch cost ripples into the hot
path of every call.

**Underlying mechanism:** Toom-3's asymptotic crossover with
Karatsuba on this hardware is at m ≈ 512, while libwignernj's typical
operating range is m ≤ ~250 (j = 4000).  The asymptotic gain doesn't
materialise in the relevant regime, and the dispatch + code-bloat cost
dominates.  FLINT covers the very-large-j case where sub-quadratic
multiplication actually wins.

**Reverted.**  Don't re-propose without bench data showing a different
hardware crossover or a substantially smaller code footprint.

### Compressed `pfrac_t.exp` / factorial-cache rows (split int8/int16/int32 tier)

Phase-1 implementation: keep `pfrac_t.exp` as `int*` for now, compress
only the cached factorial-decomposition rows into a two-tier layout
(int32 for primes p < ~158 where v_p(MAX_FACTORIAL_ARG!) can exceed
127, int8 for higher primes).  Generator emits `g_nprimes_lo` to
prime_table.inc; `pfrac_mul_factorial`/`pfrac_div_factorial` do two
loops, one per tier.  ~150 LOC across `pfrac.c`, the prime-table
generator, and `primes.h`.

**Cache profile at j = 4000:** worked as designed — cache misses
dropped 38% (LLC) / 81% (L1) on a 30-call profile.

**Wall-time at j = 4000:** ~5% slower despite the cache-miss reduction.

**Wall-time microbench:**

| j   | parent | phase 1 |
|-----|--------|---------|
| 1   | 80.8   | 105.0   | (−23%)
| 5   | 217    | 254     | (−15%)
| 15  | 462    | 568     | (−19%)
| 30  | 1076   | 1445    | (−26%)
| 100 | 6122   | 7510    | (−18%)

**Underlying mechanism:** the per-call hot path of
`pfrac_mul_factorial` got 3–4 extra struct dereferences (`r->row_lo`,
`r->row_hi`, `r->width_lo`, `r->width_hi`) compared to the previous
single `int*` row plus `g_pi_table[n]` lookup.  At small j a wigner3j
call is only ~80 ns total — a few extra L1 loads per factorial-cache
access dominate.  The cache-pressure savings only pay off when
`n_terms × max_idx` exceeds L1, i.e. j ≳ 1000.  For typical
workloads (j ≤ 200) this is a clear net negative.

A "phase 2" compression of the `pfrac_t.exp` storage itself (rather
than just the cache rows) was not implemented but is expected to
share the same shape, amplified by the larger surface — every
`exp[i]` access in the codebase becomes a two-tier dispatch.

**Reverted.**  Don't re-propose without a strategy that doesn't add
hot-path dereferences.

### Small-bigint Pass-2 specialisation

Stack-allocated `uint64_t scaled[N_SMALL]` + static-inline
`bigint_mul_u64_small` / `bigint_div_u64_small` / `bigint_add_small`
in a new `src/bigint_small.h`, with `wigner3j_exact` Pass 2
dispatching to the small path based on the seed bigint size and
falling back to the existing `bigint_t` path on growth.  Both 3j
small_pos / small_neg accumulators ran on stack arrays.

**Microbench:** at j = 100 the small path was 31% slower than the
parent; the untouched 6j and 9j routines also slowed measurably
because the `wigner3j_exact` body bloat affected whole-binary code
layout under LTO.

**Underlying mechanism:** LTO is already aggressively inlining the
existing `bigint_mul_u64` / `bigint_div_u64` / `bigint_add` into
`wigner3j_exact`'s Pass 2 — the function-call overhead the small
path was meant to eliminate had already been eliminated by the
compiler.  Adding a second specialised inline body simply doubled
the function size and pushed everything out of icache further.

**Reverted.**  Don't re-propose static-inline copies of the bigint
primitives.  If a bigint primitive really is too slow, attack it in
its own definition (e.g. with ADCX/ADOX intrinsics for the inner
mul-with-carry chain), not by copying it next to a caller.

## Out of scope

These are excluded by project policy, not by bench data:

- **Intra-call threading** (OpenMP, pthreads, threaded BLAS).
  Embeddability mandate trumps the perf gain — see `CLAUDE.md`
  Project policy.
- **Reading GPL/LGPL competitor source** (WIGXJPF, FASTWIGXJ, etc.)
  for inspiration.  libwignernj is BSD-3-Clause and stays clean-room.
- **GPU offload.**  Adds dependency surface; the library is designed
  to be a small, self-contained, embeddable C library.

## Things to consider only if a real perf gap is identified

Listed without endorsement — none of these have been benched but each
has at least a plausible argument for being net-positive on some
workload.  Bench against the parent commit before assuming any of
them helps.

- **PGO build option.**  5–10% across the board on the workload it's
  trained against.  Compiler-dependent (GCC, Clang, MSVC each have
  their own workflow), CPU-microarchitecture-dependent (a profile
  collected on Intel may not help on AMD), and training-corpus-
  dependent (a profile from `bench_term_cache` may regress workloads
  with a different j-distribution).  Probably better as a *user-side*
  recipe documented in the README than as a library build option.
- **ADCX/ADOX intrinsics in the `bigint_mul_u64` inner loop.**
  Two parallel carry chains via independent flag registers; ~2×
  throughput on hardware with ADX (Broadwell+, Zen+).  At ~24% of
  cycles on j=4000 in profile, this is the largest remaining
  per-call hot spot.  Cost: another arch-conditional path.
- **AVX-512 VPMADD52** (52-bit limb arithmetic).  Potentially 2–3×
  on bigint multiplications on Ice Lake+ Intel and Zen 4+ AMD.
  Massive rewrite (~thousands of lines).  Big-mul optimisations
  haven't translated to wall time in this library so far (Toom-3
  was a clean illustration); be cautious.
- **m = (0, 0, 0) fast path.**  When all m's are zero the symbol is
  highly symmetric; many factorials cancel.  A direct closed-form
  evaluation could be much faster than the general Pass 2.  Niche
  but real for cosmology pseudo-Cℓ workloads.
- **Cache-skip Pass 1 for 3j/6j/9j** (the analogue of the Gaunt
  Pass-1 cache-skip change for the 3-/6-/9-j evaluators).  Small bench win + real cache-pollution savings for
  callers that interleave libwignernj with other code.  Earlier
  attempts showed mixed bench data; needs re-attempting with the
  paired-alternating methodology.
