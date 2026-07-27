# Standalone GAP cyclotomic reference

This crate extracts the core algorithms from GAP's `src/cyclotom.c` behind a
small C ABI and a safe Rust ownership layer. It exists only for differential
tests and like-for-like performance work; it is not part of the published
`cyclotomic` crate.

The first milestone implements:

- construction from exponent/coefficient terms;
- conversion to GAP's Zumbroich basis;
- reduction to the minimal cyclotomic conductor;
- common-field addition;
- multiplication;
- packed, sorted sparse results;
- checked `int64_t` coefficients.

See [UPSTREAM.md](UPSTREAM.md) for the pinned GAP revision, provenance, and the
runtime substitutions made during extraction.

## Test

From the repository root:

```sh
cargo test --manifest-path reference/gap-cyclotom/Cargo.toml
```

The integration tests evaluate deterministic inputs through the extracted C
and all three Rust implementations—sparse, dense, and structure constants—then
compare exact mathematical values for addition and multiplication.

To test against the original GAP kernel, build the pinned GAP revision and
enable the optional `libgap` binding:

```sh
reference/gap-cyclotom/scripts/build-libgap.sh
export LIBGAP_ROOT="$PWD/reference/gap-cyclotom/target/libgap-d2134de71521c62512b8351c42ec16bfbac21744"
cargo test --manifest-path reference/gap-cyclotom/Cargo.toml --features libgap
```

This path does not reimplement `cyclotom.c`: it links GAP's unmodified shared
library and calls its public object API. The rational differential test covers
addition, multiplication, and division through GAP's native cyclotomic objects.

## Benchmark

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --example like_for_like > results.json
```

Add `--features libgap` with `LIBGAP_ROOT` set to include unmodified GAP in
the same integer-input benchmark.

For coefficient-for-coefficient rational timings:

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example rational_like_for_like > rational-results.json
```

For a short, deliberately harsher run covering 1% and 100% density, high
prime orders, highly composite orders, and trillion-scale numerators:

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example stress_like_for_like > stress-results.json
```

For the first representative character-table workload:

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example character_table_workloads > character-table-results.json
```

For the isolated GAP-algorithm comparison requested here:

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --example literal_gap_port > literal-c-vs-rust.json

LIBGAP_ROOT="$PWD/reference/gap-cyclotom/target/libgap-d2134de71521c62512b8351c42ec16bfbac21744" \
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example literal_gap_vs_libgap > literal-rust-vs-libgap.json
```

`literal_gap_port` compares the standalone C extraction with its direct Rust
translation, using the same checked `i64` coefficient model.
`literal_gap_vs_libgap` compares only that Rust translation with unmodified
GAP. The Rust source mirrors `ConvertToBase`, `Cyclotomic`, `FindCommonField`,
and `ProdCyc`: it keeps the dense reusable `ResultCyc` buffer, packed sorted
outputs, the smaller-right-operand loop order, the zero/±1 coefficient cases,
cached field properties, and eager minimal-conductor reduction. Exact tests
compare canonical order and every packed term with the C extraction, then
compare mathematical equality with unmodified libgap.

The translation intentionally uses checked `i64`, so this is a comparison of
the kernel algorithm for GAP immediate-integer-sized inputs. It does not yet
translate GAP's tagged `Obj` representation, garbage collector, or automatic
promotion to arbitrary integers and rationals. Consequently the two paths run
the same mathematical control flow but need not have identical allocation and
coefficient-operation costs.

For the adaptive coefficient implementation:

```sh
LIBGAP_ROOT="$PWD/reference/gap-cyclotom/target/libgap-d2134de71521c62512b8351c42ec16bfbac21744" \
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example tagged_gap_vs_libgap > tagged-rust-vs-libgap.json

LIBGAP_ROOT="$PWD/reference/gap-cyclotom/target/libgap-d2134de71521c62512b8351c42ec16bfbac21744" \
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --example tagged_character_tables > tagged-character-tables.json
```

This path dispatches an all-integer cyclotomic to the checked packed `i64`
kernel. An overflowing operation or rational input promotes transparently to
the exact kernel. Its scalar `TaggedRational` is one machine word: the low bit
marks an immediate signed integer, while the second tag bit distinguishes
boxed `rug::Integer` from boxed `rug::Rational`. Integral overflow therefore
does not enter rational arithmetic. Rational imports are staged through a
common denominator when that denominator and the scaled numerators fit `i64`.
The hot accumulator mutates promoted GMP values in place.

The Rust port also includes the low-level optimizations used by the current
benchmarks:

- modulo-free convolution indices;
- inline storage for up to four packed terms;
- adaptive generation-stamped/touched-index scratch clearing;
- cached high-order factorization, embedding, and basis-conversion plans;
- delayed normalization and direct integer extraction for character scalar
  products;
- a batched multi-output character kernel;
- exact, bound-checked CRT reconstruction for overflowing integer dot
  products, with a big-integer fallback when the fixed prime budget is
  insufficient.

`tagged_character_tables` imports GAP's exact character tables for `A5`,
`SL(2,5)`, and `PSL(2,11)` into that adaptive kernel. It precomputes the
class-size-weighted conjugates, then decomposes every unordered pair of
irreducible characters. Before timing, every multiplicity vector is checked
against GAP. Each short timing iteration is one complete sweep, and the JSON
reports both time per sweep and time per tensor decomposition. Table import and
the invariant weighted conjugates are outside both timings; result construction
and all temporary cyclotomic arithmetic remain inside.

The same workload is also a proper Cargo benchmark target:

```sh
LIBGAP_ROOT="$PWD/reference/gap-cyclotom/target/libgap-d2134de71521c62512b8351c42ec16bfbac21744" \
cargo bench \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --features libgap \
  --bench tagged_character_tables
```

It reports median and range over 15 calibrated samples per table and
implementation. `CYCLOTOMIC_BENCH_TABLE`, `CYCLOTOMIC_BENCH_IMPL`,
`CYCLOTOMIC_BENCH_SAMPLE_MS`, and `CYCLOTOMIC_BENCH_SAMPLES` filter or extend
the run for profiling.

`PreparedRepresentationRing::build` computes the full symmetric table of
irreducible tensor-product structure constants once. It then supports
allocation-free basis-product lookup and checked multiplication of arbitrary
signed virtual-character coefficient vectors. The benchmark reports these as
`rust-ring-lookup` and `rust-ring-vector`; precomputation is intentionally
outside those hot timings.

One indicative 15-sample run with a 20 ms target per sample gave:

| table | GAP | Rust arithmetic | speedup | ring basis lookup | ring vector |
| --- | ---: | ---: | ---: | ---: | ---: |
| A5 | 3.547 us | 0.844 us | 4.20x | 0.001 us | 0.062 us |
| SL(2,5) | 7.151 us | 1.119 us | 6.39x | 0.001 us | 0.133 us |
| PSL(2,11) | 5.378 us | 1.133 us | 4.75x | 0.001 us | 0.087 us |

These numbers are deliberately short and indicative. Correctness against
unmodified GAP is checked before every timed table.

This asks unmodified GAP to construct the character tables of `A5`, `SL(2,5)`,
and `PSL(2,11)`. Before timing, it transfers the exact table entries and class
sizes into both Rust representations. The timed operation decomposes the
tensor product of two non-rational irreducible characters: it pointwise
multiplies their character values, takes the weighted scalar product with
every irreducible character, and returns the integer multiplicities. Every
Rust result is checked exactly against GAP before measurements are emitted.
The Rust records expose both `prepared_arithmetic_only`, where table conversion
and invariant weighted conjugates are prepared before timing, and
`allocation_inclusive_end_to_end`, where conversion, precomputation, scratch
construction, arithmetic, and canonical integer output all occur inside each
iteration. “Arithmetic only” excludes setup allocation but still includes
temporary and result allocation intrinsic to the current arithmetic API. The
GAP record is labeled `native_representation`; its public libgap call includes
GAP's own scalar-product work and result-list construction.

The stress benchmark uses a 10 ms calibration target. It omits the structure
implementation at large orders because that implementation eagerly constructs
and retains cubic-size structure-constant tables; order 10009 also omits the
quadratic dense implementation.

The benchmark takes roughly ten seconds after compilation, constructs identical
deterministic rational inputs before timing, guarantees at least one nonzero term
per operand, and emits JSON records for:

- `gap_extracted_i64` uses checked machine integers and returns a canonical,
  minimal-conductor result on every operation.
- `gap_unmodified_libgap` uses original GAP objects, automatic small/big
  integer and rational arithmetic, and the original cyclotomic kernel.
- `rust_sparse_packed_hybrid` keeps machine-size integers inline, promotes
  exactly to GMP rationals when necessary, and reuses a packed multiplication
  scratchpad.
- `rust_dense_rational` uses the dense root-of-unity representation.
- `rust_structure_rational` uses the nested structure-constant implementation.

The libgap-enabled benchmark has coefficient-for-coefficient semantic parity
with the Rust implementations. The standalone extraction remains an
integer-only algorithm experiment.

## Representative workload plan

The synthetic order, density, and coefficient-size sweeps are retained as an
adversarial envelope, not treated as a model of typical GAP use. GAP and its
packages usually consume cyclotomics as entries in complete character rows,
class functions, and matrices. In those workloads multiplication appears
inside weighted scalar products, tensor products, traces, Galois actions, and
repeated powers. Ordinary irreducible character values are algebraic integers;
rational coefficients are more commonly introduced later by averaging,
projection, and division by a group order.

The representative suite will stay short and will be developed in this order:

1. **Character-table tensor decomposition — implemented.** Replay actual GAP
   tables, pointwise-multiply two irreducible characters, and compute all
   multiplicities with weighted scalar products.
2. **Character scalar products.** Time a single weighted
   `sum(class_size[i] * chi[i] * conjugate(psi[i])) / group_order`, plus batches
   of all pairwise row products.
3. **Repeated tensor powers.** Repeatedly multiply a character row and
   decompose it, matching representation-growth and tensor-power calculations.
4. **HeLP-style field traces.** Compute
   `Trace(chi[i] * E(k)^(-l))` for prime and highly composite `k`.
5. **Galois orbits and conjugation.** Apply `GaloisCyc` across whole character
   rows, including complex conjugation and orbit closure.
6. **Root-sparse high-order cases.** Use a single root, a conjugate pair, and
   small orbit sums at orders such as 1009, 2520, and 10009. These complement
   the less realistic “1% density” case, which already contains 101 terms at
   order 10009.
7. **Additional table corpus.** Add cyclic and dihedral families, larger
   generated tables, and optional CTblLib tables when that package is present.

Each workload should expose two separately named modes where applicable:

- **canonical semantic operation**, where every implementation produces the
  same minimal-conductor result;
- **native representation operation**, where each implementation keeps its
  normal eager or lazy normalization strategy.

Construction, conversion, allocation-inclusive execution, and arithmetic-only
execution must be reported separately rather than silently mixed. Correctness
checks run before timing and compare exact cyclotomic values or exact integer
multiplicities.

## Remaining application-level work

The extracted kernel, adaptive integer/rational runtime, exact overflow paths,
character batching, and precomputed representation-ring mode are implemented.
The next work is deliberately application-facing:

1. Continue the representative workload plan with repeated tensor powers,
   HeLP-style traces, and Galois-orbit batches.
2. Extend the table corpus to generated cyclic/dihedral families and optional
   CTblLib tables.
3. Report cold precomputation separately from hot representation-ring queries
   for larger tables.
4. Add subtraction, inversion, and general Galois action to the small public
   reference API when a benchmark workload needs them.
