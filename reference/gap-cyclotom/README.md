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

This asks unmodified GAP to construct the character tables of `A5`, `SL(2,5)`,
and `PSL(2,11)`. Before timing, it transfers the exact table entries and class
sizes into both Rust representations. The timed operation decomposes the
tensor product of two non-rational irreducible characters: it pointwise
multiplies their character values, takes the weighted scalar product with
every irreducible character, and returns the integer multiplicities. Every
Rust result is checked exactly against GAP before measurements are emitted.
Table construction and conversion are deliberately outside the timed region.
This first benchmark is labeled `native_representation`: GAP returns canonical
integer multiplicities, while the Rust implementations retain their normal
lazy cyclotomic representation. A separately labeled canonical-result mode
remains planned below.

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
- `rust_sparse_rational` uses GMP-backed `rug::Rational`, clones operands to
  satisfy its mutating API, and leaves multiplication in its normal
  noncanonical internal representation.
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

## Remaining implementation work

1. Add GMP integer and rational coefficients to the standalone extraction.
2. Add subtraction, scalar multiplication, Galois action, and inversion there.
3. Extend differential tests to more GAP edge cases and large coefficients.
4. Continue the representative workload plan above.
5. Split arithmetic-only, allocation-inclusive, and batch-workload benchmarks.
