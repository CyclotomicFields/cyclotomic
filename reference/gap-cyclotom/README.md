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

The benchmark takes roughly ten seconds after compilation, constructs identical
deterministic integer inputs before timing, guarantees at least one nonzero term
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

## Remaining work

1. Add GMP integer and rational coefficients to the standalone extraction.
2. Add subtraction, scalar multiplication, Galois action, and inversion there.
3. Extend differential tests to more GAP edge cases and large coefficients.
4. Split arithmetic-only, allocation-inclusive, and batch-workload benchmarks.
