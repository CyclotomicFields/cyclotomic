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

The integration tests evaluate deterministic inputs through both the extracted
C and the existing Rust sparse implementation, then compare their mathematical
values.

## Benchmark

```sh
cargo run --release \
  --manifest-path reference/gap-cyclotom/Cargo.toml \
  --example like_for_like > results.json
```

The benchmark constructs identical deterministic integer inputs before timing
and emits JSON records. The current labels are deliberately explicit:

- `gap_extracted_i64` uses checked machine integers and returns a canonical,
  minimal-conductor result on every operation.
- `rust_sparse_rational` uses GMP-backed `rug::Rational`, clones operands to
  satisfy its mutating API, and leaves multiplication in its normal
  noncanonical internal representation.

Those are real API costs, but they are not yet coefficient-for-coefficient
parity. GMP rational support in the C extraction and a libgap oracle are needed
before drawing a final GAP-versus-Rust conclusion.

## Remaining work

1. Add GMP integer and rational coefficients with GAP-style small-integer fast
   paths and promotion.
2. Add subtraction, scalar multiplication, Galois action, and inversion.
3. Bind unmodified libgap as a correctness and performance oracle.
4. Extend differential tests to rational coefficients and GAP's own edge cases.
5. Split arithmetic-only, allocation-inclusive, and batch-workload benchmarks.
