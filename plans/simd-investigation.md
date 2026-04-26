# SIMD Investigation

## Goal

Determine where SIMD can actually improve cyclotomic arithmetic and avoid
adding vectorized code where the data layout or coefficient type prevents a
speedup.

## Current findings

- `rug::Rational` arithmetic is not SIMD-friendly.
- Sparse maps are not SIMD-friendly.
- Nested vectors in structure constants are not ideal for vectorization.
- SIMD is most plausible for contiguous `f64`, integer, or fixed-size rational
  coefficient kernels.

## Plan

1. Establish benchmark kernels.
   - Dense vector add.
   - Dense scalar multiplication.
   - Dot products.
   - Structure-constant multiplication.
   - Matrix/vector operations if linear algebra remains a target.

2. Make data contiguous first.
   - Flatten structure constants.
   - Flatten matrices where relevant.
   - Keep coefficients in aligned contiguous vectors for SIMD candidates.

3. Choose coefficient targets.
   - `f64` for approximate or checked paths.
   - Integer coefficients for ring-of-integers experiments.
   - Fixed-size rational only if it has correct semantics and measurable use.

4. Prototype portable SIMD.
   - Prefer stable Rust APIs or a small dependency such as `wide` if needed.
   - Keep scalar fallback paths.
   - Hide architecture-specific details behind narrow internal kernels.

5. Validate speedups.
   - Compare scalar and SIMD kernels at multiple vector lengths.
   - Include realistic field operations, not only microbenchmarks.

## Tests and verification

- Exact equality tests for integer/fixed-size paths.
- Tolerance/checked tests for float paths.
- Benchmarks on at least one local machine and CI-compatible scalar fallback.

## Risks

- SIMD may not help while rational arithmetic dominates.
- Compiler auto-vectorization may already handle simple loops.
- Extra abstraction can hurt more than SIMD helps if used in the wrong layer.
