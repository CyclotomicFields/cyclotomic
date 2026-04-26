# Use f64 for Integer-Valued Paths

## Goal

Use floating-point arithmetic only in paths where the mathematical result is
known to be integral, and only when the result can be checked or bounded safely.

## Current findings

- `FloatRational` exists in `src/fields/rational.rs`, but it is not a sound
  exact rational implementation.
- It uses floating-point storage while implementing APIs that imply exact
  rational semantics, including numerator and denominator access.
- The project needs exact cyclotomic arithmetic, so unchecked `f64` substitution
  would be a correctness regression.

## Plan

1. Identify candidate integer-valued computations.
   - Structure constant construction.
   - Basis conversion/reduction coefficients.
   - Trace-like or inner-product-like computations if present later.
   - Any algorithms that mathematically end in `Z`, not just `Q`.

2. Define a separate approximate coefficient path.
   - Do not make `f64` implement the exact `Rational` trait unless the trait is
     narrowed.
   - Introduce a trait or strategy type for approximate arithmetic with explicit
     rounding/checking.

3. Add checked conversion.
   - Convert an `f64` result to an integer only if it is within a strict
     tolerance of an integer and the expected magnitude is below the exact
     integer range for `f64`.
   - Fall back to exact `rug::Rational` when the check fails.

4. Prove the safe cases locally.
   - Document why each float-backed path has an integer result.
   - Add debug assertions or optional verification against exact arithmetic.

5. Benchmark the useful cases.
   - Compare exact-only, float-checked, and float-with-exact-fallback paths.
   - Test small and large orders separately; float overhead may not pay off for
     small fields.

## Tests and verification

- Property tests comparing float-checked paths with exact rational paths.
- Stress tests near `f64` precision limits.
- Benchmarks that include failure/fallback rates.

## Risks

- Rounding can silently corrupt exact arithmetic if checks are too weak.
- The fast path may be slower if validation dominates arithmetic.
- This work depends on clearer coefficient abstractions.
