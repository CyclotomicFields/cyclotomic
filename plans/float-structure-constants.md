# Float-Backed Structure Constants

## Goal

Evaluate whether structure-constant multiplication can become competitive when
using float or other compact coefficient types.

## Current findings

- `src/fields/structure/mod.rs` stores structure constants as
  `Vec<Vec<Vec<Q>>>`.
- Multiplication is a triple loop over coefficient indices and rational
  coefficients.
- The file already contains a TODO about making the layout contiguous and
  vectorized.
- With `rug::Rational`, the representation pays a high cost per scalar
  operation and is unlikely to compete with dense or sparse multiplication.

## Plan

1. Add benchmarks first.
   - Measure structure, dense, and sparse multiplication for the same field
     orders and input densities.
   - Include precomputation time and steady-state multiplication time
     separately.

2. Flatten the tensor.
   - Replace `Vec<Vec<Vec<Q>>>` with a contiguous `Vec<C>`.
   - Use an index function such as `(i * phi + j) * phi + k`.
   - Keep the old implementation temporarily for correctness comparison.

3. Add sparse tensor storage.
   - Many structure constants may be zero.
   - Store nonzero `(k, coeff)` entries for each `(i, j)` pair and benchmark
     dense-tensor versus sparse-tensor traversal.

4. Add float-backed experiments.
   - Try `f64` constants and `f64` vectors as an explicitly approximate or
     checked backend.
   - Do not route this through the exact `Rational` trait unless that trait is
     redesigned.

5. Compare against dense polynomial multiplication.
   - Structure constants have expensive `O(phi^3)` multiplication unless sparse
     structure or vectorization helps.
   - Dense modulo arithmetic may beat it for many orders.

## Tests and verification

- Exact comparison of flattened/sparse constants against the current
  `rug::Rational` implementation.
- Benchmarks by order, Euler phi, and input sparsity.
- Optional debug validation for float paths against exact results.

## Risks

- Precomputation and memory may dominate for large `phi`.
- Float-backed constants may only be useful for approximate downstream code.
- The best storage layout may differ by field order.
