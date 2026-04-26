# Dense Polynomial Arithmetic Modulo Cyclotomic Polynomials

## Goal

Implement dense cyclotomic field multiplication as polynomial arithmetic modulo
the cyclotomic polynomial, with a clear path to faster convolution later.

## Current findings

- The dense representation stores coefficients in a vector and uses basis
  conversion/reduction.
- Dense multiplication currently performs a naive convolution-like operation and
  then reduces into the basis.
- `src/polynomial` contains standalone polynomial utilities, but they are not
  integrated into field arithmetic and are not specialized for cyclotomic
  reduction.

## Plan

1. Define the representation choice.
   - Decide whether dense elements store `phi(n)` coefficients in the basis or
     `n` coefficients before reduction.
   - Keep the public behavior unchanged during the first implementation.

2. Precompute reduction data.
   - Compute or load `Phi_n(x)`.
   - Build a map from powers `x^k` to reduced basis coordinates.
   - Cache this data per order.

3. Implement multiplication through polynomial reduction.
   - Multiply dense coefficient vectors.
   - Reduce the product modulo `Phi_n(x)`.
   - Convert into the existing basis representation.

4. Compare against existing dense code.
   - Keep the old algorithm behind a test-only or feature-gated path while the
     new one is validated.
   - Use random elements and known small fields for equality testing.

5. Optimize after correctness.
   - Improve convolution first with better loops and contiguous storage.
   - Consider FFT/NTT only if benchmarks show multiplication is the bottleneck
     and coefficient types make that practical.

## Tests and verification

- Property tests comparing old dense multiplication, sparse multiplication, and
  modulo-polynomial multiplication.
- Tests for small orders where `Phi_n(x)` is easy to inspect manually.
- Benchmarks for dense multiplication across increasing `phi(n)`.

## Risks

- Basis mismatch can cause subtle equality failures.
- Exact rational coefficient multiplication may dominate even with faster
  polynomial reduction.
- Reusing the current standalone polynomial module may add indirection without
  performance benefit.
