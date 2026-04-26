# Reduce Copy/Paste with Useful Abstractions

## Goal

Remove duplicated dense/sparse logic by extracting abstractions that match the
mathematics and the storage differences, without hiding performance-critical
details behind heavy generic layers.

## Current findings

- Dense and sparse basis conversion code duplicate similar ideas.
- Dense code is specialized around `i64` exponents and `rug::Rational`
  coefficients.
- Sparse code is more generic over exponent and coefficient types.
- The current algebra trait shape also creates duplicated clone-heavy call
  patterns.

## Plan

1. List duplicated behavior.
   - Basis membership.
   - Basis reduction.
   - Order promotion.
   - Conjugation.
   - Scalar multiplication.
   - Generic field axiom tests.

2. Extract basis/reduction logic first.
   - Introduce a small `BasisReducer` or equivalent internal abstraction.
   - Keep dense and sparse storage separate.
   - Make the abstraction return representation-friendly operations rather than
     forcing everything through one container type.

3. Extract coefficient-store behavior second.
   - Dense: contiguous vector.
   - Sparse: map from exponent to coefficient.
   - Keep insertion/removal/zero-elision policies explicit.

4. Refactor tests.
   - Replace copy/pasted field tests with shared property-test helpers that do
     not force unnecessary clones.

5. Avoid over-generalizing.
   - Do not add an abstraction unless it removes repeated code or enables a
     planned backend.
   - Keep hot arithmetic loops inspectable and benchmarked.

## Tests and verification

- Existing `cargo test`.
- Dense/sparse equivalence tests for random elements.
- Benchmarks before and after each abstraction step.

## Risks

- Generic abstractions can make compiler errors and performance worse.
- Dense and sparse code may look similar but need different zero-handling and
  iteration behavior.
- This should be incremental; a large rewrite would be hard to review.
