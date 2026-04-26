# Field of Fractions of the Ring of Integers

## Goal

Investigate an exact representation where elements are stored as fractions of
cyclotomic integers, making multiplication mostly integer arithmetic and making
some inversion workflows simpler.

## Current findings

- Current field elements are rational coefficient vectors or maps.
- This is general, but every coefficient operation pays rational arithmetic
  costs.
- For multiply-heavy workloads, storing an integer numerator plus a shared
  denominator may reduce overhead.

## Plan

1. Define the representation.
   - Numerator: integer coefficients in a cyclotomic integer basis.
   - Denominator: positive integer.
   - Optional normalization state.

2. Implement basic arithmetic.
   - Addition: cross-multiply denominators.
   - Multiplication: multiply numerators, multiply denominators.
   - Negation: negate numerator.
   - Equality: compare cross-products or normalized forms.

3. Decide normalization strategy.
   - Eager normalization reduces growth but costs gcd work.
   - Lazy normalization may be faster for multiply-heavy workloads.
   - Add explicit `normalize()` for API and benchmark control.

4. Implement inversion.
   - Use the product of Galois conjugates or another exact norm-based method.
   - Keep the output in numerator/denominator form.

5. Compare with current rational representation.
   - Multiplication-only workloads.
   - Mixed add/mul workloads.
   - Inversion-heavy workloads.
   - Large coefficient growth cases.

## Tests and verification

- Equivalence tests against existing exact dense/sparse representations.
- Normalization invariants.
- Randomized property tests for field axioms.
- Benchmarks that track coefficient growth.

## Risks

- Denominator and numerator growth may erase gains.
- Inversion can still be expensive.
- This representation may be best as a specialized backend, not the default.
