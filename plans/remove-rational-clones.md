# Remove Clones in Rational Implementations

## Goal

Remove avoidable cloning from rational arithmetic while preserving exactness
and the existing field behavior.

## Current findings

- The core algebra traits in `src/fields/mod.rs` use mutating methods such as
  `add(&mut self, z: &mut Self)` and `mul(&mut self, z: &mut Self)`.
- That API gives implementations a mutable RHS even when the RHS only needs to
  be read. For `rug::Rational`, this encourages `z.clone()` before applying
  arithmetic.
- `src/fields/rational.rs` has a local TODO documenting this exact issue.
- Sparse and dense element methods also clone coefficients or whole elements in
  scalar multiplication, basis conversion, equality checks, and tests.
- `FixedSizeRational::is_zero` currently calls itself recursively; that should
  be fixed before relying on it for performance work.

## Plan

1. Add benchmark coverage before refactoring.
   - Include rational add, multiply, negation, inversion, scalar multiplication,
     sparse element multiplication, and dense element multiplication.
   - Capture clone-heavy baseline behavior with current code.

2. Introduce borrowed assign-style operations.
   - Add methods with shapes like `add_assign_ref(&mut self, rhs: &Self)`,
     `mul_assign_ref(&mut self, rhs: &Self)`, `neg_assign(&mut self)`, and
     `inv_assign(&mut self)`.
   - Keep compatibility shims initially so the full codebase does not need to
     move in one commit.

3. Refactor `rug::Rational` implementation.
   - Use borrowed `AddAssign` and `MulAssign` implementations where `rug`
     supports them.
   - Avoid cloning the RHS in common add/multiply paths.
   - Keep exact numerator/denominator access exact, even if those methods still
     need to allocate owned integers.

4. Fix local rational implementations.
   - Fix `FixedSizeRational::is_zero`.
   - Audit `FloatRational`; it should not pretend to be an exact rational if it
     cannot supply exact numerator/denominator semantics.

5. Refactor callers.
   - Update dense, sparse, structure-constant, and test macro code to pass
     immutable RHS references where possible.
   - Remove whole-element clones in scalar multiplication where coefficient-wise
     mutation is enough.

6. Remove old compatibility methods if this is allowed as a breaking API
   change. Otherwise, keep deprecated wrappers.

## Tests and verification

- `cargo test`
- New property tests comparing old and new arithmetic paths.
- Benchmarks before and after the refactor.
- Optional: run `cargo clippy` to catch suspicious clone patterns.

## Risks

- This is likely a public API change unless compatibility wrappers remain.
- `rug` borrowed operation support may vary by operation and version.
- Removing clones from tests may require more careful value ownership.
