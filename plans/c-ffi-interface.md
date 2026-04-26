# C FFI Interface

## Goal

Expose a stable C ABI for calling cyclotomic field operations from C, C++, and
other languages that can bind to C.

## Current findings

- The Rust API uses traits, generic types, `rug` values, maps, and Rust-owned
  allocation.
- None of those should cross an `extern "C"` boundary directly.
- The README claims C/C++ use is a goal, but there is no ABI-safe layer yet.

## Plan

1. Decide crate layout.
   - Add a `cdylib` build target or a small companion crate.
   - Keep the pure Rust library API separate from ABI concerns.

2. Define opaque handles.
   - `cyclotomic_field_t`
   - `cyclotomic_elem_t`
   - Possibly `cyclotomic_error_t` or thread-local error strings.

3. Add lifecycle functions.
   - Create/free field.
   - Create/free element.
   - Clone element.
   - Clear/free returned strings or buffers.

4. Add constructors and serialization.
   - Construct zero/one/root of unity.
   - Construct from coefficient arrays.
   - Parse from a stable textual format.
   - Serialize to a stable textual or array format.

5. Add operations.
   - Add, subtract, multiply, negate, invert.
   - Equality.
   - Conjugation.
   - Optional scalar multiplication.

6. Add error handling.
   - No panics across the ABI.
   - Return status codes.
   - Provide a way to retrieve diagnostic messages.

7. Generate and test headers.
   - Maintain a checked-in C header or generate one with `cbindgen`.
   - Add C smoke tests that compile and link in CI.

## Tests and verification

- Rust tests for the ABI wrapper.
- C integration test compiled by CI.
- Leak checks where practical.
- Header compatibility check.

## Risks

- Public ABI stability is harder than Rust API stability.
- Ownership rules must be extremely clear.
- `rug` and native GMP-related linkage may complicate packaging.
