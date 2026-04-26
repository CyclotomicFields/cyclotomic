# C and C++ Usage Instructions

## Goal

Document how to call the library from C and C++ once the C ABI exists.

## Current findings

- README has a placeholder for C/C++ usage.
- There is no current FFI API to document accurately.
- Usage instructions should be written against tested examples, not aspirational
  function names.

## Plan

1. Wait for the C ABI surface.
   - This document depends on `c-ffi-interface.md`.
   - Do not publish detailed C/C++ instructions until function names, ownership,
     and error handling are stable.

2. Add minimal C example.
   - Create a field.
   - Create two elements.
   - Multiply or add them.
   - Print/serialize result.
   - Free every handle.

3. Add minimal C++ example.
   - Show direct C API usage first.
   - Optionally add a small RAII wrapper example if the API is stable enough.

4. Document building and linking.
   - macOS dynamic library naming and `DYLD_LIBRARY_PATH` or install-name notes.
   - Linux dynamic library naming and `LD_LIBRARY_PATH` notes.
   - Cargo build command for the shared library.
   - Header include path.

5. Document ownership and errors.
   - Which functions allocate.
   - Which functions free.
   - Whether returned strings/buffers need explicit release.
   - How status codes and error messages work.

6. Compile examples in CI.
   - Add a C smoke example.
   - Add a C++ smoke example if a C++ compiler is available.
   - Keep examples small enough to serve as documentation.

## Tests and verification

- Examples compile and run in CI.
- README instructions are tested manually on macOS and Linux.
- Header and library paths in docs match actual build artifacts.

## Risks

- Documentation can become stale unless examples are compiled in CI.
- Platform linker details are easy to get wrong.
- A C++ wrapper should not hide ownership problems before the C ABI is stable.
