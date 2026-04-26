# Plans

This directory contains implementation plans for the TODO items in the
project README. Each plan is intentionally scoped so it can become one or
more focused pull requests.

## TODO plans

1. [Remove clones in rational implementations](remove-rational-clones.md)
2. [Use f64 for integer-valued paths](f64-integer-paths.md)
3. [Make structure constants competitive with floats](float-structure-constants.md)
4. [Dense polynomial arithmetic modulo cyclotomic polynomials](dense-polynomial-modulo.md)
5. [Reduce copy/paste with useful abstractions](reduce-copy-paste.md)
6. [C FFI interface](c-ffi-interface.md)
7. [Field of fractions of the ring of integers](field-of-fractions.md)
8. [SIMD investigation](simd-investigation.md)
9. [Performance graphs and benchmarks](performance-graphs.md)
10. [C and C++ usage instructions](c-cpp-usage-instructions.md)

## Suggested order

1. Build a real benchmark harness and fix stale benchmark documentation.
2. Refactor the rational/algebra API enough to remove forced clones.
3. Use the benchmark harness to evaluate dense modulo arithmetic,
   structure constants, float-backed experiments, and SIMD.
4. Add the C ABI once the core Rust API shape is less likely to churn.
5. Write C/C++ usage instructions against the tested ABI.
