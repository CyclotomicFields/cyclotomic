# GAP cyclotomic reference provenance

This directory contains a standalone derivative of GAP's cyclotomic kernel
arithmetic for differential tests and like-for-like benchmarks.

- Upstream: <https://github.com/gap-system/gap>
- Pinned commit: `d2134de71521c62512b8351c42ec16bfbac21744`
- Primary source: `src/cyclotom.c`
- Upstream license: GPL-2.0-or-later

The C implementation retains GAP's `ConvertToBase`, `Cyclotomic`,
`FindCommonField`, addition, and multiplication algorithms and function
decomposition. GAP's bags, garbage collector, generic `Obj` arithmetic, error
handling, and global `ResultCyc` buffer have been replaced by:

- explicit heap-owned packed elements;
- a context-owned scratch buffer;
- checked `int64_t` coefficients;
- a narrow C ABI used by the Rust ownership wrapper.

This first milestone intentionally supports only small integer coefficients.
GMP integer/rational promotion, inversion, Galois actions, and an unmodified
libgap oracle remain follow-up work.

Because this is derived from GPL-licensed GAP code, this reference crate is
isolated from the normal LGPL library build and is not published.
