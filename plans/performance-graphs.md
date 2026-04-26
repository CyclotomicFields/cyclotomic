# Performance Graphs and Benchmarks

## Goal

Replace stale benchmark references with reproducible benchmarks and graphs that
show when each representation is useful.

## Current findings

- README still says graphs are TODO.
- The benchmark scripts reference the removed `cyclobench` binary and old ANTIC
  comparison paths.
- There is no current `benches/` Criterion harness.

## Plan

1. Remove or update stale benchmark docs.
   - Stop referring to the removed ANTIC benchmark binary unless a new optional
     comparison tool is added.
   - Make README language match the current dependency story.

2. Add Criterion benchmarks.
   - Sparse add/mul.
   - Dense add/mul.
   - Structure-constant precomputation and multiplication.
   - Scalar multiplication.
   - Basis conversion/order promotion where relevant.

3. Parameterize workloads.
   - Field order.
   - `phi(n)`.
   - Input sparsity.
   - Coefficient size.
   - Repeated multiplication versus one-off construction.

4. Export data.
   - Use Criterion reports for local inspection.
   - Export CSV or JSON for custom plotting.
   - Keep generated graph scripts in `benchmarks/`.

5. Generate README-ready graphs.
   - Plot representation comparison by order and sparsity.
   - Include enough metadata to reproduce the graph.
   - Either check in static images or publish CI artifacts.

6. Add a benchmark maintenance workflow.
   - A normal CI job should run tests only.
   - A manual or scheduled workflow can run benchmarks.
   - Document expected local runtime.

## Tests and verification

- `cargo bench` or documented benchmark command works locally.
- Plot script works from generated data.
- README graph paths resolve.

## Risks

- Benchmarks can be noisy and machine-dependent.
- Long benchmark jobs should not slow down normal CI.
- Historical ANTIC comparisons require new optional infrastructure if still
  desired.
