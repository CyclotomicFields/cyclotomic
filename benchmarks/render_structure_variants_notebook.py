#!/usr/bin/env python3
"""Render structure-constant variant benchmark results into a notebook."""

from __future__ import annotations

import datetime as _datetime
import json
import platform
from pathlib import Path

from render_notebook import code_cell, fmt_ops, line_chart, markdown_cell, ops_per_second, series_by


ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = ROOT / "benchmarks" / "results" / "structure_variants_quick.json"
NOTEBOOK_PATH = ROOT / "notebooks" / "structure_constant_variants.ipynb"


def load_records() -> list[dict[str, object]]:
    with RESULTS_PATH.open() as f:
        return json.load(f)


def structure_mul_records(records: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        record
        for record in records
        if record["operation"] == "mul" and str(record["representation"]).startswith("structure")
    ]


def table_text(records: list[dict[str, object]], density: float) -> str:
    rows = [
        record
        for record in structure_mul_records(records)
        if abs(float(record["density"]) - density) < 1e-9
    ]
    rows.sort(key=lambda r: (int(r["order"]), str(r["representation"])))

    lines = [
        "| order | representation | ops/s | speedup vs nested |",
        "|---:|---|---:|---:|",
    ]
    nested = {
        int(record["order"]): ops_per_second(record)
        for record in rows
        if record["representation"] == "structure"
    }

    for record in rows:
        order = int(record["order"])
        ops = ops_per_second(record)
        speedup = ops / nested[order]
        lines.append(
            f"| {order} | `{record['representation']}` | {fmt_ops(ops)} | {speedup:.2f}x |"
        )

    return "\n".join(lines)


def summary_text(records: list[dict[str, object]]) -> str:
    structure_records = structure_mul_records(records)
    lines = [
        f"records: {len(records)}",
        f"structure multiplication records: {len(structure_records)}",
        f"orders: {sorted({record['order'] for record in records})}",
        f"generated_at_utc: {_datetime.datetime.now(_datetime.UTC).replace(microsecond=0).isoformat()}",
        "",
    ]

    for order in sorted({int(record["order"]) for record in structure_records}):
        full_density = [
            record
            for record in structure_records
            if int(record["order"]) == order and abs(float(record["density"]) - 1.0) < 1e-9
        ]
        best = max(full_density, key=ops_per_second)
        nested = next(record for record in full_density if record["representation"] == "structure")
        lines.append(
            f"n={order}: best full-density variant is {best['representation']} "
            f"at {fmt_ops(ops_per_second(best))} ops/s "
            f"({ops_per_second(best) / ops_per_second(nested):.2f}x nested)"
        )

    return "\n".join(lines) + "\n"


def build_notebook(records: list[dict[str, object]]) -> dict[str, object]:
    structure_records = structure_mul_records(records)
    charts = [
        (
            "Full-density throughput",
            line_chart(
                "Structure-constant multiplication at full density",
                "Higher is better; log-scale operations per second",
                "cyclotomic order n",
                "operations per second",
                series_by(
                    [
                        record
                        for record in structure_records
                        if abs(float(record["density"]) - 1.0) < 1e-9
                    ],
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    ops_per_second,
                ),
            ),
        ),
        (
            "Density sensitivity",
            line_chart(
                "Structure variants by density for order 12",
                "Order 12 is where structure constants were already competitive",
                "input density",
                "operations per second",
                series_by(
                    [
                        record
                        for record in structure_records
                        if int(record["order"]) == 12
                    ],
                    lambda r: str(r["representation"]),
                    lambda r: float(r["density"]),
                    ops_per_second,
                ),
            ),
        ),
        (
            "Best sparse variant by density",
            line_chart(
                "Sparse structure variant with zero-input skipping",
                "This variant benefits most from sparse inputs",
                "input density",
                "operations per second",
                series_by(
                    [
                        record
                        for record in structure_records
                        if record["representation"] == "structure_sparse_skip_zero"
                    ],
                    lambda r: f"n={r['order']}",
                    lambda r: float(r["density"]),
                    ops_per_second,
                ),
            ),
        ),
    ]

    cells: list[dict[str, object]] = [
        markdown_cell(
            "# Structure-Constant Layout Variant Benchmarks\n\n"
            "This notebook contains the quick benchmark run for the exact "
            "structure-constant layout variants introduced in "
            "`CyclotomicField`. The committed outputs are viewable in GitHub."
        ),
        code_cell(
            "import json\n"
            "from pathlib import Path\n\n"
            "root = Path.cwd()\n"
            "if not (root / 'benchmarks').exists():\n"
            "    root = root.parent\n"
            "records = json.loads((root / 'benchmarks/results/structure_variants_quick.json').read_text())\n"
            "len(records), records[0]",
            output_text=summary_text(records),
            count=1,
        ),
        markdown_cell(
            "## Reproduction command\n\n"
            "```bash\n"
            "cargo run --release --example perf_bench -- --quick > benchmarks/results/structure_variants_quick.json\n"
            "python3 benchmarks/render_structure_variants_notebook.py\n"
            "```"
        ),
        markdown_cell("## Full-Density Table\n\n" + table_text(records, 1.0)),
        markdown_cell("## Half-Density Table\n\n" + table_text(records, 0.5)),
    ]

    for i, (heading, svg) in enumerate(charts, start=2):
        cells.append(markdown_cell(f"## {heading}"))
        cells.append(code_cell("# Rendered from benchmarks/results/structure_variants_quick.json", svg=svg, count=i))

    cells.append(
        markdown_cell(
            "## Notes\n\n"
            "- `structure` is the original nested `Vec<Vec<Vec<Q>>>` layout.\n"
            "- `structure_flat` uses one contiguous tensor.\n"
            "- `structure_sparse` stores only nonzero constants.\n"
            "- `structure_sparse_skip_zero` also skips zero input coefficients.\n"
            "- These are exact `rug::Rational` benchmarks, not approximate float benchmarks."
        )
    )

    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "pygments_lexer": "ipython3",
                "version": platform.python_version(),
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main() -> None:
    records = load_records()
    NOTEBOOK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with NOTEBOOK_PATH.open("w") as f:
        json.dump(build_notebook(records), f, indent=1)
        f.write("\n")
    print(f"wrote {NOTEBOOK_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
