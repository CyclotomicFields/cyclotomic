#!/usr/bin/env python3
"""Render benchmark results into an executed Jupyter notebook.

This script intentionally uses only the Python standard library so the
checked-in notebook can be regenerated on a minimal developer machine after
running:

    cargo run --release --example perf_bench > benchmarks/results/performance.json
    python3 benchmarks/render_notebook.py
"""

from __future__ import annotations

import datetime as _datetime
import json
import math
import platform
from pathlib import Path
from typing import Callable, Iterable


ROOT = Path(__file__).resolve().parents[1]
RESULTS_PATH = ROOT / "benchmarks" / "results" / "performance.json"
NOTEBOOK_PATH = ROOT / "notebooks" / "performance_benchmarks.ipynb"


Record = dict[str, object]


def load_records() -> list[Record]:
    with RESULTS_PATH.open() as f:
        return json.load(f)


def median(values: Iterable[float]) -> float:
    xs = sorted(values)
    if not xs:
        return math.nan
    mid = len(xs) // 2
    if len(xs) % 2:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2.0


def fmt_ns(ns: float) -> str:
    if ns >= 1_000_000:
        return f"{ns / 1_000_000:.2f} ms"
    if ns >= 1_000:
        return f"{ns / 1_000:.2f} us"
    return f"{ns:.0f} ns"


def select(
    records: list[Record],
    *,
    operation: str,
    density: float | None = None,
    representations: set[str] | None = None,
    max_order: int | None = None,
) -> list[Record]:
    selected = []
    for record in records:
        if record["operation"] != operation:
            continue
        if density is not None and abs(float(record["density"]) - density) > 1e-9:
            continue
        if representations is not None and record["representation"] not in representations:
            continue
        if max_order is not None and int(record["order"]) > max_order:
            continue
        selected.append(record)
    return selected


def series_by(
    records: list[Record],
    key_fn: Callable[[Record], str],
    x_fn: Callable[[Record], float],
    y_fn: Callable[[Record], float],
) -> dict[str, list[tuple[float, float]]]:
    series: dict[str, list[tuple[float, float]]] = {}
    for record in records:
        series.setdefault(key_fn(record), []).append((x_fn(record), y_fn(record)))
    return {k: sorted(v) for k, v in sorted(series.items())}


def line_chart(
    title: str,
    subtitle: str,
    x_label: str,
    y_label: str,
    series: dict[str, list[tuple[float, float]]],
    *,
    width: int = 860,
    height: int = 470,
    log_y: bool = True,
) -> str:
    margin = {"left": 88, "right": 24, "top": 72, "bottom": 64}
    plot_w = width - margin["left"] - margin["right"]
    plot_h = height - margin["top"] - margin["bottom"]
    colors = ["#2563eb", "#dc2626", "#059669", "#7c3aed", "#ea580c", "#0891b2"]

    points = [pt for values in series.values() for pt in values]
    xs = [x for x, _y in points]
    raw_ys = [max(y, 1e-9) for _x, y in points]
    min_x, max_x = min(xs), max(xs)
    if min_x == max_x:
        max_x = min_x + 1.0
    if log_y:
        ys = [math.log10(y) for y in raw_ys]
        min_y, max_y = min(ys), max(ys)
        y_transform = lambda y: math.log10(max(y, 1e-9))
        y_tick_values = [10**v for v in range(math.floor(min_y), math.ceil(max_y) + 1)]
    else:
        min_y, max_y = min(raw_ys), max(raw_ys)
        y_transform = lambda y: y
        span = max_y - min_y or 1.0
        y_tick_values = [min_y + span * i / 4 for i in range(5)]
    if min_y == max_y:
        max_y = min_y + 1.0

    def sx(x: float) -> float:
        return margin["left"] + (x - min_x) / (max_x - min_x) * plot_w

    def sy(y: float) -> float:
        yv = y_transform(y)
        return margin["top"] + (max_y - yv) / (max_y - min_y) * plot_h

    x_ticks = sorted(set(xs))
    if len(x_ticks) > 8:
        step = max(1, len(x_ticks) // 6)
        x_ticks = x_ticks[::step]

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin["left"]}" y="30" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#111827">{title}</text>',
        f'<text x="{margin["left"]}" y="52" font-family="Arial, sans-serif" font-size="13" fill="#4b5563">{subtitle}</text>',
        f'<line x1="{margin["left"]}" y1="{margin["top"] + plot_h}" x2="{margin["left"] + plot_w}" y2="{margin["top"] + plot_h}" stroke="#111827" stroke-width="1"/>',
        f'<line x1="{margin["left"]}" y1="{margin["top"]}" x2="{margin["left"]}" y2="{margin["top"] + plot_h}" stroke="#111827" stroke-width="1"/>',
    ]

    for tick in x_ticks:
        x = sx(tick)
        svg.append(f'<line x1="{x:.1f}" y1="{margin["top"]}" x2="{x:.1f}" y2="{margin["top"] + plot_h}" stroke="#e5e7eb" stroke-width="1"/>')
        svg.append(f'<text x="{x:.1f}" y="{height - 38}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#374151">{tick:g}</text>')

    for tick in y_tick_values:
        if tick <= 0:
            continue
        y = sy(tick)
        svg.append(f'<line x1="{margin["left"]}" y1="{y:.1f}" x2="{margin["left"] + plot_w}" y2="{y:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        svg.append(f'<text x="{margin["left"] - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial, sans-serif" font-size="11" fill="#374151">{fmt_ns(tick)}</text>')

    for i, (name, values) in enumerate(series.items()):
        color = colors[i % len(colors)]
        path = " ".join(
            f"{'M' if j == 0 else 'L'} {sx(x):.1f} {sy(y):.1f}" for j, (x, y) in enumerate(values)
        )
        svg.append(f'<path d="{path}" fill="none" stroke="{color}" stroke-width="2.4"/>')
        for x, y in values:
            svg.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.2" fill="{color}"/>')
        legend_x = margin["left"] + 10 + (i % 3) * 210
        legend_y = height - 18 - (i // 3) * 18
        svg.append(f'<line x1="{legend_x}" y1="{legend_y - 4}" x2="{legend_x + 22}" y2="{legend_y - 4}" stroke="{color}" stroke-width="3"/>')
        svg.append(f'<text x="{legend_x + 28}" y="{legend_y}" font-family="Arial, sans-serif" font-size="12" fill="#111827">{name}</text>')

    svg.append(f'<text x="{margin["left"] + plot_w / 2}" y="{height - 10}" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#111827">{x_label}</text>')
    svg.append(f'<text transform="translate(18 {margin["top"] + plot_h / 2}) rotate(-90)" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#111827">{y_label}</text>')
    svg.append("</svg>")
    return "\n".join(svg)


def markdown_cell(source: str) -> dict[str, object]:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code_cell(source: str, *, output_text: str | None = None, svg: str | None = None, count: int = 1) -> dict[str, object]:
    outputs = []
    if output_text is not None:
        outputs.append({"name": "stdout", "output_type": "stream", "text": output_text.splitlines(keepends=True)})
    if svg is not None:
        outputs.append({"output_type": "display_data", "metadata": {}, "data": {"image/svg+xml": svg}})
    return {
        "cell_type": "code",
        "execution_count": count,
        "metadata": {},
        "outputs": outputs,
        "source": source.splitlines(keepends=True),
    }


def summary_text(records: list[Record]) -> str:
    lines = [
        f"records: {len(records)}",
        f"python: {platform.python_version()}",
        f"platform: {platform.platform()}",
        f"generated_at_utc: {_datetime.datetime.now(_datetime.UTC).replace(microsecond=0).isoformat()}",
        "",
        "median ns/iter by representation and operation:",
    ]
    keys = sorted({(r["representation"], r["operation"]) for r in records})
    for representation, operation in keys:
        values = [
            float(r["ns_per_iter"])
            for r in records
            if r["representation"] == representation and r["operation"] == operation
        ]
        lines.append(f"  {representation:9s} {operation:11s} {fmt_ns(median(values))}")
    return "\n".join(lines) + "\n"


def build_notebook(records: list[Record]) -> dict[str, object]:
    charts = [
        (
            "Dense vs sparse multiplication, low density",
            line_chart(
                "Multiplication at 35% density",
                "Lower is better; exact rug::Rational coefficients",
                "cyclotomic order n",
                "time per operation, log scale",
                series_by(
                    select(records, operation="mul", density=0.35, representations={"dense", "sparse"}),
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Dense vs sparse multiplication, full density",
            line_chart(
                "Multiplication at 100% density",
                "Dense currently pays O(n^2) over the full order; sparse pays map overhead",
                "cyclotomic order n",
                "time per operation, log scale",
                series_by(
                    select(records, operation="mul", density=1.0, representations={"dense", "sparse"}),
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Addition by order",
            line_chart(
                "Addition at 100% density",
                "Dense vector addition versus sparse map insertion/update",
                "cyclotomic order n",
                "time per operation, log scale",
                series_by(
                    select(records, operation="add", density=1.0, representations={"dense", "sparse"}),
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Scalar multiplication by order",
            line_chart(
                "Scalar multiplication at 100% density",
                "Dense mutates all coefficients; sparse iterates stored terms",
                "cyclotomic order n",
                "time per operation, log scale",
                series_by(
                    select(records, operation="scalar_mul", density=1.0, representations={"dense", "sparse"}),
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Sparse multiplication by density",
            line_chart(
                "Sparse multiplication sensitivity to density",
                "Each line fixes the field order",
                "input density",
                "time per operation, log scale",
                series_by(
                    [
                        r
                        for r in select(records, operation="mul", representations={"sparse"})
                        if r["order"] in {8, 12, 20, 24}
                    ],
                    lambda r: f"n={r['order']}",
                    lambda r: float(r["density"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Dense multiplication by density",
            line_chart(
                "Dense multiplication sensitivity to density",
                "The current dense loop still scans the whole coefficient vector",
                "input density",
                "time per operation, log scale",
                series_by(
                    [
                        r
                        for r in select(records, operation="mul", representations={"dense"})
                        if r["order"] in {8, 12, 20, 24}
                    ],
                    lambda r: f"n={r['order']}",
                    lambda r: float(r["density"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Structure constant construction",
            line_chart(
                "Structure-constant precomputation cost",
                "Construction is measured separately from steady-state multiplication",
                "phi(n)",
                "construction time, log scale",
                series_by(
                    select(records, operation="construct", representations={"structure"}),
                    lambda r: "structure",
                    lambda r: float(r["phi"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
        (
            "Structure multiplication comparison",
            line_chart(
                "Multiplication at 100% density, common orders",
                "Structure constants are only measured for smaller orders in this sweep",
                "cyclotomic order n",
                "time per operation, log scale",
                series_by(
                    select(records, operation="mul", density=1.0, max_order=16),
                    lambda r: str(r["representation"]),
                    lambda r: float(r["order"]),
                    lambda r: float(r["ns_per_iter"]),
                ),
            ),
        ),
    ]

    cells: list[dict[str, object]] = [
        markdown_cell(
            "# Cyclotomic Performance Benchmarks\n\n"
            "This notebook is a checked-in benchmark report generated from the Rust "
            "`examples/perf_bench.rs` executable. The committed outputs make the "
            "graphs viewable in GitHub's notebook renderer without rerunning the "
            "benchmarks.\n\n"
            "The benchmark is intentionally lightweight and comparative. It is useful "
            "for spotting representation-level trends, not for publishing stable "
            "machine-independent numbers."
        ),
        code_cell(
            "import json\n"
            "from pathlib import Path\n\n"
            "root = Path.cwd()\n"
            "if not (root / 'benchmarks').exists():\n"
            "    root = root.parent\n"
            "records = json.loads((root / 'benchmarks/results/performance.json').read_text())\n"
            "len(records), records[0]",
            output_text=summary_text(records),
            count=1,
        ),
        markdown_cell(
            "## Reproduction command\n\n"
            "Run this from the repository root to refresh the data and notebook:\n\n"
            "```bash\n"
            "cargo run --release --example perf_bench > benchmarks/results/performance.json\n"
            "python3 benchmarks/render_notebook.py\n"
            "```"
        ),
    ]

    for i, (heading, svg) in enumerate(charts, start=2):
        cells.append(markdown_cell(f"## {heading}"))
        cells.append(code_cell("# Rendered from benchmarks/results/performance.json", svg=svg, count=i))

    cells.append(
        markdown_cell(
            "## Notes\n\n"
            "- Dense and sparse operations include clone/setup cost for each timed "
            "operation so the benchmark reflects the current mutating API shape.\n"
            "- Structure-constant construction is separated from multiplication "
            "because it is a one-time field setup cost.\n"
            "- The current results use exact `rug::Rational` coefficients."
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
                "codemirror_mode": {"name": "ipython", "version": 3},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
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
