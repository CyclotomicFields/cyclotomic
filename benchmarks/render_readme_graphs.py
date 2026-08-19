#!/usr/bin/env python3
"""Render the checked-in benchmark data as README-ready SVG graphs.

The script uses only the Python standard library. Run it from any directory:

    python3 benchmarks/render_readme_graphs.py
"""

from __future__ import annotations

import html
import json
from pathlib import Path

from render_notebook import line_chart, ops_per_second, select, series_by


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "benchmarks" / "results"
GRAPH_DIR = ROOT / "benchmarks" / "graphs"


def load_json(name: str) -> list[dict[str, object]]:
    with (RESULTS_DIR / name).open() as source:
        return json.load(source)


def gap_speedup_chart(records: list[dict[str, object]]) -> str:
    width, height = 760, 430
    margin = {"left": 78, "right": 28, "top": 72, "bottom": 72}
    plot_width = width - margin["left"] - margin["right"]
    plot_height = height - margin["top"] - margin["bottom"]
    maximum = max(float(record["speedup"]) for record in records) * 1.15
    slot_width = plot_width / len(records)
    bar_width = slot_width * 0.55

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{margin["left"]}" y="30" font-family="Arial, sans-serif" font-size="20" font-weight="700" fill="#111827">Character-table tensor decomposition</text>',
        f'<text x="{margin["left"]}" y="52" font-family="Arial, sans-serif" font-size="13" fill="#4b5563">Rust arithmetic speedup over GAP; higher is better</text>',
    ]

    for tick in range(0, 8):
        y = margin["top"] + plot_height - tick / maximum * plot_height
        svg.append(
            f'<line x1="{margin["left"]}" y1="{y:.1f}" x2="{margin["left"] + plot_width}" y2="{y:.1f}" stroke="#e5e7eb"/>'
        )
        svg.append(
            f'<text x="{margin["left"] - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial, sans-serif" font-size="11" fill="#374151">{tick}x</text>'
        )

    for index, record in enumerate(records):
        speedup = float(record["speedup"])
        x = margin["left"] + index * slot_width + (slot_width - bar_width) / 2
        bar_height = speedup / maximum * plot_height
        y = margin["top"] + plot_height - bar_height
        label = html.escape(str(record["table"]))
        svg.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" rx="4" fill="#2563eb"/>'
        )
        svg.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{y - 9:.1f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" font-weight="700" fill="#1f2937">{speedup:.2f}x</text>'
        )
        svg.append(
            f'<text x="{x + bar_width / 2:.1f}" y="{height - 42}" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#111827">{label}</text>'
        )

    svg.extend(
        [
            f'<line x1="{margin["left"]}" y1="{margin["top"] + plot_height}" x2="{margin["left"] + plot_width}" y2="{margin["top"] + plot_height}" stroke="#111827"/>',
            f'<line x1="{margin["left"]}" y1="{margin["top"]}" x2="{margin["left"]}" y2="{margin["top"] + plot_height}" stroke="#111827"/>',
            f'<text transform="translate(18 {margin["top"] + plot_height / 2}) rotate(-90)" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" fill="#111827">speedup over GAP</text>',
            "</svg>",
        ]
    )
    return "\n".join(svg)


def main() -> None:
    core = load_json("performance.json")
    structure = load_json("structure_variants_quick.json")
    gap = load_json("gap_character_tables.json")
    GRAPH_DIR.mkdir(parents=True, exist_ok=True)

    graphs = {
        "dense-vs-sparse-multiplication.svg": line_chart(
            "Multiplication throughput at 25% density",
            "Exact rug::Rational coefficients; log scale; higher is better",
            "cyclotomic order n",
            "operations per second",
            series_by(
                select(
                    core,
                    operation="mul",
                    density=0.25,
                    representations={"dense", "sparse"},
                ),
                lambda record: str(record["representation"]),
                lambda record: float(record["order"]),
                ops_per_second,
            ),
        ),
        "structure-layout-multiplication.svg": line_chart(
            "Structure-constant multiplication at full density",
            "Exact rug::Rational coefficients; log scale; higher is better",
            "cyclotomic order n",
            "operations per second",
            series_by(
                [
                    record
                    for record in structure
                    if record["operation"] == "mul"
                    and abs(float(record["density"]) - 1.0) < 1e-9
                    and str(record["representation"]).startswith("structure")
                ],
                lambda record: str(record["representation"]),
                lambda record: float(record["order"]),
                ops_per_second,
            ),
        ),
        "gap-character-table-speedup.svg": gap_speedup_chart(gap),
    }

    for name, contents in graphs.items():
        path = GRAPH_DIR / name
        path.write_text(contents + "\n")
        print(f"wrote {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
