#!/usr/bin/env python
"""Render spatial gradient backend benchmark metrics as a funky heatmap.

The script prefers ``omicverse.pl.funky_heatmap`` when OmicVerse is installed,
falls back to ``pyfunkyheatmap`` when available, and finally renders a compact
matplotlib version so benchmark artifacts can still be created locally.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd


METRIC_ORDER = [
    "Runtime",
    "Speedup",
    "Schema",
    "No S4 object",
    "TopK Jaccard",
    "Rank rho",
    "FDR rho",
    "P-value rho",
    "Curve rho",
    "Criteria",
]

METRIC_GROUPS = {
    "Runtime": "Runtime",
    "Speedup": "Runtime",
    "Schema": "Integrity",
    "No S4 object": "Integrity",
    "TopK Jaccard": "Agreement",
    "Rank rho": "Agreement",
    "FDR rho": "Agreement",
    "P-value rho": "Agreement",
    "Curve rho": "Agreement",
    "Criteria": "Criteria",
}

GLYPH_TYPES = {
    "Runtime": "bar",
    "Speedup": "bar",
    "Schema": "funkyrect",
    "No S4 object": "funkyrect",
    "TopK Jaccard": "circle",
    "Rank rho": "circle",
    "FDR rho": "circle",
    "P-value rho": "circle",
    "Curve rho": "circle",
    "Criteria": "pie",
}


def _load_funky_heatmap() -> tuple[Callable[..., Any] | None, str]:
    try:
        import omicverse as ov  # type: ignore

        if hasattr(ov, "ov_plot_set"):
            ov.ov_plot_set()
        elif hasattr(ov, "plot_set"):
            ov.plot_set()
        funky = getattr(ov.pl, "funky_heatmap")
        return funky, "omicverse.pl.funky_heatmap"
    except Exception as ov_error:
        try:
            from pyfunkyheatmap import funky_heatmap  # type: ignore

            return funky_heatmap, "pyfunkyheatmap.funky_heatmap"
        except Exception as py_error:
            return None, f"matplotlib fallback; ov={ov_error}; pyfunkyheatmap={py_error}"


def _format_number(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def _metric_slug(metric: str) -> str:
    return (
        metric.lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
        .replace("+", "plus")
    )


def _case_label(case: str) -> str:
    match = re.match(r"^(?P<reference>.+)_(?P<genes>\d+)genes_(?P<spots>\d+)spots$", case)
    if not match:
        return case
    return f"{match.group('reference')}\n{match.group('genes')} genes, {match.group('spots')} spots"


def _score_to_float(value: Any) -> float:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return np.nan
    if not np.isfinite(value):
        return np.nan
    return max(0.0, min(1.0, value))


def prepare_heatmap_data(metrics_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(metrics_path)
    required = {"case", "method", "metric", "group", "score", "label", "status"}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError(f"Missing required metrics columns: {sorted(missing)}")

    metrics["score"] = pd.to_numeric(metrics["score"], errors="coerce")
    metrics["label"] = metrics["label"].fillna("NA").astype(str)
    rows: list[dict[str, Any]] = []

    for (case, method), group in metrics.groupby(["case", "method"], sort=False):
        row: dict[str, Any] = {
            "id": f"{case} | {method}",
            "case": case,
            "case_label": _case_label(str(case)),
            "method": method,
        }
        statuses = []
        for metric in METRIC_ORDER:
            if metric == "Criteria":
                continue
            hit = group[group["metric"] == metric]
            if hit.empty:
                score = np.nan
                label = "NA"
                status = "na"
            else:
                item = hit.iloc[0]
                score = _score_to_float(item["score"])
                label = str(item["label"])
                status = str(item["status"])
            slug = _metric_slug(metric)
            row[f"{slug}_score"] = 0.0 if pd.isna(score) else score
            row[f"{slug}_label"] = label
            row[f"{slug}_status"] = status
            statuses.append(status)

        pass_count = sum(status in {"ok", "pass", "reference"} for status in statuses)
        watch_count = max(len(statuses) - pass_count, 0)
        row["criteria_pie"] = {"pass": int(pass_count), "watch": max(int(watch_count), 1e-6)}
        row["criteria_label"] = f"{pass_count}/{len(statuses)}"
        row["criteria_score"] = pass_count / max(len(statuses), 1)
        rows.append(row)

    heatmap_data = pd.DataFrame(rows)
    metric_info_rows = []
    for metric in METRIC_ORDER:
        slug = _metric_slug(metric)
        if metric == "Criteria":
            metric_info_rows.append(
                {
                    "metric": metric,
                    "score_col": "criteria_pie",
                    "label_col": "criteria_label",
                    "geom": "pie",
                    "group": "Criteria",
                }
            )
        else:
            metric_info_rows.append(
                {
                    "metric": metric,
                    "score_col": f"{slug}_score",
                    "label_col": f"{slug}_label",
                    "geom": GLYPH_TYPES[metric],
                    "group": METRIC_GROUPS[metric],
                }
            )
    return heatmap_data, pd.DataFrame(metric_info_rows)


def _column_info(metric_info: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = [
        {"id": "method", "name": "Method", "geom": "text", "group": "Design", "width": 3.2, "legend": False}
    ]
    for item in metric_info.to_dict("records"):
        rows.append(
            {
                "id": item["score_col"],
                "name": item["metric"],
                "geom": item["geom"],
                "group": item["group"],
                "palette": "criteria" if item["geom"] == "pie" else item["group"].lower(),
                "width": 2.2 if item["geom"] == "bar" else 1.8,
                "legend": item["geom"] != "text",
            }
        )
        rows.append(
            {
                "id": item["label_col"],
                "name": "",
                "geom": "text",
                "group": item["group"],
                "overlay": True,
                "width": 0.0,
                "legend": False,
            }
        )
    return pd.DataFrame(rows)


def _render_with_funky(
    funky_heatmap: Callable[..., Any],
    heatmap_data: pd.DataFrame,
    metric_info: pd.DataFrame,
    output_path: Path,
) -> None:
    column_info = _column_info(metric_info)
    group_names = ["Design", "Runtime", "Integrity", "Agreement", "Criteria"]
    column_groups = pd.DataFrame({"group": group_names, "level1": group_names})
    row_info = pd.DataFrame({"id": heatmap_data["id"], "group": heatmap_data["case"]})
    cases = list(dict.fromkeys(heatmap_data["case"]))
    row_groups = pd.DataFrame({"group": cases, "level1": cases})
    palettes = {
        "runtime": ["#f7fbff", "#c6dbef", "#6baed6", "#08519c"],
        "integrity": ["#f7fcf5", "#bae4b3", "#74c476", "#238b45"],
        "agreement": ["#fff7bc", "#fec44f", "#fd8d3c", "#d95f0e"],
        "criteria": {"pass": "#1b9e77", "watch": "#d9d9d9"},
    }
    legends = [
        {"title": "Runtime score", "palette": "runtime", "geom": "bar"},
        {"title": "Integrity", "palette": "integrity", "geom": "funkyrect"},
        {"title": "Agreement", "palette": "agreement", "geom": "circle"},
        {"title": "Criteria", "palette": "criteria", "geom": "pie"},
    ]
    figure = funky_heatmap(
        heatmap_data,
        column_info=column_info,
        column_groups=column_groups,
        row_info=row_info,
        row_groups=row_groups,
        palettes=palettes,
        legends=legends,
        position_args={
            "expand_xmin": 2.0,
            "expand_xmax": 2.6,
            "col_annot_offset": 3.4,
            "col_annot_angle": 35,
            "col_bigspace": 0.8,
            "cell_text_size": 3.2,
            "col_annot_size": 3.8,
        },
        scale_column=False,
        add_abc=False,
        fig_scale=0.30,
        dpi=200,
    )
    if hasattr(figure, "save"):
        figure.save(output_path, facecolor="white")
    elif hasattr(figure, "figure"):
        figure.figure.savefig(output_path, facecolor="white")
    else:
        figure.savefig(output_path, facecolor="white")


def _render_matplotlib(
    heatmap_data: pd.DataFrame,
    metric_info: pd.DataFrame,
    output_path: Path,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Rectangle, Wedge

    metrics = metric_info.to_dict("records")
    n_rows = len(heatmap_data)
    n_cols = len(metrics)
    width = max(11.0, 2.0 + n_cols * 1.35)
    height = max(3.2, 1.2 + n_rows * 0.65)
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlim(-2.6, n_cols - 0.35)
    ax.set_ylim(-0.8, n_rows + 0.9)
    ax.axis("off")
    cmap = plt.get_cmap("YlOrRd")
    runtime_cmap = plt.get_cmap("Blues")
    integrity_cmap = plt.get_cmap("Greens")

    def color_for(score: float, group: str) -> Any:
        if not np.isfinite(score):
            return "#d9d9d9"
        score = max(0.0, min(1.0, score))
        if group == "Runtime":
            return runtime_cmap(0.18 + 0.72 * score)
        if group == "Integrity":
            return integrity_cmap(0.18 + 0.72 * score)
        return cmap(0.18 + 0.72 * score)

    for col, item in enumerate(metrics):
        ax.text(col, n_rows + 0.16, item["metric"], ha="right", va="bottom", rotation=38, fontsize=9)
        ax.plot([col - 0.47, col + 0.47], [n_rows - 0.05, n_rows - 0.05], color="#d0d0d0", lw=0.7)

    for row_idx, row in heatmap_data.reset_index(drop=True).iterrows():
        y = n_rows - row_idx - 1
        ax.text(-3.55, y, str(row["case_label"]), ha="left", va="center", fontsize=8, color="#555555", linespacing=1.15)
        ax.text(-1.35, y, str(row["method"]), ha="left", va="center", fontsize=9, color="#111111")
        ax.plot([-3.65, n_cols - 0.42], [y - 0.38, y - 0.38], color="#eeeeee", lw=0.7)
        for col, item in enumerate(metrics):
            score_col = item["score_col"]
            label_col = item["label_col"]
            geom = item["geom"]
            group = item["group"]
            label = str(row.get(label_col, "NA"))
            if geom == "pie":
                pie = row.get(score_col, {"pass": 0, "watch": 1})
                total = float(pie.get("pass", 0)) + float(pie.get("watch", 0))
                frac = 0 if total <= 0 else float(pie.get("pass", 0)) / total
                ax.add_patch(Wedge((col, y), 0.32, 90, 90 + 360 * frac, facecolor="#1b9e77", edgecolor="white", lw=0.6))
                ax.add_patch(Wedge((col, y), 0.32, 90 + 360 * frac, 450, facecolor="#d9d9d9", edgecolor="white", lw=0.6))
                ax.add_patch(Circle((col, y), 0.32, fill=False, edgecolor="#777777", lw=0.5))
            else:
                score = _score_to_float(row.get(score_col, np.nan))
                if geom == "bar":
                    ax.add_patch(Rectangle((col - 0.42, y - 0.16), 0.84, 0.32, facecolor="#eeeeee", edgecolor="none"))
                    if np.isfinite(score):
                        ax.add_patch(Rectangle((col - 0.42, y - 0.16), 0.84 * score, 0.32, facecolor=color_for(score, group), edgecolor="none"))
                    ax.add_patch(Rectangle((col - 0.42, y - 0.16), 0.84, 0.32, fill=False, edgecolor="#777777", lw=0.5))
                elif geom == "funkyrect":
                    size = 0.16 + 0.46 * (score if np.isfinite(score) else 0)
                    ax.add_patch(Rectangle((col - size / 2, y - size / 2), size, size, facecolor=color_for(score, group), edgecolor="#777777", lw=0.5))
                else:
                    radius = 0.11 + 0.25 * (score if np.isfinite(score) else 0)
                    ax.add_patch(Circle((col, y), radius, facecolor=color_for(score, group), edgecolor="#777777", lw=0.5))
            ax.text(col, y, label, ha="center", va="center", fontsize=7.2, color="#111111")

    ax.set_xlim(-3.65, n_cols - 0.35)
    ax.text(-3.55, n_rows + 0.16, "Case", ha="left", va="bottom", fontsize=9, fontweight="bold")
    ax.text(-1.35, n_rows + 0.16, "Method", ha="left", va="bottom", fontsize=9, fontweight="bold")
    fig.savefig(output_path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_funky_heatmap(metrics_path: Path, output_dir: Path, stem: str) -> tuple[Path, Path, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    heatmap_data, metric_info = prepare_heatmap_data(metrics_path)

    input_path = output_dir / f"{stem}_input.csv"
    serializable = heatmap_data.copy()
    serializable["criteria_pie"] = serializable["criteria_pie"].map(str)
    serializable.to_csv(input_path, index=False)

    output_path = output_dir / f"{stem}.png"
    funky_heatmap, engine = _load_funky_heatmap()
    if funky_heatmap is not None:
        try:
            _render_with_funky(funky_heatmap, heatmap_data, metric_info, output_path)
            return output_path, input_path, engine
        except Exception as exc:
            engine = f"matplotlib fallback after {engine} failed: {exc}"

    _render_matplotlib(heatmap_data, metric_info, output_path)
    return output_path, input_path, engine


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--stem", default="spatial_gradient_cpp_funky")
    args = parser.parse_args()

    output_path, input_path, engine = render_funky_heatmap(
        metrics_path=args.metrics,
        output_dir=args.output_dir,
        stem=args.stem,
    )
    print(f"engine={engine}")
    print(f"png={output_path}")
    print(f"input={input_path}")


if __name__ == "__main__":
    main()
