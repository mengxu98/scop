#!/usr/bin/env python
"""Render CIBERSORT backend benchmark metrics as a funky heatmap.

The script prefers ``omicverse.pl.funky_heatmap`` when OmicVerse is installed
and falls back to ``pyfunkyheatmap`` for lightweight local rendering.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ACCEPTANCE = {
    "pearson": 0.995,
    "mae": 0.01,
    "major": 0.95,
    "pvalue_trend": 0.90,
}


def _load_funky_heatmap():
    try:
        import omicverse as ov  # type: ignore

        return ov.pl.funky_heatmap, "omicverse.pl.funky_heatmap"
    except ImportError:
        from pyfunkyheatmap import funky_heatmap  # type: ignore

        return funky_heatmap, "pyfunkyheatmap.funky_heatmap"


def _format_time(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.01:
        return f"{value:.3f}"
    if value < 10:
        return f"{value:.2f}"
    return f"{value:.1f}"


def _format_speedup(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value >= 100:
        return f"{value:.0f}x"
    if value >= 10:
        return f"{value:.1f}x"
    return f"{value:.2f}x"


def _format_number(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def _clip01(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").clip(lower=0, upper=1)


def _score_error(series: pd.Series, threshold: float) -> pd.Series:
    """Convert lower-is-better error into a higher-is-better [0, 1] score."""
    values = pd.to_numeric(series, errors="coerce")
    return (1 - (values / threshold).clip(lower=0, upper=1)).fillna(0)


def prepare_metrics(metrics_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metrics = pd.read_csv(metrics_path)
    numeric_cols = [
        "time_r",
        "time_cpp",
        "speedup",
        "fraction_pearson",
        "max_abs_diff",
        "mae",
        "rmse",
        "mean_rank_spearman",
        "top3_jaccard",
        "major_cell_consistency",
        "pvalue_spearman",
        "pvalue_bin_consistency",
    ]
    for col in numeric_cols:
        if col not in metrics.columns:
            metrics[col] = np.nan
        metrics[col] = pd.to_numeric(metrics[col], errors="coerce")

    labels = {
        "lm22_20": "LM22 x20",
        "islet_bulk": "islet bulk",
        "small": "small toy",
    }
    metrics["dataset_label"] = metrics["dataset"].replace(labels)
    metrics["setting"] = (
        "perm="
        + metrics["perm"].astype(str)
        + ", QN="
        + metrics["QN"].astype(str).str.upper()
    )
    metrics["id"] = metrics["dataset_label"] + " | " + metrics["setting"]

    max_speedup = max(float(metrics["speedup"].max(skipna=True)), 1.0)
    metrics["speedup_score"] = np.log1p(metrics["speedup"]) / np.log1p(max_speedup)
    metrics["r_time_label"] = metrics["time_r"].map(_format_time)
    metrics["cpp_time_label"] = metrics["time_cpp"].map(_format_time)
    metrics["speedup_label"] = metrics["speedup"].map(_format_speedup)

    metrics["pearson_score"] = _clip01(metrics["fraction_pearson"])
    metrics["mae_score"] = _score_error(metrics["mae"], ACCEPTANCE["mae"])
    metrics["rmse_score"] = _score_error(metrics["rmse"], ACCEPTANCE["mae"])
    metrics["rank_score"] = _clip01(metrics["mean_rank_spearman"])
    metrics["top3_score"] = _clip01(metrics["top3_jaccard"])
    metrics["major_score"] = _clip01(metrics["major_cell_consistency"])
    metrics["pvalue_score"] = _clip01(metrics["pvalue_spearman"])
    metrics["pbin_score"] = _clip01(metrics["pvalue_bin_consistency"])

    metrics["pearson_label"] = metrics["fraction_pearson"].map(lambda x: _format_number(x, 4))
    metrics["mae_label"] = metrics["mae"].map(lambda x: _format_number(x, 2))
    metrics["rmse_label"] = metrics["rmse"].map(lambda x: _format_number(x, 2))
    metrics["rank_label"] = metrics["mean_rank_spearman"].map(lambda x: _format_number(x, 3))
    metrics["top3_label"] = metrics["top3_jaccard"].map(lambda x: _format_number(x, 3))
    metrics["major_label"] = metrics["major_cell_consistency"].map(lambda x: _format_number(x, 3))
    metrics["pvalue_label"] = metrics["pvalue_spearman"].map(
        lambda x: "NA" if pd.isna(x) else _format_number(x, 3)
    )
    metrics["pbin_label"] = metrics["pvalue_bin_consistency"].map(
        lambda x: "NA" if pd.isna(x) else _format_number(x, 3)
    )

    pass_pearson = metrics["fraction_pearson"] > ACCEPTANCE["pearson"]
    pass_mae = metrics["mae"] < ACCEPTANCE["mae"]
    pass_major = metrics["major_cell_consistency"] > ACCEPTANCE["major"]
    pass_pvalue = metrics["pvalue_spearman"].isna() | (
        metrics["pvalue_spearman"] > ACCEPTANCE["pvalue_trend"]
    )
    metrics["pass_count"] = (
        pass_pearson.astype(int)
        + pass_mae.astype(int)
        + pass_major.astype(int)
        + pass_pvalue.astype(int)
    )
    metrics["criteria_pie"] = metrics["pass_count"].map(
        lambda count: {"pass": int(count), "watch": max(4 - int(count), 1e-6)}
    )
    metrics["criteria_label"] = metrics["pass_count"].astype(str) + "/4"

    heatmap_data = metrics[
        [
            "id",
            "setting",
            "r_time_label",
            "cpp_time_label",
            "speedup_score",
            "speedup_label",
            "pearson_score",
            "pearson_label",
            "mae_score",
            "mae_label",
            "rmse_score",
            "rmse_label",
            "rank_score",
            "rank_label",
            "top3_score",
            "top3_label",
            "major_score",
            "major_label",
            "pvalue_score",
            "pvalue_label",
            "pbin_score",
            "pbin_label",
            "criteria_pie",
            "criteria_label",
        ]
    ].copy()
    return metrics, heatmap_data


def _metric_columns() -> pd.DataFrame:
    """Column specification with raw-value text overlaid on glyph columns."""
    cols = [
        {"id": "setting", "name": "Setting", "geom": "text", "group": "Design", "width": 5.6, "legend": False},
        {"id": "r_time_label", "name": "R sec", "geom": "text", "group": "Runtime", "width": 2.0, "legend": False},
        {"id": "cpp_time_label", "name": "C++ sec", "geom": "text", "group": "Runtime", "width": 2.1, "legend": False},
        {"id": "speedup_score", "name": "Speedup", "geom": "bar", "group": "Runtime", "palette": "speed", "width": 3.0, "legend": True},
        {"id": "speedup_label", "name": "", "geom": "text", "group": "Runtime", "overlay": True, "width": 0.0, "legend": False},
        {"id": "pearson_score", "name": "Pearson", "geom": "circle", "group": "Fraction", "palette": "quality", "width": 1.7, "legend": True},
        {"id": "pearson_label", "name": "", "geom": "text", "group": "Fraction", "overlay": True, "width": 0.0, "legend": False},
        {"id": "mae_score", "name": "MAE", "geom": "funkyrect", "group": "Fraction", "palette": "quality", "width": 1.8, "legend": True},
        {"id": "mae_label", "name": "", "geom": "text", "group": "Fraction", "overlay": True, "width": 0.0, "legend": False},
        {"id": "rmse_score", "name": "RMSE", "geom": "funkyrect", "group": "Fraction", "palette": "quality", "width": 1.8, "legend": True},
        {"id": "rmse_label", "name": "", "geom": "text", "group": "Fraction", "overlay": True, "width": 0.0, "legend": False},
        {"id": "rank_score", "name": "Rank rho", "geom": "circle", "group": "Ranking", "palette": "rank", "width": 1.8, "legend": True},
        {"id": "rank_label", "name": "", "geom": "text", "group": "Ranking", "overlay": True, "width": 0.0, "legend": False},
        {"id": "top3_score", "name": "Top3", "geom": "circle", "group": "Ranking", "palette": "rank", "width": 1.6, "legend": True},
        {"id": "top3_label", "name": "", "geom": "text", "group": "Ranking", "overlay": True, "width": 0.0, "legend": False},
        {"id": "major_score", "name": "Top1", "geom": "circle", "group": "Ranking", "palette": "rank", "width": 1.6, "legend": True},
        {"id": "major_label", "name": "", "geom": "text", "group": "Ranking", "overlay": True, "width": 0.0, "legend": False},
        {"id": "pvalue_score", "name": "P rho", "geom": "circle", "group": "Permutation", "palette": "pvalue", "width": 1.8, "legend": True},
        {"id": "pvalue_label", "name": "", "geom": "text", "group": "Permutation", "overlay": True, "width": 0.0, "legend": False},
        {"id": "pbin_score", "name": "P bins", "geom": "circle", "group": "Permutation", "palette": "pvalue", "width": 1.8, "legend": True},
        {"id": "pbin_label", "name": "", "geom": "text", "group": "Permutation", "overlay": True, "width": 0.0, "legend": False},
        {"id": "criteria_pie", "name": "Criteria", "geom": "pie", "group": "Criteria", "palette": "criteria", "width": 1.7, "legend": True},
        {"id": "criteria_label", "name": "", "geom": "text", "group": "Criteria", "overlay": True, "width": 0.0, "legend": False},
    ]
    return pd.DataFrame(cols)


def render_funky_heatmap(metrics_path: Path, output_dir: Path, stem: str) -> tuple[Path, Path, str]:
    funky_heatmap, engine = _load_funky_heatmap()
    metrics, heatmap_data = prepare_metrics(metrics_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_path = output_dir / f"{stem}_input.csv"
    heatmap_data.to_csv(input_path, index=False)

    column_info = _metric_columns()
    group_names = ["Design", "Runtime", "Fraction", "Ranking", "Permutation", "Criteria"]
    column_groups = pd.DataFrame({"group": group_names, "level1": group_names})

    dataset_labels = list(dict.fromkeys(metrics["dataset_label"]))
    row_info = pd.DataFrame({"id": heatmap_data["id"], "group": metrics["dataset_label"]})
    row_groups = pd.DataFrame({"group": dataset_labels, "level1": dataset_labels})
    palettes = {
        "speed": ["#f7fbff", "#c6dbef", "#6baed6", "#08519c"],
        "quality": ["#fff7bc", "#fec44f", "#fd8d3c", "#d95f0e"],
        "rank": ["#f7fcf5", "#bae4b3", "#74c476", "#238b45"],
        "pvalue": ["#f7fbff", "#bdd7e7", "#6baed6", "#2171b5"],
        "criteria": {"pass": "#1b9e77", "watch": "#d9d9d9"},
    }
    legends = [
        {"title": "Speedup", "palette": "speed", "geom": "bar"},
        {"title": "Fraction score", "palette": "quality", "geom": "circle"},
        {"title": "Ranking score", "palette": "rank", "geom": "circle"},
        {"title": "P-value score", "palette": "pvalue", "geom": "circle"},
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
            "expand_xmax": 2.5,
            "col_annot_offset": 3.6,
            "col_annot_angle": 35,
            "col_bigspace": 0.8,
            "cell_text_size": 3.4,
            "col_annot_size": 3.8,
        },
        scale_column=False,
        add_abc=False,
        fig_scale=0.26,
        dpi=180,
    )

    output_path = output_dir / f"{stem}.png"
    if hasattr(figure, "save"):
        figure.save(output_path, facecolor="white")
    elif hasattr(figure, "figure"):
        figure.figure.savefig(output_path, facecolor="white")
    else:
        figure.savefig(output_path, facecolor="white")
    return output_path, input_path, engine


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--stem", default="cibersort_cpp_funky")
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
