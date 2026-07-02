from __future__ import annotations

import argparse
import os
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import AutoMinorLocator

try:
    from .common import ensure_dir
except ImportError:
    from common import ensure_dir

RED = "#D41645"
BLUE = "#193F90"
GREEN = "#18974C"
BLACK = "#111111"
GREY = "#666666"
AXIS_TITLE_SIZE = 12
LABEL_SIZE = 9
TITLE_SIZE = 11
DEFAULT_FORMATS = ("png", "pdf", "svg")
CONFIG_COLUMNS = (
    "mode",
    "k",
    "threshold_mode",
    "min_gene_threshold",
    "min_report_unit_threshold",
)


def font_prop(env_name: str) -> FontProperties | None:
    path = os.environ.get(env_name, "")
    if path and Path(path).is_file():
        return FontProperties(fname=path)
    return None


def load_fonts() -> tuple[FontProperties | None, FontProperties | None, FontProperties | None]:
    regular = font_prop("IBM_PLEX_SANS_REGULAR")
    italic = font_prop("IBM_PLEX_SANS_ITALIC")
    bold = font_prop("IBM_PLEX_SANS_BOLD")
    if not regular:
        print("Warning: IBM Plex Sans regular font not configured; using matplotlib default.")
    return regular, italic or regular, bold or regular


def set_text_font(axis: plt.Axes, regular: FontProperties | None) -> None:
    if not regular:
        return
    for label in axis.get_xticklabels() + axis.get_yticklabels():
        label.set_fontproperties(regular)


def add_headers(axis: plt.Axes, left: str, right: str, regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool) -> None:
    axis.text(
        0,
        1.01,
        left,
        fontproperties=italic or regular,
        horizontalalignment="left",
        verticalalignment="bottom",
        transform=axis.transAxes,
    )
    axis.text(
        1,
        1.01,
        right,
        fontproperties=regular,
        horizontalalignment="right",
        verticalalignment="bottom",
        transform=axis.transAxes,
    )
    if preliminary:
        axis.text(
            0.5,
            1.01,
            "Preliminary",
            fontproperties=bold or regular,
            horizontalalignment="center",
            verticalalignment="bottom",
            transform=axis.transAxes,
        )


def style_axis(axis: plt.Axes, regular: FontProperties | None, minor: bool = True) -> None:
    axis.tick_params(which="major", direction="in")
    axis.tick_params(which="minor", direction="in")
    axis.xaxis.set_ticks_position("both")
    axis.yaxis.set_ticks_position("both")
    if minor:
        axis.xaxis.set_minor_locator(AutoMinorLocator())
        axis.yaxis.set_minor_locator(AutoMinorLocator())
    if axis.get_xscale() == "linear":
        axis.ticklabel_format(useMathText=True, axis="x", scilimits=(-4, 4))
    if regular:
        axis.get_xaxis().get_offset_text().set_fontproperties(regular)
    set_text_font(axis, regular)


def save_figure(fig: plt.Figure, out_dir: Path, basename: str, formats: list[str]) -> None:
    ensure_dir(out_dir)
    for fmt in formats:
        fig.savefig(out_dir / f"{basename}.{fmt}", bbox_inches="tight")
    plt.close(fig)


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def threshold_column(dataframe: pd.DataFrame) -> str:
    if "min_report_unit_threshold" in dataframe.columns:
        return "min_report_unit_threshold"
    if "min_gene_group_threshold" in dataframe.columns:
        return "min_gene_group_threshold"
    raise SystemExit("Metrics file is missing min_report_unit_threshold/min_gene_group_threshold")


def normalise_config_columns(dataframe: pd.DataFrame) -> pd.DataFrame:
    out = dataframe.copy()
    if "min_report_unit_threshold" not in out.columns and "min_gene_group_threshold" in out.columns:
        out["min_report_unit_threshold"] = out["min_gene_group_threshold"]
    return out


def choose_config(aggregate: pd.DataFrame | None, mode: str | None, k: str | None) -> dict[str, str]:
    if aggregate is None or aggregate.empty:
        config = {}
        if mode:
            config["mode"] = mode
        if k:
            config["k"] = str(k)
        return config

    aggregate = normalise_config_columns(aggregate)
    candidates = aggregate.copy()
    if mode:
        candidates = candidates[candidates["mode"].astype(str) == str(mode)]
    if k:
        candidates = candidates[candidates["k"].astype(str) == str(k)]
    if candidates.empty:
        raise SystemExit("No aggregate metrics row matches the requested mode/k")
    idx = candidates["report_unit_f1"].astype(float).idxmax()
    row = candidates.loc[idx]
    return {column: str(row[column]) for column in CONFIG_COLUMNS if column in row.index}


def filter_config(dataframe: pd.DataFrame, config: dict[str, str]) -> pd.DataFrame:
    dataframe = normalise_config_columns(dataframe)
    out = dataframe.copy()
    for column, value in config.items():
        if column not in out.columns:
            continue
        out = out[out[column].astype(str) == str(value)]
    if out.empty:
        config_text = ", ".join(f"{key}={value}" for key, value in config.items()) or "first available config"
        raise SystemExit(f"No metrics rows match {config_text}")
    return out


def maybe_limit(dataframe: pd.DataFrame, label_column: str, value_column: str, max_labels: int) -> pd.DataFrame:
    out = dataframe.sort_values([value_column, label_column], ascending=[False, True])
    if max_labels > 0:
        out = out.head(max_labels)
    return out.iloc[::-1]


def metric_plot(dataframe: pd.DataFrame, label_column: str, title: str, basename: str, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, max_labels: int) -> None:
    required = {label_column, "report_unit_sensitivity", "report_unit_specificity", "report_unit_f1"}
    missing = required - set(dataframe.columns)
    if missing:
        raise SystemExit(f"Metrics file is missing required columns: {', '.join(sorted(missing))}")
    plot_df = dataframe.copy()
    plot_df["report_unit_sensitivity"] = plot_df["report_unit_sensitivity"].map(safe_float)
    plot_df["report_unit_specificity"] = plot_df["report_unit_specificity"].map(safe_float)
    plot_df["report_unit_f1"] = plot_df["report_unit_f1"].map(safe_float)
    plot_df = maybe_limit(plot_df, label_column, "report_unit_f1", max_labels)
    labels = plot_df[label_column].astype(str).tolist()
    y = np.arange(len(labels))
    height = max(4.0, 0.31 * len(labels) + 1.6)
    fig, ax = plt.subplots(figsize=(8.5, height), dpi=150)
    ax.barh(y, plot_df["report_unit_specificity"], height=0.72, color=RED, alpha=0.45, label="Specificity", zorder=1)
    ax.barh(y, plot_df["report_unit_sensitivity"], height=0.42, color=BLUE, alpha=0.80, label="Sensitivity", zorder=2)
    ax.scatter(plot_df["report_unit_f1"], y, color=BLACK, s=18, label="F1", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=LABEL_SIZE)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Value", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    ax.legend(loc="lower right", frameon=False, prop=regular, handlelength=0.7, handletextpad=0.6, labelspacing=0.3)
    add_headers(ax, title, "Report-unit metrics", regular, italic, bold, preliminary)
    style_axis(ax, regular)
    save_figure(fig, out_dir, basename, formats)


def split_classes(raw: object) -> list[str]:
    return [part.strip() for part in str(raw or "").split(";") if part.strip()]


def species_counts(selected: pd.DataFrame) -> pd.Series:
    return selected["species"].fillna("Unknown species").replace("", "Unknown species").value_counts()


def class_counts(selected: pd.DataFrame) -> pd.Series:
    counter: Counter[str] = Counter()
    for raw in selected.get("classes", pd.Series(dtype=str)):
        classes = split_classes(raw)
        if not classes:
            counter["Unclassified"] += 1
        for class_name in classes:
            counter[class_name] += 1
    return pd.Series(counter).sort_values(ascending=False)


def count_plot(counts: pd.Series, title: str, basename: str, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, max_labels: int) -> None:
    if counts.empty:
        raise SystemExit(f"No values available for {title}")
    counts = counts.sort_values(ascending=False)
    if max_labels > 0:
        counts = counts.head(max_labels)
    counts = counts.iloc[::-1]
    labels = [str(label) for label in counts.index]
    y = np.arange(len(labels))
    height = max(4.0, 0.31 * len(labels) + 1.6)
    fig, ax = plt.subplots(figsize=(8.5, height), dpi=150)
    ax.barh(y, counts.values, height=0.62, color=GREEN, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=LABEL_SIZE)
    ax.set_xlim(0, max(float(counts.max()) * 1.08, 1.0))
    ax.set_xlabel("Assemblies", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    add_headers(ax, title, "Assembly counts", regular, italic, bold, preliminary)
    style_axis(ax, regular)
    save_figure(fig, out_dir, basename, formats)


def species_class_table(selected: pd.DataFrame) -> pd.DataFrame:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for _, row in selected.iterrows():
        species = str(row.get("species", "") or "Unknown species")
        classes = split_classes(row.get("classes", "")) or ["Unclassified"]
        for class_name in classes:
            counts[species][class_name] += 1
    return pd.DataFrame(counts).T.fillna(0).astype(float)


def maybe_limit_table(table: pd.DataFrame, max_labels: int) -> pd.DataFrame:
    if max_labels <= 0 or len(table.index) <= max_labels:
        return table
    return table.loc[table.sum(axis=1).sort_values(ascending=False).head(max_labels).index]


def stacked_plot(table: pd.DataFrame, title: str, contribution_label: str, basename: str, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, max_labels: int) -> None:
    table = maybe_limit_table(table, max_labels)
    table = table.loc[table.sum(axis=1).sort_values(ascending=True).index]
    if table.empty:
        raise SystemExit(f"No values available for {title}")
    labels = [str(label) for label in table.index]
    contributions = list(table.columns)
    y = np.arange(len(labels))
    height = max(4.5, 0.33 * len(labels) + 1.8)
    fig, ax = plt.subplots(figsize=(9.5, height), dpi=150)
    palette = plt.get_cmap("tab20")
    left = np.zeros(len(labels))
    for idx, contribution in enumerate(contributions):
        values = table[contribution].to_numpy(dtype=float)
        if np.all(values == 0):
            continue
        ax.barh(y, values, left=left, height=0.68, color=palette(idx % 20), label=str(contribution), alpha=0.88)
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=LABEL_SIZE)
    ax.set_xscale("log")
    positive_total = table.sum(axis=1)
    positive_total = positive_total[positive_total > 0]
    ax.set_xlim(0.8, max(float(positive_total.max()) * 1.5, 1.2))
    ax.set_xlabel("Assemblies", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    add_headers(ax, title, f"Stacked by {contribution_label}", regular, italic, bold, preliminary)
    style_axis(ax, regular, minor=False)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, prop=regular, fontsize=8, handlelength=0.7, handletextpad=0.5, labelspacing=0.25)
    save_figure(fig, out_dir, basename, formats)


def parse_formats(raw: str) -> list[str]:
    formats = [part.strip().lower().lstrip(".") for part in raw.split(",") if part.strip()]
    return formats or list(DEFAULT_FORMATS)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create plots from AMR benchmark outputs")
    parser.add_argument("--selected-manifest", type=Path, required=True)
    parser.add_argument("--aggregate-metrics", type=Path)
    parser.add_argument("--species-metrics", type=Path, required=True)
    parser.add_argument("--class-metrics", type=Path, required=True)
    parser.add_argument("--species-class-metrics", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--mode")
    parser.add_argument("--k")
    parser.add_argument("--out-formats", default=",".join(DEFAULT_FORMATS))
    parser.add_argument("--max-labels", type=int, default=0)
    parser.add_argument("--no-preliminary", action="store_true")
    args = parser.parse_args()

    formats = parse_formats(args.out_formats)
    regular, italic, bold = load_fonts()
    aggregate = pd.read_csv(args.aggregate_metrics) if args.aggregate_metrics and args.aggregate_metrics.exists() else None
    config = choose_config(aggregate, args.mode, args.k)
    selected = pd.read_csv(args.selected_manifest)
    species_metrics = filter_config(pd.read_csv(args.species_metrics), config)
    class_metrics = filter_config(pd.read_csv(args.class_metrics), config)
    preliminary = not args.no_preliminary

    ensure_dir(args.out_dir)
    metric_plot(species_metrics, "species", "Species", "report_unit_species_metrics", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    metric_plot(class_metrics, "class_name", "Antibiotic class", "report_unit_class_metrics", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    count_plot(species_counts(selected), "Species", "species_assembly_counts", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    count_plot(class_counts(selected), "Antibiotic class", "class_assembly_counts", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    table = species_class_table(selected)
    stacked_plot(table, "Species", "antibiotic class", "species_counts_by_class_stacked", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    stacked_plot(table.T, "Antibiotic class", "species", "class_counts_by_species_stacked", args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)


if __name__ == "__main__":
    main()
