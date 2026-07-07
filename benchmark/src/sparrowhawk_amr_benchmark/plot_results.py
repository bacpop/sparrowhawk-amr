from __future__ import annotations

import argparse
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable

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
STACKED_PALETTE = [
    "#C5D5F0",
    "#4F0427",
    "#44F270",
    "#F612A8",
    "#83AA3E",
    "#760796",
    "#C9F07F",
    "#270FE2",
    "#E9C338",
    "#2E2274",
    "#87D6BC",
    "#AD2A58",
    "#0CA82E",
    "#CE3BF6",
    "#1B511D",
    "#DD8EEB",
    "#111111",
    "#FA8D80",
    "#1F9383",
    "#F8381B",
    "#3D84E3",
    "#9E4302",
    "#074D65",
    "#FD8F20",
    "#874E57",
    "#7A8683",
    "#957206",
]
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
CLASS_ALIASES = {
    "FUSICDIC ACID": "FUSIDIC_ACID",
    "FUSICDIC_ACID": "FUSIDIC_ACID",
    "FUSIDIC ACID": "FUSIDIC_ACID",
}
BIOCIDE_LABELS = {
    "BIOCIDE",
    "BIOCIDES",
    "QUATERNARY AMMONIUM",
}
METAL_LABELS = {
    "METAL",
    "METALS",
    "ARSENIC",
    "ARSENATE",
    "ARSENITE",
    "ORGANOARSENIC",
    "CADMIUM",
    "CADMIUM/COBALT/NICKEL",
    "CADMIUM/LEAD/ZINC",
    "CHROMATE",
    "COPPER",
    "COPPER/GOLD",
    "COPPER/NICKEL",
    "COPPER/SILVER",
    "FLUORIDE",
    "GOLD",
    "MERCURY",
    "ORGANOMERCURY",
    "PHENYLMERCURY",
    "NICKEL",
    "SILVER",
    "TELLURIUM",
}


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


def set_text_font(axis: plt.Axes, regular: FontProperties | None, include_y: bool = True) -> None:
    if not regular:
        return
    labels = axis.get_xticklabels()
    if include_y:
        labels += axis.get_yticklabels()
    for label in labels:
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


def style_axis(axis: plt.Axes, regular: FontProperties | None, minor: bool = True, include_y_font: bool = True) -> None:
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
    set_text_font(axis, regular, include_y=include_y_font)


def save_figure(fig: plt.Figure, out_dir: Path, basename: str, formats: list[str]) -> None:
    ensure_dir(out_dir)
    for fmt in formats:
        fig.savefig(out_dir / f"{basename}.{fmt}", bbox_inches="tight")
    plt.close(fig)


def write_metric_csv(path: Path, dataframe: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(path, index=False)


def read_plot_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False)


def reject_unresolved_labels(dataframe: pd.DataFrame, column: str, path: Path) -> None:
    if column not in dataframe.columns:
        return
    labels = dataframe[column].astype(str).str.strip().str.upper()
    bad = labels.isin({"", "NA", "UNCLASSIFIED"})
    if bad.any():
        examples = dataframe.loc[bad, column].head(5).astype(str).tolist()
        raise SystemExit(
            f"{path} contains unresolved {column} labels {examples}; "
            "rerun comparison after subtype fallback"
        )


def require_columns(dataframe: pd.DataFrame, path: Path, columns: set[str]) -> None:
    missing = columns - set(dataframe.columns)
    if missing:
        raise SystemExit(
            f"{path} is missing required columns {sorted(missing)}; "
            "rerun amr-compare-amrfinder-batch with the updated benchmark code"
        )


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def safe_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    denominator = denominator.replace(0, np.nan)
    return (numerator / denominator).fillna(0)


def is_missing(value: object) -> bool:
    return value is None or bool(pd.isna(value))


def split_classes(raw: object) -> list[str]:
    if is_missing(raw):
        return []
    return [part.strip() for part in str(raw).split(";") if part.strip()]


def normalise_class_name(value: object) -> str:
    if is_missing(value):
        return "Unclassified"
    text = str(value).strip()
    if not text:
        return "Unclassified"
    return CLASS_ALIASES.get(text.upper(), text)


def normalise_species_name(value: object) -> str:
    if is_missing(value):
        return "Unknown species"
    text = str(value).strip()
    return text if text else "Unknown species"


def aggregate_species_name(value: object) -> str:
    species = normalise_species_name(value)
    genus = species.split(" ", 1)[0]
    if genus == "Shigella":
        return "Shigella spp."
    if genus == "Enterobacter":
        return "Enterobacter spp."
    if genus == "Providencia":
        return "Providencia spp."
    return species


def aggregate_class_name(value: object) -> str:
    class_name = normalise_class_name(value)
    class_upper = class_name.upper()
    if class_upper in METAL_LABELS:
        return "Stress/Metal"
    if class_upper in BIOCIDE_LABELS:
        return "Stress/Biocide"
    if class_upper == "ACID":
        return "Stress/Acid"
    if class_upper == "HEAT":
        return "Stress/Heat"
    return class_name


def aggregate_concern_group(row: pd.Series) -> str:
    class_name = normalise_class_name(row.get("class_name", ""))
    type_name = str(row.get("type_name", "") or "").strip().upper()
    subtype = str(row.get("subtype", "") or "").strip().upper()
    class_upper = class_name.upper()

    if type_name == "VIRULENCE":
        return "Virulence"
    if type_name == "STRESS":
        if subtype == "METAL" or class_upper in METAL_LABELS:
            return "Stress/Metal"
        if subtype == "BIOCIDE" or class_upper in BIOCIDE_LABELS:
            return "Stress/Biocide"
        if subtype == "ACID" or class_upper == "ACID":
            return "Stress/Acid"
        if subtype == "HEAT" or class_upper == "HEAT":
            return "Stress/Heat"
        raise ValueError(f"Unexpected STRESS subtype for aggregation: {subtype!r}")
    if class_upper in METAL_LABELS:
        return "Stress/Metal"
    if class_upper in BIOCIDE_LABELS:
        return "Stress/Biocide"
    if class_upper == "ACID":
        return "Stress/Acid"
    if class_upper == "HEAT":
        return "Stress/Heat"
    return aggregate_class_name(class_name)


def normalise_selected(dataframe: pd.DataFrame) -> pd.DataFrame:
    out = dataframe.copy()
    if "species" in out.columns:
        out["species"] = out["species"].map(normalise_species_name)
    else:
        out["species"] = "Unknown species"
    classes = out.get("classes", pd.Series([""] * len(out), index=out.index))
    out["classes"] = classes.map(lambda raw: ";".join(normalise_class_name(item) for item in split_classes(raw)))
    return out


def normalise_metric_labels(dataframe: pd.DataFrame) -> pd.DataFrame:
    out = dataframe.copy()
    if "species" in out.columns:
        out["species"] = out["species"].map(normalise_species_name)
    if "class_name" in out.columns:
        out["class_name"] = out["class_name"].map(normalise_class_name)
    if "type_name" in out.columns:
        out["type_name"] = out["type_name"].fillna("Unclassified").replace("", "Unclassified")
    return out


def metric_count_columns(dataframe: pd.DataFrame) -> list[str]:
    return [column for column in dataframe.columns if column.endswith(("_tp", "_fp", "_fn", "_tn"))]


def recompute_metric_fields(dataframe: pd.DataFrame) -> pd.DataFrame:
    out = dataframe.copy()
    for prefix in ("exact", "report_unit"):
        required = {f"{prefix}_tp", f"{prefix}_fp", f"{prefix}_fn", f"{prefix}_tn"}
        if not required <= set(out.columns):
            continue
        precision = safe_ratio(out[f"{prefix}_tp"], out[f"{prefix}_tp"] + out[f"{prefix}_fp"])
        recall = safe_ratio(out[f"{prefix}_tp"], out[f"{prefix}_tp"] + out[f"{prefix}_fn"])
        out[f"{prefix}_precision"] = precision
        out[f"{prefix}_recall"] = recall
        out[f"{prefix}_sensitivity"] = recall
        out[f"{prefix}_specificity"] = safe_ratio(out[f"{prefix}_tn"], out[f"{prefix}_tn"] + out[f"{prefix}_fp"])
        out[f"{prefix}_f1"] = safe_ratio(2 * precision * recall, precision + recall)
    return out.fillna(0)


def merge_metric_rows(dataframe: pd.DataFrame, label_columns: str | list[str]) -> pd.DataFrame:
    labels = [label_columns] if isinstance(label_columns, str) else list(label_columns)
    if dataframe.empty or any(label not in dataframe.columns for label in labels):
        return dataframe
    count_columns = metric_count_columns(dataframe)
    if not count_columns:
        return dataframe
    group_columns = [column for column in CONFIG_COLUMNS if column in dataframe.columns] + labels
    sum_columns = count_columns + (["assemblies_compared"] if "assemblies_compared" in dataframe.columns else [])
    grouped = dataframe.groupby(group_columns, dropna=False)[sum_columns].sum().reset_index()
    return recompute_metric_fields(grouped)


def aggregate_metric_rows(dataframe: pd.DataFrame, label_column: str, mapper: Callable[[object], str]) -> pd.DataFrame:
    if dataframe.empty or label_column not in dataframe.columns:
        return dataframe
    out = dataframe.copy()
    out[label_column] = out[label_column].map(mapper)
    return merge_metric_rows(out, label_column)


def aggregate_metric_rows_by_row(dataframe: pd.DataFrame, label_column: str, mapper: Callable[[pd.Series], str]) -> pd.DataFrame:
    if dataframe.empty or label_column not in dataframe.columns:
        return dataframe
    out = dataframe.copy()
    out[label_column] = out.apply(mapper, axis=1)
    return merge_metric_rows(out, label_column)


def aggregate_species_class_metric_rows(dataframe: pd.DataFrame) -> pd.DataFrame:
    if dataframe.empty:
        return dataframe
    out = dataframe.copy()
    out["species"] = out["species"].map(aggregate_species_name)
    out["class_name"] = out.apply(aggregate_concern_group, axis=1)
    return merge_metric_rows(out, ["species", "class_name"])


def metric_label_counts(dataframe: pd.DataFrame, label_column: str) -> pd.Series:
    if dataframe.empty or label_column not in dataframe.columns or "assemblies_compared" not in dataframe.columns:
        return pd.Series(dtype=float)
    return dataframe.groupby(label_column, dropna=False)["assemblies_compared"].sum().sort_values(ascending=False)


def species_class_table_from_metrics(dataframe: pd.DataFrame) -> pd.DataFrame:
    if dataframe.empty or not {"species", "class_name", "assemblies_compared"} <= set(dataframe.columns):
        return pd.DataFrame()
    return dataframe.pivot_table(
        index="species",
        columns="class_name",
        values="assemblies_compared",
        aggfunc="sum",
        fill_value=0,
    )


def aggregate_selected(dataframe: pd.DataFrame) -> pd.DataFrame:
    out = dataframe.copy()
    out["species"] = out["species"].map(aggregate_species_name)
    out["classes"] = out["classes"].map(
        lambda raw: ";".join(sorted({aggregate_class_name(item) for item in split_classes(raw)}))
    )
    return out


def warn_unclassified_sources(selected: pd.DataFrame, class_metrics: pd.DataFrame, species_class_metrics: pd.DataFrame | None) -> None:
    selected_empty = int(selected.get("classes", pd.Series(dtype=str)).map(lambda raw: len(split_classes(raw)) == 0).sum())
    metric_unclassified = int((class_metrics.get("class_name", pd.Series(dtype=str)) == "Unclassified").sum())
    if selected_empty:
        print(f"Warning: {selected_empty} selected assemblies have no class metadata and will be plotted as Unclassified.")
    if metric_unclassified:
        print(
            "Warning: class_metrics.csv contains Unclassified rows. "
            "This usually means AMRFinderPlus Class or detector class_name was blank/missing for at least one call."
        )
    if species_class_metrics is not None and "class_name" in species_class_metrics.columns:
        species_class_unclassified = int((species_class_metrics["class_name"] == "Unclassified").sum())
        if species_class_unclassified:
            print(f"Warning: species_class_metrics.csv contains {species_class_unclassified} Unclassified rows.")


def ordered_labels_from_counts(counts: pd.Series, max_labels: int) -> list[str]:
    ordered = counts.sort_values(ascending=False)
    if max_labels > 0:
        ordered = ordered.head(max_labels)
    return [str(label) for label in ordered.index]


def append_metric_only_labels(order: list[str], dataframe: pd.DataFrame, label_column: str, max_labels: int) -> list[str]:
    seen = set(order)
    extra = []
    if label_column not in dataframe.columns:
        return order
    for label in dataframe[label_column].dropna().astype(str):
        if label not in seen:
            seen.add(label)
            extra.append(label)
    if max_labels > 0:
        available = max(0, max_labels - len(order))
        extra = extra[:available]
    return order + extra


def apply_order_to_metric_rows(dataframe: pd.DataFrame, label_column: str, order: list[str]) -> pd.DataFrame:
    if dataframe.empty:
        return dataframe
    rows = dataframe.set_index(label_column, drop=False)
    ordered_rows = [rows.loc[label] for label in order if label in rows.index]
    return pd.DataFrame(ordered_rows) if ordered_rows else dataframe.iloc[0:0]


def apply_order_to_counts(counts: pd.Series, order: list[str]) -> pd.Series:
    return counts.reindex([label for label in order if label in counts.index]).dropna()


def apply_order_to_table(table: pd.DataFrame, row_order: list[str], column_order: list[str]) -> pd.DataFrame:
    rows = [label for label in row_order if label in table.index]
    cols = [label for label in column_order if label in table.columns]
    return table.reindex(index=rows, columns=cols).fillna(0)


def set_axis_ticklabels(axis: plt.Axes, labels: list[str], regular: FontProperties | None, italic: FontProperties | None, italic_labels: bool) -> None:
    axis.set_yticklabels(labels, fontsize=LABEL_SIZE)
    if italic_labels and italic:
        for label in axis.get_yticklabels():
            label.set_fontproperties(italic)
    elif regular:
        for label in axis.get_yticklabels():
            label.set_fontproperties(regular)


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


def metric_plot(
    dataframe: pd.DataFrame,
    label_column: str,
    title: str,
    basename: str,
    out_dir: Path,
    formats: list[str],
    regular: FontProperties | None,
    italic: FontProperties | None,
    bold: FontProperties | None,
    preliminary: bool,
    counts: pd.Series | None = None,
    italic_labels: bool = False,
) -> None:
    required = {label_column, "report_unit_sensitivity", "report_unit_specificity", "report_unit_precision", "report_unit_f1"}
    missing = required - set(dataframe.columns)
    if missing:
        raise SystemExit(f"Metrics file is missing required columns: {', '.join(sorted(missing))}")
    plot_df = dataframe.copy()
    plot_df["report_unit_sensitivity"] = plot_df["report_unit_sensitivity"].map(safe_float)
    plot_df["report_unit_specificity"] = plot_df["report_unit_specificity"].map(safe_float)
    plot_df["report_unit_precision"] = plot_df["report_unit_precision"].map(safe_float)
    plot_df["report_unit_f1"] = plot_df["report_unit_f1"].map(safe_float)
    labels = plot_df[label_column].astype(str).tolist()
    display_labels = []
    for label, (_, row) in zip(labels, plot_df.iterrows()):
        if counts is not None and label in counts.index:
            display_labels.append(f"{label} (n={int(counts[label])})")
        elif "assemblies_compared" in row.index:
            display_labels.append(f"{label} (n={int(row['assemblies_compared'])})")
        else:
            display_labels.append(label)
    y = np.arange(len(labels))
    height = max(4.0, 0.31 * len(labels) + 1.6)
    fig, ax = plt.subplots(figsize=(8.5, height), dpi=150)
    ax.barh(y, plot_df["report_unit_specificity"], height=0.72, color=RED, alpha=0.45, label="Specificity", zorder=1)
    ax.barh(y, plot_df["report_unit_sensitivity"], height=0.42, color=BLUE, alpha=0.80, label="Sensitivity", zorder=2)
    ax.scatter(plot_df["report_unit_precision"], y, color=GREY, marker="D", s=18, label="Precision", zorder=4)
    ax.scatter(plot_df["report_unit_f1"], y, color=BLACK, s=18, label="F1", zorder=5)
    ax.set_yticks(y)
    set_axis_ticklabels(ax, display_labels, regular, italic, italic_labels)
    ax.set_xlim(0.0, 1.0)
    if len(labels):
        ax.set_ylim(-0.5, len(labels) - 0.5)
    ax.margins(y=0)
    ax.set_xlabel("Value", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    ax.legend(
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="white",
        framealpha=0.92,
        prop=regular,
        handlelength=0.7,
        handletextpad=0.6,
        labelspacing=0.3,
    )
    add_headers(ax, title, "Report-unit metrics", regular, italic, bold, preliminary)
    style_axis(ax, regular, include_y_font=not italic_labels)
    save_figure(fig, out_dir, basename, formats)


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


def count_plot(counts: pd.Series, title: str, basename: str, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, italic_labels: bool = False) -> None:
    if counts.empty:
        raise SystemExit(f"No values available for {title}")
    counts = counts.iloc[::-1]
    labels = [str(label) for label in counts.index]
    y = np.arange(len(labels))
    height = max(4.0, 0.31 * len(labels) + 1.6)
    fig, ax = plt.subplots(figsize=(8.5, height), dpi=150)
    ax.barh(y, counts.values, height=0.62, color=GREEN, alpha=0.85)
    ax.set_yticks(y)
    set_axis_ticklabels(ax, labels, regular, italic, italic_labels)
    ax.set_xlim(0, max(float(counts.max()) * 1.08, 1.0))
    ax.set_xlabel("Assemblies", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    add_headers(ax, title, "Assembly counts", regular, italic, bold, preliminary)
    style_axis(ax, regular, include_y_font=not italic_labels)
    save_figure(fig, out_dir, basename, formats)


def species_class_table(selected: pd.DataFrame) -> pd.DataFrame:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for _, row in selected.iterrows():
        species = normalise_species_name(row.get("species", ""))
        classes = split_classes(row.get("classes", "")) or ["Unclassified"]
        for class_name in classes:
            counts[species][normalise_class_name(class_name)] += 1
    return pd.DataFrame(counts).T.fillna(0).astype(float)


def stacked_plot(table: pd.DataFrame, title: str, contribution_label: str, basename: str, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, italic_labels: bool = False, italic_contributions: bool = False) -> None:
    if table.empty:
        raise SystemExit(f"No values available for {title}")
    table = table.iloc[::-1]
    labels = [str(label) for label in table.index]
    contributions = list(table.columns)
    y = np.arange(len(labels))
    height = max(4.5, 0.33 * len(labels) + 1.8)
    fig, ax = plt.subplots(figsize=(9.5, height), dpi=150)
    left = np.zeros(len(labels))
    warned_palette_wrap = False
    for idx, contribution in enumerate(contributions):
        values = table[contribution].to_numpy(dtype=float)
        if np.all(values == 0):
            continue
        if idx >= len(STACKED_PALETTE) and not warned_palette_wrap:
            print(
                f"Warning: stacked plot palette has {len(STACKED_PALETTE)} colours; "
                "additional contributions will reuse colours."
            )
            warned_palette_wrap = True
        colour = STACKED_PALETTE[idx % len(STACKED_PALETTE)]
        ax.barh(y, values, left=left, height=0.68, color=colour, label=str(contribution), alpha=0.88)
        left += values
    ax.set_yticks(y)
    set_axis_ticklabels(ax, labels, regular, italic, italic_labels)
    ax.set_xscale("log")
    positive_total = table.sum(axis=1)
    positive_total = positive_total[positive_total > 0]
    ax.set_xlim(0.8, max(float(positive_total.max()) * 1.5, 1.2))
    ax.set_xlabel("Assemblies", fontproperties=regular, loc="right", fontsize=AXIS_TITLE_SIZE)
    ax.set_ylabel(title, fontproperties=regular, loc="top", fontsize=AXIS_TITLE_SIZE)
    add_headers(ax, title, f"Stacked by {contribution_label}", regular, italic, bold, preliminary)
    style_axis(ax, regular, minor=False, include_y_font=not italic_labels)
    legend = ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, prop=regular, fontsize=8, handlelength=0.7, handletextpad=0.5, labelspacing=0.25)
    if italic_contributions and italic:
        for text in legend.get_texts():
            text.set_fontproperties(italic)
    save_figure(fig, out_dir, basename, formats)


def plot_suite(
    selected: pd.DataFrame,
    species_metrics: pd.DataFrame,
    class_metrics: pd.DataFrame,
    out_dir: Path,
    formats: list[str],
    regular: FontProperties | None,
    italic: FontProperties | None,
    bold: FontProperties | None,
    preliminary: bool,
    max_labels: int,
    suffix: str = "",
    class_count_values: pd.Series | None = None,
    species_class_count_table: pd.DataFrame | None = None,
) -> None:
    species_count_values = species_counts(selected)
    if class_count_values is None:
        class_count_values = class_counts(selected)
    species_order = ordered_labels_from_counts(species_count_values, max_labels)
    class_order = ordered_labels_from_counts(class_count_values, max_labels)
    species_order = append_metric_only_labels(species_order, species_metrics, "species", max_labels)
    class_order = append_metric_only_labels(class_order, class_metrics, "class_name", max_labels)
    species_table = apply_order_to_metric_rows(species_metrics, "species", species_order)
    class_table = apply_order_to_metric_rows(class_metrics, "class_name", class_order)
    metric_plot(species_table.iloc[::-1], "species", "Species", f"report_unit_species_metrics{suffix}", out_dir, formats, regular, italic, bold, preliminary, species_count_values, italic_labels=True)
    metric_plot(class_table.iloc[::-1], "class_name", "Antibiotic class", f"report_unit_class_metrics{suffix}", out_dir, formats, regular, italic, bold, preliminary, class_count_values)
    count_plot(apply_order_to_counts(species_count_values, species_order), "Species", f"species_assembly_counts{suffix}", out_dir, formats, regular, italic, bold, preliminary, italic_labels=True)
    count_plot(apply_order_to_counts(class_count_values, class_order), "Antibiotic class", f"class_assembly_counts{suffix}", out_dir, formats, regular, italic, bold, preliminary)
    table = species_class_count_table if species_class_count_table is not None else species_class_table(selected)
    stacked_plot(apply_order_to_table(table, species_order, class_order), "Species", "antibiotic class", f"species_counts_by_class_stacked{suffix}", out_dir, formats, regular, italic, bold, preliminary, italic_labels=True)
    stacked_plot(apply_order_to_table(table.T, class_order, species_order), "Antibiotic class", "species", f"class_counts_by_species_stacked{suffix}", out_dir, formats, regular, italic, bold, preliminary, italic_contributions=True)


def type_metric_plot(type_metrics: pd.DataFrame, out_dir: Path, formats: list[str], regular: FontProperties | None, italic: FontProperties | None, bold: FontProperties | None, preliminary: bool, max_labels: int) -> None:
    if type_metrics.empty:
        return
    ordered = type_metrics.sort_values(["report_unit_f1", "type_name"], ascending=[False, True])
    if max_labels > 0:
        ordered = ordered.head(max_labels)
    metric_plot(
        ordered.iloc[::-1],
        "type_name",
        "AMRFinderPlus type",
        "report_unit_type_metrics",
        out_dir,
        formats,
        regular,
        italic,
        bold,
        preliminary,
    )


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
    parser.add_argument("--type-metrics", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--mode")
    parser.add_argument("--k")
    parser.add_argument("--out-formats", default=",".join(DEFAULT_FORMATS))
    parser.add_argument("--max-labels", type=int, default=0)
    parser.add_argument("--no-preliminary", action="store_true")
    args = parser.parse_args()

    formats = parse_formats(args.out_formats)
    regular, italic, bold = load_fonts()
    aggregate = read_plot_csv(args.aggregate_metrics) if args.aggregate_metrics and args.aggregate_metrics.exists() else None
    config = choose_config(aggregate, args.mode, args.k)
    selected = normalise_selected(read_plot_csv(args.selected_manifest))
    species_metrics = merge_metric_rows(normalise_metric_labels(filter_config(read_plot_csv(args.species_metrics), config)), "species")
    class_raw = read_plot_csv(args.class_metrics)
    require_columns(class_raw, args.class_metrics, {"type_name", "subtype", "class_name"})
    reject_unresolved_labels(class_raw, "class_name", args.class_metrics)
    class_filtered = normalise_metric_labels(filter_config(class_raw, config))
    class_metrics_for_aggregation = class_filtered
    class_metrics = merge_metric_rows(class_filtered, "class_name")
    species_class_metrics = None
    if args.species_class_metrics and args.species_class_metrics.exists():
        species_class_raw = read_plot_csv(args.species_class_metrics)
        require_columns(species_class_raw, args.species_class_metrics, {"species", "type_name", "subtype", "class_name"})
        reject_unresolved_labels(species_class_raw, "class_name", args.species_class_metrics)
        species_class_metrics = normalise_metric_labels(filter_config(species_class_raw, config))
    type_metrics = None
    if args.type_metrics and args.type_metrics.exists():
        type_metrics = normalise_metric_labels(filter_config(read_plot_csv(args.type_metrics), config))
    preliminary = not args.no_preliminary

    warn_unclassified_sources(selected, class_metrics, species_class_metrics)
    ensure_dir(args.out_dir)
    plot_suite(selected, species_metrics, class_metrics, args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)
    if type_metrics is not None:
        type_metrics_for_plot = type_metrics
        if aggregate is not None and not aggregate.empty:
            total_row = filter_config(aggregate, config).iloc[[0]].copy()
            total_row["type_name"] = "Total"
            type_metrics_for_plot = pd.concat([type_metrics, total_row], ignore_index=True, sort=False)
        type_metric_plot(type_metrics_for_plot, args.out_dir, formats, regular, italic, bold, preliminary, args.max_labels)

    aggregated_selected = aggregate_selected(selected)
    aggregated_species_metrics = aggregate_metric_rows(species_metrics, "species", aggregate_species_name)
    aggregated_class_metrics = aggregate_metric_rows_by_row(
        class_metrics_for_aggregation,
        "class_name",
        aggregate_concern_group,
    )
    aggregated_species_class_metrics = None
    if species_class_metrics is not None:
        aggregated_species_class_metrics = aggregate_species_class_metric_rows(species_class_metrics)
        write_metric_csv(args.out_dir.parent / "species_class_metrics_aggregated.csv", aggregated_species_class_metrics)

    write_metric_csv(args.out_dir.parent / "species_metrics_aggregated.csv", aggregated_species_metrics)
    write_metric_csv(args.out_dir.parent / "class_metrics_aggregated.csv", aggregated_class_metrics)

    aggregated_class_counts = metric_label_counts(aggregated_class_metrics, "class_name")
    aggregated_species_class_table = (
        species_class_table_from_metrics(aggregated_species_class_metrics)
        if aggregated_species_class_metrics is not None
        else None
    )
    plot_suite(
        aggregated_selected,
        aggregated_species_metrics,
        aggregated_class_metrics,
        args.out_dir,
        formats,
        regular,
        italic,
        bold,
        preliminary,
        args.max_labels,
        "_aggregated",
        class_count_values=aggregated_class_counts,
        species_class_count_table=aggregated_species_class_table,
    )


if __name__ == "__main__":
    main()
