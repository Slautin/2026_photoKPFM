from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = ROOT / "data"
DEFAULT_OUTPUT = ROOT / "reproduced_results" / "publication_figures"

COLORS = {
    "blue": "#3B6FA5",
    "orange": "#E6832A",
    "teal": "#2A9D8F",
    "gold": "#D9A62E",
    "red": "#D95D45",
    "gray": "#A9ADB3",
    "dark": "#252A31",
    "muted": "#56616F",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate the experimental PL/photoKPFM publication figures."
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.4,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )

def save_figure(figure: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(figure)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

def ternary_xy(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    fa = pd.to_numeric(frame["FAPbI3_pct"], errors="coerce").to_numpy(float) / 100.0
    pea = pd.to_numeric(frame["PEA2PbI4_pct"], errors="coerce").to_numpy(float) / 100.0
    return pea + 0.5 * fa, (math.sqrt(3.0) / 2.0) * fa

def draw_ternary(axis: plt.Axes) -> None:
    height = math.sqrt(3.0) / 2.0
    axis.plot([0, 1, 0.5, 0], [0, 0, height, 0], color=COLORS["dark"], lw=1.0)
    for fraction in (0.2, 0.4, 0.6, 0.8):
        axis.plot(
            [0.5 * fraction, 1 - 0.5 * fraction],
            [height * fraction] * 2,
            color="#E4E7EB",
            lw=0.5,
            zorder=0,
        )
        axis.plot(
            [fraction, 0.5 + 0.5 * fraction],
            [0, height * (1 - fraction)],
            color="#E4E7EB",
            lw=0.5,
            zorder=0,
        )
        axis.plot(
            [1 - fraction, 0.5 * (1 - fraction)],
            [0, height * (1 - fraction)],
            color="#E4E7EB",
            lw=0.5,
            zorder=0,
        )
    axis.text(-0.015, -0.055, r"BDAPbI$_4$", ha="center", va="top", fontsize=8.7)
    axis.text(1.015, -0.055, r"PEA$_2$PbI$_4$", ha="center", va="top", fontsize=8.7)
    axis.text(0.5, height + 0.035, r"FAPbI$_3$", ha="center", va="bottom", fontsize=8.7)
    axis.set_aspect("equal")
    axis.set_xlim(-0.075, 1.075)
    axis.set_ylim(-0.085, height + 0.075)
    axis.axis("off")

def selected_well_order(selection: pd.DataFrame) -> list[str]:
    data = selection.copy()
    data["retention"] = pd.to_numeric(data["retention"], errors="coerce")
    data["relative_range"] = pd.to_numeric(data["relative_range"], errors="coerce")
    data["peak_time_min"] = pd.to_numeric(data["peak_time_min"], errors="coerce")
    chosen: list[str] = []

    def add(value: object) -> None:
        if isinstance(value, str) and value and value not in chosen:
            chosen.append(value)

    stable = data.sort_values(["retention", "relative_range"], ascending=[False, True])
    if not stable.empty:
        add(stable.iloc[0]["well_id"])
    lower = data[data["temporal_behavior_class"].eq("Lower observed temporal variation")]
    if not lower.empty:
        add(lower.sort_values("relative_range").iloc[0]["well_id"])
    evolved = data[data["temporal_behavior_class"].eq("Pronounced temporal evolution")]
    if not evolved.empty:
        add(evolved.sort_values("relative_range", ascending=False).iloc[0]["well_id"])
    rise_decline = data[data["peak_time_min"].between(18, 171)]
    if not rise_decline.empty:
        add(rise_decline.sort_values("retention").iloc[0]["well_id"])
    if not stable.empty:
        add(stable.sort_values("retention").iloc[0]["well_id"])
    for well in data.sort_values("relative_range", ascending=False)["well_id"]:
        if len(chosen) >= 5:
            break
        add(well)
    return chosen[:5]

def annotate_selected_wells(
    axis: plt.Axes,
    frame: pd.DataFrame,
    selected: list[str],
) -> None:
    subset = frame[frame["well_id"].isin(selected)].copy()
    order = {well: index for index, well in enumerate(selected, start=1)}
    subset["selection_number"] = subset["well_id"].map(order)
    subset = subset.sort_values("selection_number")
    x, y = ternary_xy(subset)
    axis.scatter(
        x,
        y,
        s=72,
        facecolors="none",
        edgecolors=COLORS["dark"],
        linewidths=0.9,
        zorder=5,
    )
    offsets = {
        "G1": (0.018, 0.020),
        "B11": (0.022, 0.020),
        "C12": (0.026, 0.003),
        "H11": (-0.027, 0.015),
        "D5": (0.020, 0.020),
    }
    for row, xpos, ypos in zip(subset.itertuples(index=False), x, y):
        dx, dy = offsets.get(str(row.well_id), (0.02, 0.02))
        axis.text(
            xpos + dx,
            ypos + dy,
            str(int(row.selection_number)),
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="bold",
            color=COLORS["dark"],
            zorder=6,
        )

def plot_full_library_ternary(
    data: pd.DataFrame,
    matched_wells: set[str],
    selected: list[str],
    value_column: str,
    colorbar_label: str,
    stem: Path,
    cmap: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    figure, axis = plt.subplots(figsize=(4.45, 3.75))
    draw_ternary(axis)
    x, y = ternary_xy(data)
    values = pd.to_numeric(data[value_column], errors="coerce").to_numpy(float)
    scatter = axis.scatter(
        x,
        y,
        c=values,
        cmap=cmap,
        vmin=np.nanmin(values) if vmin is None else vmin,
        vmax=np.nanmax(values) if vmax is None else vmax,
        s=31,
        edgecolors="none",
        zorder=2,
    )
    matched = data["well_id"].isin(matched_wells).to_numpy()
    axis.scatter(
        x[matched],
        y[matched],
        s=50,
        facecolors="none",
        edgecolors=COLORS["dark"],
        linewidths=0.75,
        zorder=4,
    )
    annotate_selected_wells(axis, data, selected)
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.043, pad=0.02, shrink=0.83)
    colorbar.set_label(colorbar_label, fontsize=8.5)
    colorbar.ax.tick_params(labelsize=7.6, length=2.5)
    legend = Line2D(
        [],
        [],
        marker="o",
        linestyle="",
        markerfacecolor="none",
        markeredgecolor=COLORS["dark"],
        markeredgewidth=0.8,
        markersize=5.6,
        label="photoKPFM-matched subset",
    )
    axis.legend(
        handles=[legend],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.045),
        frameon=False,
        fontsize=7.7,
        borderaxespad=0,
    )
    figure.subplots_adjust(left=0.02, right=0.91, bottom=0.13, top=0.98)
    save_figure(figure, stem)

def validate_trace_endpoints(trace_data: pd.DataFrame) -> float:
    endpoints = (
        trace_data.sort_values("elapsed_time_min")
        .groupby("well_id", as_index=False)
        .tail(1)
    )
    error = np.abs(
        pd.to_numeric(endpoints["normalized_PL_peak_intensity"], errors="coerce")
        - pd.to_numeric(endpoints["retention"], errors="coerce")
    )
    return float(np.nanmax(error))

def composition_label(row: pd.Series | object) -> str:
    return (
        f"{int(float(getattr(row, 'FAPbI3_pct')))} / "
        f"{int(float(getattr(row, 'BDAPbI4_pct')))} / "
        f"{int(float(getattr(row, 'PEA2PbI4_pct')))}"
    )

def plot_trace_panel(
    data: pd.DataFrame,
    well_order: list[str],
    stem: Path,
    full_set: bool,
) -> None:
    if full_set:
        figure, axis = plt.subplots(figsize=(7.1, 4.0))
        palette = plt.get_cmap("tab20")(np.linspace(0, 0.95, len(well_order)))
        line_styles = ["-", "--", "-."]
    else:
        figure, axis = plt.subplots(figsize=(5.25, 3.35))
        palette = ["#1B9E77", "#7570B3", "#66A61E", "#A6761D", "#666666"]
        line_styles = ["-"]

    lookup = data.drop_duplicates("well_id").set_index("well_id")
    for index, well in enumerate(well_order):
        curve = data[data["well_id"].eq(well)].sort_values("elapsed_time_min")
        row = lookup.loc[well]
        axis.plot(
            curve["elapsed_time_min"],
            curve["normalized_PL_peak_intensity"],
            color=palette[index],
            linestyle=line_styles[index % len(line_styles)],
            marker="o" if not full_set else None,
            markersize=2.5,
            linewidth=1.35 if full_set else 1.55,
            label=f"{well}  {composition_label(row)}",
        )
    axis.axvline(189, color=COLORS["dark"], linewidth=0.85, linestyle=(0, (3, 2)))
    axis.text(186.0, 0.055, "189 min", ha="right", va="bottom", fontsize=8.0, color=COLORS["muted"])
    axis.set_xlim(0, 198)
    axis.set_ylim(0, 1.04)
    axis.set_xticks([0, 50, 100, 150])
    axis.set_yticks(np.linspace(0, 1, 6))
    axis.set_xlabel("Elapsed time after initial PL measurement (min)")
    axis.set_ylabel(r"Normalized PL peak intensity, $I(t)/I_{max}$")
    axis.grid(color="#E3E6EA", linewidth=0.55)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.text(
        0.02,
        0.035,
        "Each trace normalized independently",
        transform=axis.transAxes,
        fontsize=7.6,
        color=COLORS["muted"],
        ha="left",
        va="bottom",
    )
    axis.legend(
        title="Well   FA / BDA / PEA (%)",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        borderaxespad=0,
        ncol=2 if full_set else 1,
        fontsize=7.1 if full_set else 7.8,
        title_fontsize=7.4 if full_set else 8.0,
        handlelength=2.0,
        columnspacing=1.0,
    )
    figure.subplots_adjust(
        left=0.12,
        right=0.70 if full_set else 0.72,
        bottom=0.16,
        top=0.98,
    )
    save_figure(figure, stem)

def marker_sizes(
    values: pd.Series,
    minimum: float | None = None,
    maximum: float | None = None,
) -> np.ndarray:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
    lower = float(np.nanmin(numeric)) if minimum is None else float(minimum)
    upper = float(np.nanmax(numeric)) if maximum is None else float(maximum)
    span = upper - lower
    if span <= 0:
        return np.full(len(numeric), 70.0)
    return 42.0 + 105.0 * (numeric - lower) / span

def plot_matched_ternary(data: pd.DataFrame, stem: Path) -> None:
    figure, axis = plt.subplots(figsize=(6.5, 3.9))
    draw_ternary(axis)
    x, y = ternary_xy(data)
    intensities = pd.to_numeric(data["matched_PL_intensity"], errors="coerce")
    intensity_min = float(intensities.min())
    intensity_max = float(intensities.max())
    sizes = marker_sizes(intensities, intensity_min, intensity_max)
    values = pd.to_numeric(data["abs_delta_SP"], errors="coerce")
    scatter = axis.scatter(
        x,
        y,
        s=np.maximum(sizes, 48.0),
        marker="o",
        c=values,
        cmap="magma",
        vmin=float(values.min()),
        vmax=float(values.max()),
        edgecolors=COLORS["dark"],
        linewidths=0.65,
        zorder=3,
    )
    colorbar_axis = figure.add_axes([0.635, 0.18, 0.025, 0.67])
    colorbar = figure.colorbar(scatter, cax=colorbar_axis)
    colorbar.ax.set_title("photoKPFM\n" + r"$|\Delta SP|$ (V)", fontsize=7.5, pad=6)
    colorbar.ax.tick_params(labelsize=7.4, length=2.5)

    intensity_levels = np.quantile(intensities, [0.2, 0.5, 0.8])
    size_handles = []
    for value in intensity_levels:
        size = max(
            float(
                marker_sizes(
                    pd.Series([value]),
                    intensity_min,
                    intensity_max,
                )[0]
            ),
            48.0,
        )
        size_handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                markerfacecolor="#C8CBD0",
                markeredgecolor=COLORS["dark"],
                markersize=math.sqrt(size) * 0.65,
                label=rf"{value / 1000.0:.0f}$\times 10^3$",
            )
        )
    figure.legend(
        handles=size_handles,
        title="PL peak intensity at 189 min\n(a.u.)",
        loc="upper left",
        bbox_to_anchor=(0.72, 0.74),
        frameon=False,
        fontsize=7.2,
        title_fontsize=7.5,
        borderaxespad=0,
        handletextpad=0.5,
    )
    figure.subplots_adjust(left=0.02, right=0.59, bottom=0.05, top=0.99)
    save_figure(figure, stem)

def build_rank_table(data: pd.DataFrame) -> pd.DataFrame:
    ranked = data.copy()
    ranked["PL_rank"] = pd.to_numeric(ranked["matched_PL_intensity"], errors="coerce").rank(
        ascending=False, method="first"
    )
    ranked["photoKPFM_rank"] = pd.to_numeric(ranked["abs_delta_SP"], errors="coerce").rank(
        ascending=False, method="first"
    )
    ranked = ranked.dropna(subset=["PL_rank", "photoKPFM_rank"]).copy()
    ranked[["PL_rank", "photoKPFM_rank"]] = ranked[["PL_rank", "photoKPFM_rank"]].astype(int)
    ranked["rank_shift"] = ranked["photoKPFM_rank"] - ranked["PL_rank"]
    ranked["rank_category"] = "Other"
    ranked.loc[
        ranked[["PL_rank", "photoKPFM_rank"]].max(axis=1).le(4), "rank_category"
    ] = "Both in top 4"
    ranked.loc[ranked["rank_shift"].ge(5), "rank_category"] = "PL at least 5 ranks higher"
    ranked.loc[ranked["rank_shift"].le(-5), "rank_category"] = "photoKPFM at least 5 ranks higher"
    return ranked.sort_values(["PL_rank", "well_id"]).reset_index(drop=True)

def plot_rank_comparison(data: pd.DataFrame, stem: Path) -> pd.DataFrame:
    ranked = build_rank_table(data)
    y = np.arange(len(ranked))
    figure, axis = plt.subplots(figsize=(5.45, 4.35))
    for index in range(len(ranked)):
        if index % 2 == 0:
            axis.axhspan(index - 0.5, index + 0.5, color="#F6F7F9", linewidth=0, zorder=0)
    category_style = {
        "Both in top 4": (COLORS["teal"], 2.2),
        "PL at least 5 ranks higher": (COLORS["gold"], 2.2),
        "photoKPFM at least 5 ranks higher": (COLORS["red"], 2.2),
        "Other": (COLORS["gray"], 1.1),
    }
    for index, row in ranked.iterrows():
        color, width = category_style[str(row["rank_category"])]
        axis.plot(
            [row["PL_rank"], row["photoKPFM_rank"]],
            [index, index],
            color=color,
            linewidth=width,
            solid_capstyle="round",
            zorder=1,
        )
    axis.scatter(
        ranked["PL_rank"],
        y,
        marker="s",
        s=34,
        color=COLORS["blue"],
        edgecolor="white",
        linewidth=0.5,
        label="PL peak intensity rank",
        zorder=3,
    )
    axis.scatter(
        ranked["photoKPFM_rank"],
        y,
        marker="o",
        s=34,
        color=COLORS["orange"],
        edgecolor="white",
        linewidth=0.5,
        label=r"photoKPFM $|\Delta SP|$ rank",
        zorder=3,
    )
    labels = [
        f"{int(float(row.FAPbI3_pct))} / {int(float(row.BDAPbI4_pct))} / {int(float(row.PEA2PbI4_pct))}"
        for row in ranked.itertuples(index=False)
    ]
    axis.set_yticks(y)
    axis.set_yticklabels(labels, fontsize=7.7, fontfamily="DejaVu Sans Mono")
    axis.invert_yaxis()
    axis.set_xlim(0.5, len(ranked) + 0.5)
    axis.set_xticks(range(1, len(ranked) + 1, 2))
    axis.set_xlabel("Response rank (1 = strongest)")
    axis.tick_params(axis="y", length=0, pad=4)
    axis.grid(axis="x", color="#DDE1E6", linewidth=0.6)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.text(
        -0.015,
        1.015,
        "FA / BDA / PEA (%)",
        transform=axis.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.4,
        fontweight="bold",
        color=COLORS["dark"],
    )
    measurement_legend = axis.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.005),
        ncol=2,
        frameon=False,
        fontsize=7.5,
        handletextpad=0.4,
        columnspacing=1.0,
        borderaxespad=0,
    )
    axis.add_artist(measurement_legend)
    category_handles = [
        Line2D([], [], color=COLORS["teal"], linewidth=2.2, label="Both in top 4"),
        Line2D([], [], color=COLORS["gold"], linewidth=2.2, label=r"PL higher by $\geq$5 ranks"),
        Line2D([], [], color=COLORS["red"], linewidth=2.2, label=r"photoKPFM higher by $\geq$5 ranks"),
        Line2D([], [], color=COLORS["gray"], linewidth=1.2, label="Other"),
    ]
    axis.legend(
        handles=category_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=2,
        frameon=False,
        fontsize=7.1,
        handlelength=2.0,
        handletextpad=0.4,
        columnspacing=1.1,
        borderaxespad=0,
    )
    axis.text(
        0.99,
        0.985,
        f"n = {len(ranked)} matched compositions",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=7.2,
        color=COLORS["muted"],
    )
    figure.subplots_adjust(left=0.275, right=0.985, bottom=0.22, top=0.88)
    save_figure(figure, stem)
    return ranked

def plot_plate_map(source: pd.DataFrame, stem: Path) -> None:
    rows = list("ABCDEFGH")
    columns = list(range(1, 13))
    source = source.copy()
    source = source.set_index("Well")
    endpoints = {
        "PEA": np.asarray([213, 94, 94], dtype=float) / 255.0,
        "BDA": np.asarray([89, 161, 79], dtype=float) / 255.0,
        "FA": np.asarray([76, 120, 168], dtype=float) / 255.0,
    }
    colors = np.zeros((len(rows), len(columns), 3), dtype=float)
    for row_index, row_name in enumerate(rows):
        for column_index, column in enumerate(columns):
            record = source.loc[f"{row_name}{column}"]
            fa = float(record["FA_pct"])
            bda = float(record["BDA_pct"])
            pea = float(record["PEA_pct"])
            colors[row_index, column_index] = (fa * endpoints["FA"] + bda * endpoints["BDA"] + pea * endpoints["PEA"]) / 100.0
    figure, axis = plt.subplots(figsize=(7.0, 4.0))
    axis.imshow(colors, aspect="auto", interpolation="nearest")
    axis.set_xticks(np.arange(len(columns)), [str(value) for value in columns])
    axis.set_yticks(np.arange(len(rows)), rows)
    axis.xaxis.tick_top()
    axis.tick_params(axis="both", which="major", length=0, pad=4)
    axis.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=0.65, alpha=0.75)
    axis.tick_params(which="minor", bottom=False, left=False)
    for row_index, row_name in enumerate(rows):
        for column_index, column in enumerate(columns):
            record = source.loc[f"{row_name}{column}"]
            text = f"{int(float(record['FA_pct']))}/{int(float(record['BDA_pct']))}/{int(float(record['PEA_pct']))}"
            rgb = colors[row_index, column_index]
            luminance = float(0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2])
            axis.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                fontsize=6.0,
                fontweight="bold",
                color="white" if luminance < 0.49 else "#111111",
            )
    for spine in axis.spines.values():
        spine.set_visible(False)
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=endpoints[key], edgecolor="none", label=f"{key} endpoint")
        for key in ("FA", "BDA", "PEA")
    ]
    axis.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.045),
        ncol=3,
        frameon=False,
        fontsize=8.0,
        handlelength=1.3,
        columnspacing=1.6,
        borderaxespad=0,
    )
    axis.text(
        0.5,
        -0.145,
        "Cell labels: FA / BDA / PEA (%)",
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=7.8,
        color=COLORS["muted"],
    )
    figure.subplots_adjust(left=0.055, right=0.995, bottom=0.18, top=0.94)
    save_figure(figure, stem)

def build_layout_qa(main_dir: Path, supplementary_dir: Path) -> tuple[pd.DataFrame, int]:
    rows = []
    hashes: dict[str, list[str]] = {}
    for role, directory in (("main", main_dir), ("supplementary", supplementary_dir)):
        for svg in sorted(directory.glob("*.svg")):
            png = svg.with_suffix(".png")
            pdf = svg.with_suffix(".pdf")
            with Image.open(png) as image:
                width, height = image.size
            png_hash = sha256(png)
            hashes.setdefault(png_hash, []).append(svg.stem)
            rows.append(
                {
                    "stem": svg.stem,
                    "role": role,
                    "width_px": width,
                    "height_px": height,
                    "aspect_ratio": round(width / height, 4),
                    "png_600dpi": True,
                    "svg_present": svg.is_file(),
                    "pdf_present": pdf.is_file(),
                    "triplet_complete": png.is_file() and svg.is_file() and pdf.is_file(),
                    "minimum_dimension_ge_1200px": min(width, height) >= 1200,
                    "png_sha256": png_hash,
                }
            )
    duplicate_groups = sum(1 for stems in hashes.values() if len(stems) > 1)
    return pd.DataFrame(rows), duplicate_groups

def main() -> None:
    args = parse_args()
    data_dir = args.data_dir.resolve()
    output = args.out_dir.resolve()
    if not data_dir.is_dir():
        raise FileNotFoundError(data_dir)
    if output.exists():
        if not args.overwrite:
            raise FileExistsError(f"Output exists: {output}. Use --overwrite to replace it.")
        root = ROOT.resolve()
        if output == root or root not in output.parents:
            raise RuntimeError(f"Refusing to remove an unsafe output path: {output}")
        shutil.rmtree(output)
    main_dir = output / "main"
    supplementary_dir = output / "supplementary"
    source_dir = output / "source_data"
    audit_dir = output / "audit"
    for directory in (main_dir, supplementary_dir, source_dir, audit_dir):
        directory.mkdir(parents=True, exist_ok=True)

    configure_plotting()
    processed = data_dir / "processed"
    retention_path = processed / "pl_retention_metrics.csv"
    selection_path = processed / "representative_pl_trace_selection.csv"
    matched_path = processed / "matched_pl_photokpfm_metrics.csv"
    association_path = processed / "exploratory_association_statistics.csv"
    plate_path = processed / "plate_composition_map.csv"
    representative_trace_path = processed / "representative_pl_traces.csv"
    full_trace_path = processed / "all_matched_pl_traces.csv"
    source_paths = [
        retention_path,
        selection_path,
        matched_path,
        association_path,
        plate_path,
        representative_trace_path,
        full_trace_path,
    ]
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing source files:\n" + "\n".join(missing))

    retention = pd.read_csv(retention_path)
    selection = pd.read_csv(selection_path)
    matched = pd.read_csv(matched_path).merge(
        selection[["well_id", "temporal_behavior_class"]],
        on="well_id",
        how="left",
        validate="one_to_one",
    )
    associations = pd.read_csv(association_path)
    plate = pd.read_csv(plate_path)
    representative_traces = pd.read_csv(representative_trace_path)
    full_traces = pd.read_csv(full_trace_path)

    selected_wells = selected_well_order(selection)
    matched_wells = set(matched["well_id"].astype(str))
    plot_full_library_ternary(
        retention,
        matched_wells,
        selected_wells,
        "I_189",
        "PL peak intensity at 189 min (a.u.)",
        main_dir / "Figure4A_PL_189min_ternary",
        "viridis",
    )
    plot_full_library_ternary(
        retention,
        matched_wells,
        selected_wells,
        "PL_retention_189_over_Imax",
        r"PL retention, $I_{189}/I_{max}$",
        main_dir / "Figure4B_PL_retention_ternary",
        "viridis",
        0.0,
        1.0,
    )
    representative_error = validate_trace_endpoints(representative_traces)
    plot_trace_panel(
        representative_traces,
        selected_wells,
        supplementary_dir / "FigureS4_representative_PL_time_traces",
        full_set=False,
    )
    full_error = validate_trace_endpoints(full_traces)
    plot_matched_ternary(matched, main_dir / "Figure5A_matched_screening_space")
    ranks = plot_rank_comparison(matched, main_dir / "Figure5B_candidate_rank_comparison")
    plot_plate_map(plate, supplementary_dir / "FigureS3_plate_composition_map")

    copied_sources = {
        "pl_retention_metrics.csv": retention_path,
        "representative_pl_trace_selection.csv": selection_path,
        "matched_pl_photokpfm_metrics.csv": matched_path,
        "exploratory_association_statistics.csv": association_path,
        "plate_composition_map.csv": plate_path,
        "representative_pl_traces.csv": representative_trace_path,
        "all_matched_pl_traces.csv": full_trace_path,
    }
    for name, path in copied_sources.items():
        shutil.copy2(path, source_dir / name)
    ranks.to_csv(source_dir / "candidate_rank_comparison.csv", index=False)

    layout_qa, duplicate_png_groups = build_layout_qa(main_dir, supplementary_dir)
    layout_qa.to_csv(audit_dir / "figure_layout_qa.csv", index=False)
    retention_values = pd.to_numeric(retention["PL_retention_189_over_Imax"], errors="coerce")
    composition_sums = retention[["FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]].apply(
        pd.to_numeric, errors="coerce"
    ).sum(axis=1)
    association_row = associations[
        associations["kpfm_metric"].eq("photovoltage_abs_V")
        & associations["pl_metric"].eq("Top_PL_peak_intensity")
    ].iloc[0]
    try:
        data_label = data_dir.relative_to(ROOT).as_posix()
    except ValueError:
        data_label = str(data_dir)
    qa = {
        "data_directory": data_label,
        "retention_rows": int(len(retention)),
        "unique_nominal_compositions": int(
            retention[["FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]]
            .drop_duplicates()
            .shape[0]
        ),
        "composition_sum_100_all_rows": bool(np.allclose(composition_sums, 100.0)),
        "retention_missing": int(retention_values.isna().sum()),
        "retention_min": float(retention_values.min()),
        "retention_max": float(retention_values.max()),
        "retention_within_0_1": bool(retention_values.dropna().between(0, 1).all()),
        "matched_wells": int(len(matched)),
        "representative_wells": selected_wells,
        "representative_trace_endpoint_max_abs_error": representative_error,
        "all_13_trace_endpoint_max_abs_error": full_error,
        "trace_endpoint_validation_pass": bool(max(representative_error, full_error) < 5e-5),
        "association_n": int(association_row["n"]),
        "association_raw_spearman_rho": float(association_row["spearman_rho"]),
        "association_raw_spearman_q_bh": float(association_row["spearman_q_bh"]),
        "association_adjusted_spearman_rho": float(association_row["partial_spearman_rho"]),
        "association_adjusted_spearman_q_bh": float(association_row["partial_spearman_q_bh"]),
        "publication_figure_stems": int(len(layout_qa)),
        "publication_triplets_complete": bool(layout_qa["triplet_complete"].all()),
        "publication_png_minimum_dimension_ge_1200px": bool(
            layout_qa["minimum_dimension_ge_1200px"].all()
        ),
        "duplicate_publication_png_hash_groups": int(duplicate_png_groups),
    }
    (audit_dir / "publication_figure_qa.json").write_text(
        json.dumps(qa, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Publication figures: {output}")
    print(
        f"Figure stems: main={len(list(main_dir.glob('*.svg')))}, "
        f"supplementary={len(list(supplementary_dir.glob('*.svg')))}"
    )
    print(f"QA: {audit_dir / 'publication_figure_qa.json'}")


if __name__ == "__main__":
    main()
