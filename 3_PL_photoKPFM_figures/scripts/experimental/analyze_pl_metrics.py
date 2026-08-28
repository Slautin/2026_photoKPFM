"""Additional photoluminescence metric analysis.

This script is intentionally additive. It reads the existing raw PL export and
photoKPFM processed data, computes additional PL metrics and diagnostics, and
writes all outputs to a separate folder without touching the manuscript or the
current Figure 4/5 outputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
import csv
import json
import math
import re
import shutil
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[2]
RAW_PL_CSV = ROOT / "data" / "raw" / "pl_plate_reader_export.csv"
OUT = ROOT / "results" / "validation" / "pl_metrics"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from data_processing import (  # noqa: E402
    MATCHED_TIME_MIN,
    MATCHED_TIMEPOINT,
    build_all_data,
    generate_composition_table,
    ternary_xy,
)


ROWS = list("ABCDEFGH")
WELLS_ROW_MAJOR = [f"{row}{col}" for row in ROWS for col in range(1, 13)]

PL_RANGE = (455.0, 870.0)
METRIC_RANGE = (455.0, 870.0)
PEAK_SEARCH_RANGE = (500.0, 850.0)
BASELINE_LOW_RANGE = (455.0, 485.0)
BASELINE_HIGH_RANGE = (845.0, 870.0)
SMOOTH_WINDOW_NM = 9
LARGE_DROP_FRACTION = -0.25
TEMPORAL_RANGE_THRESHOLD_PCT = 30.0
PL_CEILING_FLAG = 99000.0


plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 10.5,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


@dataclass
class SpectrumBlock:
    read_number: int
    geometry: str
    timepoint: int
    elapsed_time_min: float
    timestamp_approx: str
    wavelengths: np.ndarray
    values: np.ndarray


def prepare_data(pl: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    initial = (
        pl.loc[pl["timepoint"].eq(1), ["well_id", "matched_PL_intensity"]]
        .rename(columns={"matched_PL_intensity": "PL_initial"})
    )
    result = matched.merge(initial, on="well_id", validate="one_to_one")
    result["PL_change_pct"] = (
        100.0
        * (result["matched_PL_intensity"] - result["PL_initial"])
        / result["PL_initial"]
    )
    pl_cut = pl.loc[
        pl["well_id"].isin(result["well_id"])
        & pl["elapsed_time_min"].le(MATCHED_TIME_MIN)
    ]
    variation = pl_cut.groupby("well_id")["matched_PL_intensity"].agg(
        PL_min="min", PL_max="max", PL_mean="mean", PL_std="std"
    )
    variation["PL_relative_range_pct"] = (
        100.0 * (variation["PL_max"] - variation["PL_min"]) / variation["PL_mean"]
    )
    variation["PL_cv_pct"] = 100.0 * variation["PL_std"] / variation["PL_mean"]
    result = result.merge(
        variation,
        left_on="well_id",
        right_index=True,
        validate="one_to_one",
    )
    result["temporal_behavior_class"] = np.where(
        result["PL_relative_range_pct"].gt(TEMPORAL_RANGE_THRESHOLD_PCT),
        "Pronounced temporal evolution",
        "Lower observed temporal variation",
    )
    result["PL_ceiling_flag"] = result["PL_max"].ge(PL_CEILING_FLAG)
    result["PL_rank"] = result["matched_PL_intensity"].rank(
        method="min", ascending=False
    ).astype(int)
    result["delta_SP_rank"] = result["abs_delta_SP"].rank(
        method="min", ascending=False
    ).astype(int)
    result["rank_difference"] = result["delta_SP_rank"] - result["PL_rank"]
    return result


def safe_float(value: str) -> float:
    text = str(value).strip()
    if not text or text.upper() in {"OVRFLW", "OVERFLOW", "OVER", "NAN"}:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def parse_start_datetime(rows: list[list[str]]) -> datetime | None:
    date_text = None
    time_text = None
    for row in rows[:20]:
        if row and row[0].strip() == "Date" and len(row) > 1:
            date_text = row[1].strip()
        if row and row[0].strip() == "Time" and len(row) > 1:
            time_text = row[1].strip()
    if not date_text or not time_text:
        return None
    for fmt in ("%m/%d/%Y %I:%M:%S %p", "%m/%d/%y %I:%M:%S %p"):
        try:
            return datetime.strptime(f"{date_text} {time_text}", fmt)
        except ValueError:
            continue
    return None


def read_geometry_from_number(read_number: int) -> str:
    mod = (read_number - 1) % 3
    if mod == 0:
        return "top"
    if mod == 1:
        return "bottom"
    return "absorbance"


def timepoint_from_read(read_number: int, geometry: str) -> int:
    if geometry == "top":
        return (read_number - 1) // 3 + 1
    if geometry == "bottom":
        return (read_number - 2) // 3 + 1
    return (read_number - 3) // 3 + 1


def parse_fluorescence_blocks() -> tuple[list[SpectrumBlock], datetime | None, int]:
    with RAW_PL_CSV.open(newline="", errors="replace") as handle:
        rows = list(csv.reader(handle))
    start_dt = parse_start_datetime(rows)
    starts: list[tuple[int, int]] = []
    for index, row in enumerate(rows):
        first = row[0].strip() if row else ""
        match = re.match(r"Read\s+(\d+):", first)
        if match:
            starts.append((index, int(match.group(1))))

    blocks: list[SpectrumBlock] = []
    skipped_empty = 0
    for block_index, (start, read_number) in enumerate(starts):
        geometry = read_geometry_from_number(read_number)
        if geometry == "absorbance":
            continue
        end = starts[block_index + 1][0] if block_index + 1 < len(starts) else len(rows)
        header_index = next(
            (j for j in range(start, end) if len(rows[j]) > 1 and rows[j][1].strip() == "Wavelength"),
            None,
        )
        if header_index is None:
            continue
        wells = [cell.strip() for cell in rows[header_index][2:] if cell.strip()]
        if wells != WELLS_ROW_MAJOR:
            raise ValueError(f"Unexpected well header for read {read_number}.")

        wavelengths = []
        values = []
        for row in rows[header_index + 1 : end]:
            if len(row) < 98:
                continue
            wl = safe_float(row[1])
            if not np.isfinite(wl):
                continue
            vals = [safe_float(v) for v in row[2:98]]
            wavelengths.append(wl)
            values.append(vals)
        if not wavelengths:
            continue
        wl_arr = np.asarray(wavelengths, dtype=float)
        val_arr = np.asarray(values, dtype=float)
        if np.isnan(val_arr).all():
            skipped_empty += 1
            continue
        timepoint = timepoint_from_read(read_number, geometry)
        elapsed = float((timepoint - 1) * 9)
        timestamp = ""
        if start_dt is not None:
            timestamp = (start_dt + timedelta(minutes=elapsed)).strftime("%Y-%m-%d %H:%M:%S")
        blocks.append(
            SpectrumBlock(
                read_number=read_number,
                geometry=geometry,
                timepoint=timepoint,
                elapsed_time_min=elapsed,
                timestamp_approx=timestamp,
                wavelengths=wl_arr,
                values=val_arr,
            )
        )
    return blocks, start_dt, skipped_empty


def baseline_correct(wl: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, list[str]]:
    flags: list[str] = []
    y = np.asarray(y, dtype=float)
    finite = np.isfinite(y)
    if finite.sum() < 10:
        return np.full_like(y, np.nan), np.nan, ["missing_or_too_few_points"]
    if finite.sum() < len(y):
        flags.append("missing_or_overflow_points")
    edge = (
        ((wl >= BASELINE_LOW_RANGE[0]) & (wl <= BASELINE_LOW_RANGE[1]))
        | ((wl >= BASELINE_HIGH_RANGE[0]) & (wl <= BASELINE_HIGH_RANGE[1]))
    )
    baseline = float(np.nanpercentile(y[edge & finite], 10)) if np.any(edge & finite) else float(np.nanpercentile(y[finite], 5))
    yc = y - baseline
    yc[yc < 0] = 0
    if np.nanmax(yc) <= 0:
        flags.append("nonpositive_after_baseline")
    if np.nanmax(y) >= 65000:
        flags.append("possible_saturation")
    return yc, baseline, flags


def odd_window(n: int, requested: int) -> int:
    requested = max(5, int(requested))
    if requested % 2 == 0:
        requested += 1
    max_window = n if n % 2 else n - 1
    return max(5, min(requested, max_window))


def spectrum_metrics(wl: np.ndarray, y: np.ndarray) -> dict[str, object]:
    yc, baseline, flags = baseline_correct(wl, y)
    metric_mask = (wl >= METRIC_RANGE[0]) & (wl <= METRIC_RANGE[1])
    peak_mask = (wl >= PEAK_SEARCH_RANGE[0]) & (wl <= PEAK_SEARCH_RANGE[1])
    if not np.any(metric_mask) or not np.any(peak_mask) or np.isnan(yc).all():
        return {
            "existing_peak_intensity_recomputed": np.nan,
            "integrated_pl_area": np.nan,
            "peak_wavelength_nm": np.nan,
            "fwhm_nm": np.nan,
            "baseline": baseline,
            "qc_flags": ";".join(sorted(set(flags + ["metric_failed"]))),
        }

    y_peak = yc[peak_mask]
    wl_peak = wl[peak_mask]
    if np.isfinite(y_peak).sum() < 10:
        flags.append("too_few_peak_points")
        return {
            "existing_peak_intensity_recomputed": np.nan,
            "integrated_pl_area": np.nan,
            "peak_wavelength_nm": np.nan,
            "fwhm_nm": np.nan,
            "baseline": baseline,
            "qc_flags": ";".join(sorted(set(flags))),
        }

    win = odd_window(len(y_peak), SMOOTH_WINDOW_NM)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_smooth = savgol_filter(np.nan_to_num(y_peak, nan=0.0), window_length=win, polyorder=2)
    peak_idx = int(np.nanargmax(y_smooth))
    peak_intensity = float(y_smooth[peak_idx])
    peak_wl = float(wl_peak[peak_idx])
    area = float(np.trapezoid(np.nan_to_num(yc[metric_mask], nan=0.0), wl[metric_mask]))
    fwhm = np.nan
    if peak_intensity > 0:
        try:
            widths = peak_widths(y_smooth, [peak_idx], rel_height=0.5)[0]
            step = float(np.nanmedian(np.diff(wl_peak)))
            fwhm = float(widths[0] * step)
            if not (1.0 <= fwhm <= 250.0):
                flags.append("fwhm_outside_reliable_range")
                fwhm = np.nan
        except Exception:
            flags.append("fwhm_failed")

    if peak_intensity < 100:
        flags.append("low_signal")
    return {
        "existing_peak_intensity_recomputed": peak_intensity,
        "integrated_pl_area": area,
        "peak_wavelength_nm": peak_wl,
        "fwhm_nm": fwhm,
        "baseline": baseline,
        "qc_flags": ";".join(sorted(set(flags))) if flags else "ok",
    }


def compute_metrics(blocks: list[SpectrumBlock], compositions: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    comp = compositions.set_index("well_id")
    for block in blocks:
        for well_index, well in enumerate(WELLS_ROW_MAJOR):
            metrics = spectrum_metrics(block.wavelengths, block.values[:, well_index])
            row = comp.loc[well]
            records.append(
                {
                    "well_id": well,
                    "FAPbI3_pct": row["FAPbI3_pct"],
                    "BDAPbI4_pct": row["BDAPbI4_pct"],
                    "PEA2PbI4_pct": row["PEA2PbI4_pct"],
                    "acquisition_number": block.timepoint,
                    "elapsed_time_min": block.elapsed_time_min,
                    "timestamp_approx": block.timestamp_approx,
                    "read_geometry": block.geometry,
                    "read_number": block.read_number,
                    **metrics,
                }
            )
    return pd.DataFrame(records)


def savefig(fig: plt.Figure, stem: str) -> None:
    for ext, kwargs in (("png", {"dpi": 450}), ("svg", {}), ("pdf", {})):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", pad_inches=0.08, facecolor="white", **kwargs)
    plt.close(fig)


def draw_ternary(ax) -> None:
    height = math.sqrt(3) / 2
    ax.plot([0, 1, 0.5, 0], [0, 0, height, 0], color="black", lw=1.1)
    for f in (0.2, 0.4, 0.6, 0.8):
        ax.plot([0.5 * f, 1 - 0.5 * f], [height * f] * 2, color="0.90", lw=0.5)
        ax.plot([f, 0.5 + 0.5 * f], [0, height * (1 - f)], color="0.90", lw=0.5)
        ax.plot([1 - f, 0.5 * (1 - f)], [0, height * (1 - f)], color="0.90", lw=0.5)
    ax.text(-0.04, -0.065, r"BDAPbI$_4$", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.text(1.04, -0.065, r"PEA$_2$PbI$_4$", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.text(0.5, height + 0.04, r"FAPbI$_3$", ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-0.09, 1.09)
    ax.set_ylim(-0.11, height + 0.11)
    ax.axis("off")


def ternary_map(frame: pd.DataFrame, value_col: str, stem: str, label: str, cmap: str = "viridis", vmin=None, vmax=None) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    draw_ternary(ax)
    x, y = ternary_xy(frame)
    vals = frame[value_col].to_numpy(float)
    sc = ax.scatter(x, y, c=vals, cmap=cmap, s=36, edgecolors="none", vmin=vmin, vmax=vmax)
    cb = fig.colorbar(sc, ax=ax, fraction=0.040, pad=0.025, shrink=0.82)
    cb.set_label(label)
    savefig(fig, stem)


def fisher_ci(r: float, n: int) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def correlation_summary(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 3:
        return {}
    pr, pp = pearsonr(x, y)
    sr, sp = spearmanr(x, y)
    lo, hi = fisher_ci(float(pr), n)
    loo = []
    for i in range(n):
        if n - 1 >= 3:
            loo.append(float(pearsonr(np.delete(x, i), np.delete(y, i))[0]))
    return {
        "n": n,
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "pearson_ci_low": lo,
        "pearson_ci_high": hi,
        "spearman_rho": float(sr),
        "spearman_p": float(sp),
        "leave_one_out_pearson_min": float(np.nanmin(loo)) if loo else np.nan,
        "leave_one_out_pearson_max": float(np.nanmax(loo)) if loo else np.nan,
    }


def plot_scatter(frame: pd.DataFrame, x_col: str, y_col: str, stem: str, x_label: str, y_label: str) -> None:
    fig, ax = plt.subplots(figsize=(4.9, 4.1))
    ax.scatter(frame[x_col], frame[y_col], s=54, c="#2F6B9A", edgecolors="black", linewidths=0.6)
    for row in frame.itertuples():
        ax.text(getattr(row, x_col), getattr(row, y_col), f" {row.well_id}", fontsize=7.5, va="center")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(color="0.90", lw=0.5)
    savefig(fig, stem)


def plot_rank_comparison(frame: pd.DataFrame, metric_rank_col: str, stem: str, metric_label: str) -> None:
    ordered = frame.sort_values(metric_rank_col).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(5.3, 4.2))
    ypos = np.arange(len(ordered))
    ax.hlines(ypos, ordered[metric_rank_col], ordered["delta_SP_rank"], color="0.75", lw=0.8)
    ax.scatter(ordered[metric_rank_col], ypos, marker="s", s=36, color="#4C78A8", label=metric_label, zorder=3)
    ax.scatter(ordered["delta_SP_rank"], ypos, marker="o", s=36, color="#F28E2B", label="|Delta SP| rank", zorder=3)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.well_id} ({int(r.FAPbI3_pct)}/{int(r.BDAPbI4_pct)}/{int(r.PEA2PbI4_pct)})" for r in ordered.itertuples()], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Response rank (1 = strongest)")
    ax.grid(axis="x", color="0.90", lw=0.5)
    ax.legend(frameon=True)
    savefig(fig, stem)


def make_timecourse_qc(metrics: pd.DataFrame, blocks: list[SpectrumBlock]) -> tuple[pd.DataFrame, dict[str, object]]:
    q_rows = []
    summary: dict[str, object] = {}
    for geom, g in metrics.groupby("read_geometry"):
        piv = g.pivot_table(index="well_id", columns="acquisition_number", values="existing_peak_intensity_recomputed")
        median_raw = piv.median(axis=0)
        norm = piv.div(piv.max(axis=1).replace(0, np.nan), axis=0)
        median_norm = norm.median(axis=0)
        frac_drop = piv.pct_change(axis=1).lt(LARGE_DROP_FRACTION).mean(axis=0)
        for tp in piv.columns:
            q_rows.append(
                {
                    "read_geometry": geom,
                    "acquisition_number": int(tp),
                    "elapsed_time_min": float((tp - 1) * 9),
                    "median_raw_peak_intensity": float(median_raw.loc[tp]),
                    "median_normalized_peak_intensity": float(median_norm.loc[tp]),
                    "fraction_large_drop_from_previous": float(frac_drop.loc[tp]) if tp in frac_drop else np.nan,
                }
            )
        max_drop_tp = int(frac_drop.idxmax())
        summary[f"{geom}_max_drop_timepoint"] = max_drop_tp
        summary[f"{geom}_max_drop_fraction"] = float(frac_drop.loc[max_drop_tp])
        summary[f"{geom}_median_drop_pct_at_max"] = float(piv.pct_change(axis=1)[max_drop_tp].median())

        fig, ax = plt.subplots(figsize=(6.2, 3.5))
        ax.plot(median_raw.index, median_raw.values, marker="o", lw=1.3)
        ax.set_xlabel("PL acquisition number")
        ax.set_ylabel("Median raw peak intensity")
        ax.set_title(f"{geom.title()} read: median raw PL peak intensity")
        ax.grid(color="0.90", lw=0.5)
        savefig(fig, f"timecourse_qc_{geom}_median_raw")

        fig, ax = plt.subplots(figsize=(6.2, 3.5))
        ax.plot(median_norm.index, median_norm.values, marker="o", lw=1.3)
        ax.set_xlabel("PL acquisition number")
        ax.set_ylabel("Median normalized peak intensity")
        ax.set_title(f"{geom.title()} read: median normalized PL")
        ax.grid(color="0.90", lw=0.5)
        savefig(fig, f"timecourse_qc_{geom}_median_normalized")

        fig, ax = plt.subplots(figsize=(6.2, 3.5))
        ax.plot(frac_drop.index, frac_drop.values, marker="o", lw=1.3, color="#A23B3B")
        ax.set_xlabel("PL acquisition number")
        ax.set_ylabel(f"Fraction wells with >{abs(LARGE_DROP_FRACTION):.0%} drop")
        ax.set_title(f"{geom.title()} read: simultaneous drop diagnostic")
        ax.grid(color="0.90", lw=0.5)
        savefig(fig, f"timecourse_qc_{geom}_drop_fraction")

        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        im = ax.imshow(norm.to_numpy(float), aspect="auto", interpolation="nearest", vmin=0, vmax=1, cmap="viridis")
        ax.set_xlabel("PL acquisition number")
        ax.set_ylabel("Well index")
        ax.set_xticks(np.arange(len(norm.columns))[::5])
        ax.set_xticklabels(norm.columns[::5])
        cb = fig.colorbar(im, ax=ax)
        cb.set_label("Normalized peak intensity")
        ax.set_title(f"{geom.title()} read: full-library normalized PL heatmap")
        savefig(fig, f"timecourse_qc_{geom}_normalized_heatmap")
    qc = pd.DataFrame(q_rows)
    return qc, summary


def plot_matched_trajectories(metrics: pd.DataFrame, matched_wells: list[str]) -> None:
    for geom in ("top", "bottom"):
        g = metrics[(metrics["read_geometry"].eq(geom)) & (metrics["well_id"].isin(matched_wells))]
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        for well, sub in g.groupby("well_id"):
            s = sub.sort_values("acquisition_number")
            y = s["existing_peak_intensity_recomputed"].to_numpy(float)
            if np.nanmax(y) > 0:
                y = y / np.nanmax(y)
            ax.plot(s["elapsed_time_min"], y, lw=1.0, alpha=0.85, label=well)
        ax.axvline(MATCHED_TIME_MIN, color="0.25", ls="--", lw=0.9)
        ax.set_xlabel("Elapsed time after first PL measurement (min)")
        ax.set_ylabel("Peak intensity / max")
        ax.set_title(f"{geom.title()} read: all 13 matched normalized PL trajectories")
        ax.legend(ncol=3, fontsize=6.5, frameon=True)
        ax.grid(color="0.90", lw=0.5)
        savefig(fig, f"timecourse_qc_{geom}_matched_13_trajectories")


def representative_spectra_before_after(blocks: list[SpectrumBlock], metrics: pd.DataFrame, qc: pd.DataFrame) -> None:
    for geom in ("top", "bottom"):
        qg = qc[qc["read_geometry"].eq(geom)]
        if qg.empty:
            continue
        tp = int(qg.sort_values("fraction_large_drop_from_previous", ascending=False).iloc[0]["acquisition_number"])
        if tp <= 1:
            continue
        candidates = metrics[
            (metrics["read_geometry"].eq(geom))
            & (metrics["acquisition_number"].isin([tp - 1, tp]))
        ].pivot_table(index="well_id", columns="acquisition_number", values="existing_peak_intensity_recomputed")
        if tp not in candidates or (tp - 1) not in candidates:
            continue
        candidates["drop"] = (candidates[tp] - candidates[tp - 1]) / candidates[tp - 1].replace(0, np.nan)
        wells = candidates.sort_values("drop").head(4).index.tolist()
        block_prev = next((b for b in blocks if b.geometry == geom and b.timepoint == tp - 1), None)
        block_now = next((b for b in blocks if b.geometry == geom and b.timepoint == tp), None)
        if block_prev is None or block_now is None:
            continue
        for well in wells:
            fig, ax = plt.subplots(figsize=(5.8, 4.2))
            wi = WELLS_ROW_MAJOR.index(well)
            ax.plot(block_prev.wavelengths, block_prev.values[:, wi], label=f"Before t{tp-1}", lw=1.1)
            ax.plot(block_now.wavelengths, block_now.values[:, wi], label=f"After t{tp}", lw=1.1)
            ax.set_title(
                f"{geom.title()} read {well}: spectra before and after largest plate-level drop",
                fontsize=10,
            )
            ax.set_xlim(455, 870)
            ax.grid(color="0.92", lw=0.5)
            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Raw PL signal")
            ax.legend(frameon=True, fontsize=7)
            savefig(fig, f"timecourse_qc_{geom}_{well}_spectra_before_after")


def top_bottom_comparisons(metrics: pd.DataFrame, matched_wells: list[str]) -> pd.DataFrame:
    matched_metrics = metrics[
        (metrics["well_id"].isin(matched_wells))
        & (metrics["acquisition_number"].eq(MATCHED_TIMEPOINT))
    ].copy()
    wide = matched_metrics.pivot_table(
        index=["well_id", "FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"],
        columns="read_geometry",
        values=["existing_peak_intensity_recomputed", "integrated_pl_area", "peak_wavelength_nm", "fwhm_nm"],
    )
    wide.columns = [f"{metric}_{geom}" for metric, geom in wide.columns]
    wide = wide.reset_index()
    for metric, label in [
        ("existing_peak_intensity_recomputed", "Peak intensity"),
        ("integrated_pl_area", "Integrated PL area"),
        ("peak_wavelength_nm", "Peak wavelength (nm)"),
    ]:
        xcol = f"{metric}_top"
        ycol = f"{metric}_bottom"
        if xcol in wide and ycol in wide:
            fig, ax = plt.subplots(figsize=(4.4, 4.2))
            ax.scatter(wide[xcol], wide[ycol], s=48, edgecolors="black", linewidths=0.6)
            lo = np.nanmin([wide[xcol].min(), wide[ycol].min()])
            hi = np.nanmax([wide[xcol].max(), wide[ycol].max()])
            ax.plot([lo, hi], [lo, hi], color="0.5", ls="--", lw=0.9)
            ax.set_xlabel(f"Top read {label}")
            ax.set_ylabel(f"Bottom read {label}")
            ax.grid(color="0.90", lw=0.5)
            savefig(fig, f"additional_top_vs_bottom_{metric}")
    return wide


def spectral_feature_inventory(blocks: list[SpectrumBlock], compositions: pd.DataFrame, matched_wells: list[str]) -> pd.DataFrame:
    records = []
    comp = compositions.set_index("well_id")
    selected_blocks = [b for b in blocks if b.timepoint in {1, MATCHED_TIMEPOINT, max(bb.timepoint for bb in blocks if bb.geometry == b.geometry)}]
    for block in selected_blocks:
        for well in matched_wells:
            wi = WELLS_ROW_MAJOR.index(well)
            yc, _, flags = baseline_correct(block.wavelengths, block.values[:, wi])
            mask = (block.wavelengths >= PEAK_SEARCH_RANGE[0]) & (block.wavelengths <= PEAK_SEARCH_RANGE[1])
            y = yc[mask]
            wl = block.wavelengths[mask]
            if np.nanmax(y) <= 0 or np.isfinite(y).sum() < 10:
                continue
            smooth = savgol_filter(np.nan_to_num(y, nan=0), odd_window(len(y), 9), 2)
            prom = max(np.nanmax(smooth) * 0.08, 50)
            peaks, props = find_peaks(smooth, prominence=prom, distance=18)
            row = comp.loc[well]
            for p, prom_val in zip(peaks, props.get("prominences", [])):
                records.append(
                    {
                        "well_id": well,
                        "FAPbI3_pct": row["FAPbI3_pct"],
                        "BDAPbI4_pct": row["BDAPbI4_pct"],
                        "PEA2PbI4_pct": row["PEA2PbI4_pct"],
                        "read_geometry": block.geometry,
                        "acquisition_number": block.timepoint,
                        "elapsed_time_min": block.elapsed_time_min,
                        "approx_peak_position_nm": float(wl[p]),
                        "relative_peak_height": float(smooth[p] / np.nanmax(smooth)),
                        "prominence": float(prom_val),
                        "confidence_flag": "descriptive_only_no_phase_assignment",
                        "qc_flags": ";".join(flags) if flags else "ok",
                    }
                )
    return pd.DataFrame(records)


def write_timepoint_report(start_dt: datetime | None, matched: pd.DataFrame) -> None:
    matched_rows = matched[["well_id", "matched_PL_time_min", "acquisition_iteration"]].copy()
    timestamp = ""
    if start_dt is not None:
        timestamp = (start_dt + timedelta(minutes=MATCHED_TIME_MIN)).strftime("%Y-%m-%d %H:%M:%S")
    text = f"""# Timepoint 22 Verification

- Timepoint 22 is the 22nd Top PL acquisition in the existing analysis convention.
- PL timepoint indexing begins at 1, not 0.
- The elapsed time assigned by the existing analysis is `(22 - 1) x 9 min = {MATCHED_TIME_MIN} min` after the first Top PL measurement.
- Experiment start timestamp from the Cytation export: {start_dt.strftime('%Y-%m-%d %H:%M:%S') if start_dt else 'not available'}.
- Approximate timestamp for Top PL timepoint 22 using the 9 min spacing convention: {timestamp or 'not available'}.
- The exact per-read timestamp is not present in the CSV export; the timestamp above is therefore approximate.
- The existing workflow matches photoKPFM to Top PL timepoint 22 at {MATCHED_TIME_MIN} min. The photoKPFM HDF5 contains acquisition order and file IDs, but not absolute clock timestamps, so exact temporal overlap cannot be independently verified from the present files.

## Matched photoKPFM wells

{markdown_table(matched_rows)}
"""
    (OUT / "timepoint_22_verification.md").write_text(text, encoding="utf-8")


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "_No rows._"
    display = frame.copy()
    for col in display.columns:
        display[col] = display[col].map(lambda v: "" if pd.isna(v) else f"{v:.4g}" if isinstance(v, float) else str(v))
    header = "| " + " | ".join(display.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(display.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in display.astype(str).to_numpy()]
    return "\n".join([header, sep, *rows])


def write_qc_summary(qc: pd.DataFrame, summary: dict[str, object]) -> None:
    lines = ["# PL Timecourse QC Summary", ""]
    for geom in ("top", "bottom"):
        lines.append(f"## {geom.title()} read")
        if f"{geom}_max_drop_timepoint" in summary:
            tp = summary[f"{geom}_max_drop_timepoint"]
            frac = summary[f"{geom}_max_drop_fraction"]
            med = summary[f"{geom}_median_drop_pct_at_max"]
            lines.append(f"- Largest plate-wide large-drop fraction occurs at acquisition {tp}, elapsed time {(int(tp)-1)*9} min.")
            lines.append(f"- Fraction of wells with >25% drop from previous acquisition: {frac:.3f}.")
            lines.append(f"- Median previous-step change at that acquisition: {med:.3f}.")
            if frac > 0.40:
                lines.append("- This pattern is more consistent with a plate-wide step or measurement-condition change than a purely composition-specific response.")
            else:
                lines.append("- The large-drop diagnostic is not strongly plate-wide by this threshold; composition-specific behavior remains plausible.")
        lines.append("")
    lines += [
        "## Interpretation limits",
        "- These diagnostics do not automatically label the behavior as degradation.",
        "- A plate-wide step may indicate measurement conditions, focus/position, illumination history, or another unresolved acquisition factor.",
        "- Raw instrument logs would be required to verify interruptions, plate movement, or focus changes.",
    ]
    (OUT / "pl_timecourse_qc_summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_spectral_summary(inv: pd.DataFrame) -> None:
    if inv.empty:
        text = "# Additional Spectral Feature Summary\n\nNo robust multi-feature inventory entries were generated by the conservative descriptive peak finder.\n"
    else:
        counts = inv.groupby(["read_geometry", "acquisition_number"]).size().reset_index(name="feature_count")
        text = f"""# Additional Spectral Feature Summary

This is a descriptive inventory only. No PL features are assigned to structural phases or n-values.

## Inventory size

- Total descriptive feature entries: {len(inv)}
- Wells surveyed: {inv['well_id'].nunique()}

## Feature-count summary

{markdown_table(counts)}

## Interpretation limits

- Peak positions are approximate and derived from baseline-corrected raw spectra.
- Multiple features or shoulders should not be assigned to specific phases without structural validation.
"""
    (OUT / "additional_spectral_feature_summary.md").write_text(text, encoding="utf-8")


def write_final_summary(
    qc_summary: dict[str, object],
    rank_diff: pd.DataFrame,
    corr_stats: pd.DataFrame,
    top_bottom: pd.DataFrame,
    feature_inv: pd.DataFrame,
) -> None:
    top_drop = qc_summary.get("top_max_drop_fraction", np.nan)
    area_stats = corr_stats[corr_stats["comparison"].eq("integrated_area_vs_abs_delta_SP")]
    peak_stats = corr_stats[corr_stats["comparison"].eq("peak_intensity_vs_abs_delta_SP")]
    rank_changed = int((rank_diff["abs_rank_difference_area_vs_peak"] >= 3).sum()) if not rank_diff.empty else 0
    text = f"""# Additional PL Analysis Summary

This analysis package supplements the current manuscript figures and does not replace the reported conclusions.

## 1. Timepoint 22

Timepoint 22 is the 22nd Top PL acquisition. The existing analysis indexes PL timepoints starting at 1. Using the existing 9 min spacing convention, timepoint 22 corresponds to {MATCHED_TIME_MIN} min after the first Top PL measurement.

## 2. Late-time PL decrease

The largest Top-read large-drop fraction was {top_drop:.3f}. See `pl_timecourse_qc_summary.md` and `pl_timecourse_qc.csv` for the acquisition-specific diagnostics. The present data can identify whether the drop is plate-wide by PL signal statistics, but cannot by itself prove the underlying cause.

## 3. Measurement artifact evidence

The additional analysis checks for simultaneous plate-wide steps and representative before/after spectra. The current files do not contain instrument event logs, focus metadata, or plate-movement records, so any artifact assignment remains unresolved.

## 4. Integrated PL area versus peak intensity

Integrated area was calculated with a baseline-subtracted spectrum over {METRIC_RANGE[0]:.0f}-{METRIC_RANGE[1]:.0f} nm. {rank_changed} wells change rank by at least 3 positions between peak intensity and integrated area at the matched Top-read timepoint.

## 5. Integrated area and PL-photoKPFM statistics

Peak-intensity versus |Delta SP| and integrated-area versus |Delta SP| statistics are reported in `additional_pl_metric_photoKPFM_statistics.csv`.

Peak-intensity comparison:
{markdown_table(peak_stats) if not peak_stats.empty else 'not available'}

Integrated-area comparison:
{markdown_table(area_stats) if not area_stats.empty else 'not available'}

## 6. Peak wavelength trends

Peak wavelength maps and peak-wavelength versus |Delta SP| comparison plots are provided. These are optical descriptors only and are not phase assignments.

## 7. Top-read versus bottom-read differences

Top/bottom matched metrics are reported in `additional_top_vs_bottom_pl_metrics.csv`. The comparison should be described as read-geometry-dependent optical differences, not proof of vertical phase segregation.

## 8. Multiple spectral features

The descriptive feature inventory contains {len(feature_inv)} entries across {feature_inv['well_id'].nunique() if not feature_inv.empty else 0} wells. No features are assigned to specific phases or n-values.

## 9. Reliable observations

- Timepoint indexing and elapsed-time convention are reproducible from the existing script.
- Peak intensity, integrated area, peak wavelength, and FWHM descriptors are reproducibly calculated from raw spectra with documented baseline correction.
- PhotoKPFM comparisons use the existing matched 13-well subset and |Delta SP| definition.

## 10. Ambiguous observations

- The cause of any late-time PL decrease cannot be proven without instrument logs or independent measurement metadata.
- Optical peak shifts and shoulders cannot be structurally assigned from PL alone.

## 11. Supporting-information outputs

- `timecourse_qc_*`
- `integrated_area_map.*`
- `peak_wavelength_map.*`
- `additional_integrated_area_vs_photoKPFM.*`
- `additional_top_vs_bottom_*.png`

## 12. Questions not fully answerable from the present data

- Whether a late-time step is definitively caused by plate motion, focus changes, or another instrument event.
- Whether spectral shoulders correspond to specific structural phases or n-values.
- Whether top/bottom differences prove vertical phase segregation.
"""
    (OUT / "additional_analysis_summary.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    compositions, pl_existing, raw_matched = build_all_data()
    matched = prepare_data(pl_existing, raw_matched)
    matched_wells = matched["well_id"].tolist()

    blocks, start_dt, skipped_empty = parse_fluorescence_blocks()
    metrics = compute_metrics(blocks, compositions)
    metrics.to_csv(OUT / "additional_pl_metrics.csv", index=False)


    top_match = metrics[(metrics["read_geometry"].eq("top")) & (metrics["acquisition_number"].eq(MATCHED_TIMEPOINT))].copy()
    top_match = top_match.merge(
        pl_existing[pl_existing["timepoint"].eq(MATCHED_TIMEPOINT)][["well_id", "matched_PL_intensity"]],
        on="well_id",
        how="left",
    )
    ternary_map(top_match, "integrated_pl_area", "integrated_area_map", "Integrated PL area (baseline-subtracted)", "viridis")
    ternary_map(top_match, "peak_wavelength_nm", "peak_wavelength_map", "Peak wavelength (nm)", "plasma")
    reliable_fwhm = top_match["fwhm_nm"].notna().mean() > 0.60
    if reliable_fwhm:
        ternary_map(top_match, "fwhm_nm", "fwhm_map", "FWHM (nm)", "magma")


    corr_area_peak = correlation_summary(
        top_match["matched_PL_intensity"].to_numpy(float),
        top_match["integrated_pl_area"].to_numpy(float),
    )
    plot_scatter(
        top_match,
        "matched_PL_intensity",
        "integrated_pl_area",
        "integrated_area_scatter",
        "Existing peak intensity at 189 min",
        "Integrated PL area",
    )
    top_match["peak_intensity_rank"] = top_match["matched_PL_intensity"].rank(ascending=False, method="min")
    top_match["integrated_area_rank"] = top_match["integrated_pl_area"].rank(ascending=False, method="min")
    top_match["rank_difference_area_vs_peak"] = top_match["integrated_area_rank"] - top_match["peak_intensity_rank"]
    top_match["abs_rank_difference_area_vs_peak"] = top_match["rank_difference_area_vs_peak"].abs()
    rank_diff = top_match.sort_values("abs_rank_difference_area_vs_peak", ascending=False)
    rank_diff.to_csv(OUT / "additional_pl_metric_rank_comparison.csv", index=False)


    matched_metrics = top_match[top_match["well_id"].isin(matched_wells)].merge(
        matched[["well_id", "abs_delta_SP", "delta_SP", "delta_SP_se", "FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]],
        on=["well_id", "FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"],
        how="left",
    )
    stats_rows = []
    for comp_name, col in [
        ("peak_intensity_vs_abs_delta_SP", "matched_PL_intensity"),
        ("integrated_area_vs_abs_delta_SP", "integrated_pl_area"),
        ("peak_wavelength_vs_abs_delta_SP", "peak_wavelength_nm"),
    ]:
        stats = correlation_summary(matched_metrics[col].to_numpy(float), matched_metrics["abs_delta_SP"].to_numpy(float))
        stats_rows.append({"comparison": comp_name, **stats})
    corr_stats = pd.DataFrame(stats_rows)
    corr_stats["peak_intensity_area_pearson_r"] = corr_area_peak.get("pearson_r", np.nan)
    corr_stats["peak_intensity_area_spearman_rho"] = corr_area_peak.get("spearman_rho", np.nan)
    corr_stats.to_csv(OUT / "additional_pl_metric_photoKPFM_statistics.csv", index=False)

    matched_metrics["integrated_area_rank"] = matched_metrics["integrated_pl_area"].rank(ascending=False, method="min")
    matched_metrics["peak_intensity_rank"] = matched_metrics["matched_PL_intensity"].rank(ascending=False, method="min")
    matched_metrics["delta_SP_rank"] = matched_metrics["abs_delta_SP"].rank(ascending=False, method="min")
    matched_metrics["integrated_area_rank_disagreement"] = matched_metrics["delta_SP_rank"] - matched_metrics["integrated_area_rank"]
    matched_metrics["peak_intensity_rank_disagreement"] = matched_metrics["delta_SP_rank"] - matched_metrics["peak_intensity_rank"]
    matched_metrics.to_csv(OUT / "additional_matched_pl_photokpfm_metrics.csv", index=False)
    plot_scatter(
        matched_metrics,
        "integrated_pl_area",
        "abs_delta_SP",
        "additional_integrated_area_vs_photoKPFM",
        "Integrated PL area at matched timepoint",
        "Relative |Delta SP| (V)",
    )
    plot_scatter(
        matched_metrics,
        "peak_wavelength_nm",
        "abs_delta_SP",
        "additional_peak_wavelength_vs_photoKPFM",
        "Peak wavelength (nm)",
        "Relative |Delta SP| (V)",
    )
    plot_rank_comparison(matched_metrics, "integrated_area_rank", "additional_integrated_area_rank_comparison", "Integrated area rank")


    qc, qc_summary = make_timecourse_qc(metrics, blocks)
    qc.to_csv(OUT / "pl_timecourse_qc.csv", index=False)
    plot_matched_trajectories(metrics, matched_wells)
    representative_spectra_before_after(blocks, metrics, qc)


    top_bottom = top_bottom_comparisons(metrics, matched_wells)
    top_bottom.to_csv(OUT / "additional_top_vs_bottom_pl_metrics.csv", index=False)


    feature_inv = spectral_feature_inventory(blocks, compositions, matched_wells)
    feature_inv.to_csv(OUT / "additional_spectral_feature_inventory.csv", index=False)

    write_timepoint_report(start_dt, matched)
    write_qc_summary(qc, qc_summary)
    write_spectral_summary(feature_inv)
    write_final_summary(qc_summary, rank_diff, corr_stats, top_bottom, feature_inv)

    assumptions = {
        "raw_pl_csv": str(RAW_PL_CSV),
        "project_root": str(ROOT),
        "output_folder": str(OUT),
        "matched_timepoint": MATCHED_TIMEPOINT,
        "matched_time_min_existing_convention": MATCHED_TIME_MIN,
        "timepoint_indexing": "starts at 1",
        "baseline_low_range_nm": BASELINE_LOW_RANGE,
        "baseline_high_range_nm": BASELINE_HIGH_RANGE,
        "metric_range_nm": METRIC_RANGE,
        "peak_search_range_nm": PEAK_SEARCH_RANGE,
        "smoothing_window_nm": SMOOTH_WINDOW_NM,
        "large_drop_fraction_threshold": LARGE_DROP_FRACTION,
        "skipped_fully_empty_pl_blocks": skipped_empty,
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
    }
    (OUT / "additional_analysis_reproducibility.json").write_text(json.dumps(assumptions, indent=2), encoding="utf-8")
    shutil.copy2(Path(__file__), OUT / "additional_pl_metrics_analysis.py")
    print(f"Saved additional PL analysis package to: {OUT}")


if __name__ == "__main__":
    main()
