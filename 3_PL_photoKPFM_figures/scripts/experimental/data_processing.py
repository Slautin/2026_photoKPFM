"""Shared, data-driven helpers for manuscript Figures 4 and 5.

Inputs are the original Cytation5 CSV and raw photoKPFM HDF5 file.  No
measurement values are embedded in this module.  The 96-well composition
lookup is regenerated with the exact 8% ternary-grid algorithm stored in the
original photoKPFM notebook, avoiding a known F12 transcription error in a
later copied table.
"""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
PL_CSV = ROOT / "data" / "raw" / "pl_plate_reader_export.csv"
KPFM_H5 = ROOT / "data" / "raw" / "photokpfm_measurements.h5"
OUTDIR = ROOT / "results" / "intermediate"

ROWS = list("ABCDEFGH")
WELLS_ROW_MAJOR = [f"{row}{col}" for row in ROWS for col in range(1, 13)]
MATCHED_TIMEPOINT = 22
TIME_STEP_MIN = 9
MATCHED_TIME_MIN = (MATCHED_TIMEPOINT - 1) * TIME_STEP_MIN
PL_WL_MIN = 450.0
PL_WL_MAX = 850.0


def require_inputs() -> None:
    missing = [str(path) for path in (PL_CSV, KPFM_H5) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required source file(s) missing:\n" + "\n".join(missing))
    OUTDIR.mkdir(parents=True, exist_ok=True)


def generate_composition_table() -> pd.DataFrame:
    """Reproduce the original notebook's 8% grid and column-major plate order."""
    compositions: list[tuple[int, int, int]] = []
    for k in range(100 // 8 + 1):
        pea = 100 - 8 * k
        for bda_step in range(k, -1, -1):
            compositions.append((pea, 8 * bda_step, 8 * (k - bda_step)))
    if len(compositions) < 96:
        compositions.extend(compositions[-(96 - len(compositions)) :])
    compositions = compositions[:96]

    wells_column_major = [f"{row}{col}" for col in range(1, 13) for row in ROWS]
    frame = pd.DataFrame(compositions, columns=["PEA2PbI4_pct", "BDAPbI4_pct", "FAPbI3_pct"])
    frame.insert(0, "well_id", wells_column_major)
    totals = frame[["FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]].sum(axis=1)
    if len(frame) != 96 or not np.allclose(totals, 100.0, atol=1e-8):
        raise ValueError("Composition validation failed: expected 96 rows summing to 100%.")
    if frame["well_id"].duplicated().any():
        raise ValueError("Composition validation failed: duplicate well IDs.")
    return frame


def _safe_float(value: str) -> float:
    text = str(value).strip()
    if not text or text.upper() in {"OVRFLW", "OVERFLOW", "OVER", "NAN"}:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def parse_top_pl_peaks() -> pd.DataFrame:
    """Parse top-read PL blocks and retain the maximum signal from 450–850 nm."""
    with PL_CSV.open(newline="", errors="replace") as handle:
        rows = list(csv.reader(handle))
    starts = []
    for index, row in enumerate(rows):
        first = row[0].strip() if row else ""
        match = re.match(r"Read\s+(\d+):(.*)", first)
        if match:
            starts.append((index, int(match.group(1)), match.group(2).strip()))
    if not starts:
        raise ValueError("No Cytation5 read blocks found in the PL CSV.")

    records: list[dict[str, object]] = []
    top_number = 0
    for block_index, (start, read_number, title) in enumerate(starts):
        if (read_number - 1) % 3 != 0:
            continue
        end = starts[block_index + 1][0] if block_index + 1 < len(starts) else len(rows)
        header_index = next(
            (j for j in range(start, end) if len(rows[j]) > 1 and rows[j][1].strip() == "Wavelength"),
            None,
        )
        if header_index is None:
            raise ValueError(f"Missing wavelength header in top-read block {read_number}.")
        wells = [cell.strip() for cell in rows[header_index][2:] if cell.strip()]
        if wells != WELLS_ROW_MAJOR:
            raise ValueError(f"Unexpected or incomplete well header in top-read block {read_number}.")
        peaks = np.full(96, np.nan)
        for row in rows[header_index + 1 : end]:
            if len(row) < 2:
                continue
            wavelength = _safe_float(row[1])
            if not np.isfinite(wavelength) or not PL_WL_MIN <= wavelength <= PL_WL_MAX:
                continue
            values = np.asarray([_safe_float(v) for v in row[2:98]], dtype=float)
            peaks = np.fmax(peaks, values)
        if np.isnan(peaks).all():
            if top_number < MATCHED_TIMEPOINT:
                raise ValueError(f"Fully empty top-read block encountered before timepoint {MATCHED_TIMEPOINT}.")
            print(f"Excluding fully empty trailing top-read block {read_number}; no measured values were present.")
            continue
        top_number += 1
        if np.isnan(peaks).any():
            bad = [w for w, value in zip(wells, peaks) if not np.isfinite(value)]
            raise ValueError(f"Missing PL peak values in timepoint {top_number}: {bad}")
        for well, peak in zip(wells, peaks):
            records.append(
                {
                    "well_id": well,
                    "timepoint": top_number,
                    "elapsed_time_min": (top_number - 1) * TIME_STEP_MIN,
                    "matched_PL_intensity": float(peak),
                    "read_number": read_number,
                    "read_title": title,
                }
            )
    frame = pd.DataFrame(records)
    if frame.empty or frame.duplicated(["well_id", "timepoint"]).any():
        raise ValueError("PL parsing produced no rows or duplicate well/timepoint rows.")
    if frame["timepoint"].max() < MATCHED_TIMEPOINT:
        raise ValueError(f"PL data contain only {frame['timepoint'].max()} top-read timepoints.")
    return frame


def _gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * sigma**2))


def _histogram_center(values: np.ndarray) -> tuple[float, float]:
    """Return fitted Gaussian center and its covariance-derived standard error."""
    counts, edges = np.histogram(values, bins=30, range=(-2.5, -0.5))
    centers = 0.5 * (edges[:-1] + edges[1:])
    if counts.max() <= 0:
        raise ValueError("photoKPFM histogram is empty in the specified analysis range.")
    params, covariance = curve_fit(
        _gaussian,
        centers,
        counts,
        p0=[counts.max(), float(np.mean(values)), max(float(np.std(values)), 1e-6)],
        maxfev=5000,
    )
    center_se = float(np.sqrt(max(covariance[1, 1], 0.0)))
    if not np.isfinite(center_se):
        raise ValueError("Non-finite Gaussian-center uncertainty.")
    return float(params[1]), center_se


def plate_index_to_well(index: int) -> str:
    if not 0 <= index < 96:
        raise ValueError(f"Invalid zero-based plate index: {index}")
    return f"{ROWS[index % 8]}{index // 8 + 1}"


def load_kpfm_acquisitions() -> pd.DataFrame:
    """Fit SP centers for every raw dark/light pair; retain repeated acquisitions."""
    with h5py.File(KPFM_H5, "r") as handle:
        required = {"idx", "dark_data", "light_data", "dark_fn", "light_fn"}
        missing = required.difference(handle.keys())
        if missing:
            raise KeyError(f"Required HDF5 datasets missing: {sorted(missing)}")
        idx = np.asarray(handle["idx"], dtype=int)
        dark = np.asarray(handle["dark_data"], dtype=float)
        light = np.asarray(handle["light_data"], dtype=float)
        dark_names = np.asarray(handle["dark_fn"]).astype(str)
        light_names = np.asarray(handle["light_fn"]).astype(str)
    lengths = {len(idx), len(dark), len(light), len(dark_names), len(light_names)}
    if len(lengths) != 1:
        raise ValueError("HDF5 acquisition arrays have inconsistent lengths.")

    rows = []
    for iteration, (plate_idx, dark_image, light_image, dark_name, light_name) in enumerate(
        zip(idx, dark, light, dark_names, light_names), start=1
    ):
        dark_values = dark_image[10:-10, 10:-10].ravel()
        light_values = light_image[10:-10, 10:-10].ravel()
        sp_dark, sp_dark_se = _histogram_center(dark_values)
        sp_light, sp_light_se = _histogram_center(light_values)
        delta_sp_se = math.sqrt(sp_dark_se**2 + sp_light_se**2)
        rows.append(
            {
                "well_id": plate_index_to_well(int(plate_idx)),
                "plate_index": int(plate_idx),
                "SP_dark": sp_dark,
                "SP_light": sp_light,
                "SP_dark_se": sp_dark_se,
                "SP_light_se": sp_light_se,
                "delta_SP": sp_light - sp_dark,
                "delta_SP_se": delta_sp_se,
                "abs_delta_SP": abs(sp_light - sp_dark),
                "acquisition_iteration": iteration,
                "dark_acquisition": dark_name,
                "light_acquisition": light_name,
            }
        )
    frame = pd.DataFrame(rows)
    if len(frame) != 19:
        raise ValueError(f"Expected 19 raw acquisition pairs, found {len(frame)}.")
    if frame["well_id"].nunique() != 13:
        raise ValueError(f"Expected exactly 13 unique matched wells, found {frame['well_id'].nunique()}.")
    return frame


def build_all_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    require_inputs()
    compositions = generate_composition_table()
    pl = parse_top_pl_peaks()
    acquisitions = load_kpfm_acquisitions()

    matched_pl = pl.loc[pl["timepoint"] == MATCHED_TIMEPOINT].copy()
    if len(matched_pl) != 96 or matched_pl["elapsed_time_min"].nunique() != 1:
        raise ValueError("Matched PL timepoint does not contain exactly 96 unique wells.")
    actual_time = float(matched_pl["elapsed_time_min"].iloc[0])
    if not math.isclose(actual_time, MATCHED_TIME_MIN):
        raise ValueError(f"Matched timepoint is {actual_time} min, expected {MATCHED_TIME_MIN} min.")

    averaged = acquisitions.groupby("well_id", sort=False).agg(
        SP_dark=("SP_dark", "mean"),
        SP_light=("SP_light", "mean"),
        delta_SP=("delta_SP", "mean"),
        delta_SP_se=("delta_SP_se", lambda s: float(np.sqrt(np.sum(np.square(s))) / len(s))),
        acquisition_iteration=("acquisition_iteration", lambda s: ";".join(map(str, s))),
    ).reset_index()


    averaged["abs_delta_SP"] = averaged["delta_SP"].abs()
    averaged["acquisition_type"] = "Gaussian histogram center; mean across repeated acquisitions"

    matched = (
        averaged.merge(compositions, on="well_id", how="left", validate="one_to_one")
        .merge(
            matched_pl[["well_id", "elapsed_time_min", "matched_PL_intensity"]],
            on="well_id",
            how="left",
            validate="one_to_one",
        )
        .rename(columns={"elapsed_time_min": "matched_PL_time_min"})
    )
    if len(matched) != 13 or matched.isna().any().any() or matched["well_id"].duplicated().any():
        raise ValueError("Matched-table validation failed: missing, duplicate, or non-13 rows.")
    matched = matched[
        [
            "well_id",
            "FAPbI3_pct",
            "BDAPbI4_pct",
            "PEA2PbI4_pct",
            "matched_PL_time_min",
            "matched_PL_intensity",
            "SP_dark",
            "SP_light",
            "delta_SP",
            "delta_SP_se",
            "abs_delta_SP",
            "acquisition_type",
            "acquisition_iteration",
        ]
    ]
    matched.to_csv(OUTDIR / "matched_PL_photoKPFM_processed.csv", index=False)
    acquisitions.to_csv(OUTDIR / "photoKPFM_acquisitions_all_19.csv", index=False)
    compositions.to_csv(OUTDIR / "ternary_96well_composition_lookup_validated.csv", index=False)
    return compositions, pl, matched


def correlation_statistics(matched: pd.DataFrame) -> dict[str, object]:
    x = matched["matched_PL_intensity"].to_numpy(float)
    y = matched["abs_delta_SP"].to_numpy(float)
    pear = pearsonr(x, y)
    confidence = pear.confidence_interval(confidence_level=0.95)
    spear = spearmanr(x, y)
    leave_one_out = []
    for omitted in range(len(matched)):
        keep = np.arange(len(matched)) != omitted
        result = pearsonr(x[keep], y[keep])
        leave_one_out.append(
            {
                "omitted_well_id": matched.iloc[omitted]["well_id"],
                "pearson_r": float(result.statistic),
                "p_value_two_sided": float(result.pvalue),
            }
        )
    report = {
        "n": len(matched),
        "matched_timepoint": MATCHED_TIMEPOINT,
        "matched_time_min": MATCHED_TIME_MIN,
        "pearson_r": float(pear.statistic),
        "pearson_p_value_two_sided": float(pear.pvalue),
        "pearson_r_95_ci": [float(confidence.low), float(confidence.high)],
        "spearman_rho": float(spear.statistic),
        "spearman_p_value_two_sided": float(spear.pvalue),
        "leave_one_out_r_min": float(min(row["pearson_r"] for row in leave_one_out)),
        "leave_one_out_r_max": float(max(row["pearson_r"] for row in leave_one_out)),
        "leave_one_out": leave_one_out,
    }
    (OUTDIR / "matched_PL_photoKPFM_statistics.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    return report


def ternary_xy(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    fa = frame["FAPbI3_pct"].to_numpy(float) / 100.0
    pea = frame["PEA2PbI4_pct"].to_numpy(float) / 100.0
    return pea + 0.5 * fa, (math.sqrt(3) / 2.0) * fa


def draw_manuscript_ternary(ax) -> None:
    """Draw the shared compact manuscript ternary geometry and typography."""
    height = math.sqrt(3) / 2
    ax.plot([0, 1, 0.5, 0], [0, 0, height, 0], color="0.2", lw=1.2)
    for fraction in (0.2, 0.4, 0.6, 0.8):
        ax.plot([0.5 * fraction, 1 - 0.5 * fraction], [height * fraction] * 2, color="0.92", lw=0.5)
        ax.plot([fraction, 0.5 + 0.5 * fraction], [0, height * (1 - fraction)], color="0.92", lw=0.5)
        ax.plot([1 - fraction, 0.5 * (1 - fraction)], [0, height * (1 - fraction)], color="0.92", lw=0.5)
    ax.text(-0.02, -0.055, r"BDAPbI$_4$", ha="center", va="top", fontsize=10)
    ax.text(1.02, -0.055, r"PEA$_2$PbI$_4$", ha="center", va="top", fontsize=10)
    ax.text(0.5, height + 0.045, r"FAPbI$_3$", ha="center", va="bottom", fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.10, height + 0.10)
    ax.axis("off")


draw_ternary_frame = draw_manuscript_ternary


def save_figure(fig, stem: str) -> None:
    for suffix, kwargs in (
        (".png", {"dpi": 600}),
        (".pdf", {}),
        (".svg", {}),
    ):
        fig.savefig(OUTDIR / f"{stem}{suffix}", bbox_inches="tight", facecolor="white", **kwargs)
