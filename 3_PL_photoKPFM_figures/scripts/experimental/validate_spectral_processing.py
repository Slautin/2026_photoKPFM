"""Validation of photoluminescence spectral processing.

This script is intentionally additive. It creates a new validation folder and
does not edit manuscript figures, existing additional-analysis outputs, or raw
data. It treats PL and photoKPFM as nominally matched replicate libraries, not
as spatially colocated or simultaneous measurements.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import math
import shutil
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import spsolve
from scipy.stats import pearsonr, spearmanr


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import analyze_pl_metrics as base  # noqa: E402


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "validation" / "spectral_processing"
FIG = OUT / "figures"
SPEC_FIG = FIG / "spectral_processing_validation"
SUMMARY_FIG = FIG / "summary"

MATCHED_TIMEPOINT = base.MATCHED_TIMEPOINT
MATCHED_TIME_MIN = base.MATCHED_TIME_MIN

REGIONS = {
    "short": (500.0, 620.0),
    "intermediate": (620.0, 720.0),
    "long": (720.0, 850.0),
}
PEAK_SEARCH_RANGE = (500.0, 850.0)
VALID_METRIC_RANGE = (500.0, 850.0)
LATE_STEP = (20, 21)

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 10.5,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
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
class ProcessedSpectrum:
    well_id: str
    geometry: str
    timepoint: int
    read_number: int
    elapsed_time_min: float
    wavelengths: np.ndarray
    raw: np.ndarray
    linear_baseline: np.ndarray
    linear_corrected: np.ndarray
    als_baseline: np.ndarray
    als_corrected: np.ndarray
    flags: list[str]


def ensure_dirs() -> None:
    for path in [OUT, FIG, SPEC_FIG, SUMMARY_FIG]:
        path.mkdir(parents=True, exist_ok=True)


def savefig(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def linear_edge_baseline(wl: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, list[str]]:
    flags: list[str] = []
    finite = np.isfinite(y)
    if finite.sum() < 10:
        return np.full_like(y, np.nan, dtype=float), ["missing_or_too_few_points"]
    edge = (
        ((wl >= base.BASELINE_LOW_RANGE[0]) & (wl <= base.BASELINE_LOW_RANGE[1]))
        | ((wl >= base.BASELINE_HIGH_RANGE[0]) & (wl <= base.BASELINE_HIGH_RANGE[1]))
    ) & finite
    if edge.sum() < 6:
        edge = finite
        flags.append("edge_baseline_sparse")
    coeff = np.polyfit(wl[edge], y[edge], 1)
    baseline = np.polyval(coeff, wl)
    return baseline.astype(float), flags


def als_baseline(y: np.ndarray, lam: float = 1e6, p: float = 0.002, niter: int = 10) -> np.ndarray:
    finite = np.isfinite(y)
    if finite.sum() < 10:
        return np.full_like(y, np.nan, dtype=float)
    filled = np.asarray(y, dtype=float).copy()
    if not finite.all():
        x = np.arange(len(y))
        filled[~finite] = np.interp(x[~finite], x[finite], filled[finite])
    length = len(filled)
    dmat = csc_matrix(np.diff(np.eye(length), 2, axis=0))
    weights = np.ones(length)
    for _ in range(niter):
        wmat = diags(weights, 0, shape=(length, length))
        z = spsolve(wmat + lam * dmat.T @ dmat, weights * filled)
        weights = p * (filled > z) + (1 - p) * (filled < z)
    return np.asarray(z, dtype=float)


def smooth_for_peaks(wl: np.ndarray, y: np.ndarray) -> np.ndarray:
    finite = np.isfinite(y)
    if finite.sum() < 7:
        return y.copy()
    yy = y.copy()
    if not finite.all():
        x = np.arange(len(y))
        yy[~finite] = np.interp(x[~finite], x[finite], yy[finite])
    step = float(np.nanmedian(np.diff(wl)))
    window = max(7, int(round(base.SMOOTH_WINDOW_NM / max(step, 1e-9))))
    window = base.odd_window(len(yy), window)
    return savgol_filter(yy, window, 2)


def half_max_width(wl: np.ndarray, y: np.ndarray, peak_idx: int) -> tuple[float, float, float, bool]:
    ymax = float(y[peak_idx])
    if not np.isfinite(ymax) or ymax <= 0:
        return np.nan, np.nan, np.nan, False
    half = ymax / 2
    left_candidates = np.where(y[: peak_idx + 1] <= half)[0]
    right_candidates = np.where(y[peak_idx:] <= half)[0]
    if len(left_candidates) == 0 or len(right_candidates) == 0:
        return np.nan, np.nan, np.nan, False
    li = int(left_candidates[-1])
    ri = int(peak_idx + right_candidates[0])
    if li >= peak_idx or ri <= peak_idx:
        return np.nan, np.nan, np.nan, False

    def interp_cross(i0: int, i1: int) -> float:
        x0, x1 = float(wl[i0]), float(wl[i1])
        y0, y1 = float(y[i0]), float(y[i1])
        if abs(y1 - y0) < 1e-12:
            return x0
        return x0 + (half - y0) * (x1 - x0) / (y1 - y0)

    left = interp_cross(li, li + 1)
    right = interp_cross(ri - 1, ri)
    return right - left, float(wl[peak_idx] - left), float(right - wl[peak_idx]), True


def detect_and_classify(wl: np.ndarray, corrected: np.ndarray) -> dict[str, object]:
    mask = (wl >= PEAK_SEARCH_RANGE[0]) & (wl <= PEAK_SEARCH_RANGE[1]) & np.isfinite(corrected)
    if mask.sum() < 10:
        return {
            "spectrum_class": "too_weak_or_missing",
            "flags": ["missing_or_too_few_points"],
            "peak_indices": [],
            "dominant_idx": None,
            "fwhm_nm": np.nan,
            "fwhm_defensible": False,
            "dominant_peak_wavelength_nm": np.nan,
            "dominant_peak_intensity": np.nan,
            "secondary_peak_wavelength_nm": np.nan,
            "secondary_peak_intensity": np.nan,
            "num_peaks": 0,
            "asymmetry_factor": np.nan,
        }

    x = wl[mask]
    y = corrected[mask]
    y_s = smooth_for_peaks(x, y)
    max_y = float(np.nanmax(y_s))
    flags: list[str] = []
    if max_y <= 0 or not np.isfinite(max_y):
        flags.append("nonpositive_after_baseline")
    if max_y < 100:
        flags.append("low_signal")

    step = float(np.nanmedian(np.diff(x)))
    distance = max(1, int(round(18.0 / max(step, 1e-9))))
    prominence = max(0.03 * max_y, 50.0)
    peak_rel, props = find_peaks(y_s, prominence=prominence, distance=distance)
    if len(peak_rel) == 0 and np.isfinite(max_y) and max_y > 0:
        peak_rel = np.array([int(np.nanargmax(y_s))])
        props = {"prominences": np.array([max_y])}
    if len(peak_rel) == 0:
        flags.append("no_detected_peak")
        return {
            "spectrum_class": "too_weak_or_missing",
            "flags": flags,
            "peak_indices": [],
            "dominant_idx": None,
            "fwhm_nm": np.nan,
            "fwhm_defensible": False,
            "dominant_peak_wavelength_nm": np.nan,
            "dominant_peak_intensity": np.nan,
            "secondary_peak_wavelength_nm": np.nan,
            "secondary_peak_intensity": np.nan,
            "num_peaks": 0,
            "asymmetry_factor": np.nan,
        }

    order = np.argsort(y_s[peak_rel])[::-1]
    peak_rel = peak_rel[order]
    peak_abs = np.where(mask)[0][peak_rel]
    dominant_rel = int(peak_rel[0])
    dom_wl = float(x[dominant_rel])
    dom_int = float(y_s[dominant_rel])
    near_boundary = dom_wl < PEAK_SEARCH_RANGE[0] + 20 or dom_wl > PEAK_SEARCH_RANGE[1] - 20
    if near_boundary:
        flags.append("dominant_peak_near_boundary")

    secondary_wl = np.nan
    secondary_int = np.nan
    strong_secondary = False
    shoulder = False
    if len(peak_rel) > 1:
        secondary_rel = int(peak_rel[1])
        secondary_wl = float(x[secondary_rel])
        secondary_int = float(y_s[secondary_rel])
        sep = abs(secondary_wl - dom_wl)
        rel_height = secondary_int / max(dom_int, 1e-12)
        strong_secondary = sep >= 30 and rel_height >= 0.25
        shoulder = rel_height >= 0.12

    fwhm, left_hw, right_hw, crossings = half_max_width(x, y_s, dominant_rel)
    asym = right_hw / left_hw if crossings and left_hw > 0 else np.nan
    if crossings and np.isfinite(asym) and (asym > 2.0 or asym < 0.5):
        flags.append("strongly_asymmetric")
        shoulder = True

    if max_y < 100:
        spectrum_class = "too_weak_or_missing"
    elif near_boundary:
        spectrum_class = "boundary_peak"
    elif strong_secondary:
        spectrum_class = "multi_peaked"
    elif shoulder:
        spectrum_class = "strongly_shouldered_or_asymmetric"
    else:
        spectrum_class = "single_peaked"

    fwhm_defensible = bool(spectrum_class == "single_peaked" and crossings and np.isfinite(fwhm))
    if not fwhm_defensible:
        flags.append("fwhm_not_defensible")

    return {
        "spectrum_class": spectrum_class,
        "flags": flags,
        "peak_indices": peak_abs.tolist(),
        "dominant_idx": int(peak_abs[0]),
        "fwhm_nm": float(fwhm) if fwhm_defensible else np.nan,
        "fwhm_defensible": fwhm_defensible,
        "dominant_peak_wavelength_nm": dom_wl,
        "dominant_peak_intensity": dom_int,
        "secondary_peak_wavelength_nm": secondary_wl,
        "secondary_peak_intensity": secondary_int,
        "num_peaks": int(len(peak_rel)),
        "asymmetry_factor": float(asym) if np.isfinite(asym) else np.nan,
    }


def process_spectrum(block: base.SpectrumBlock, well_id: str) -> ProcessedSpectrum:
    idx = base.WELLS_ROW_MAJOR.index(well_id)
    wl = block.wavelengths.astype(float)
    raw = block.values[:, idx].astype(float)
    flags: list[str] = []
    if np.isnan(raw).any():
        flags.append("missing_or_overflow_points")
    if np.nanmax(raw) >= 65000:
        flags.append("possible_saturation_or_clipping")
    lin_base, lin_flags = linear_edge_baseline(wl, raw)
    lin_corr = raw - lin_base
    lin_corr[lin_corr < 0] = 0
    als_base = np.full_like(raw, np.nan, dtype=float)
    als_corr = np.full_like(raw, np.nan, dtype=float)
    return ProcessedSpectrum(
        well_id=well_id,
        geometry=block.geometry,
        timepoint=block.timepoint,
        read_number=block.read_number,
        elapsed_time_min=block.elapsed_time_min,
        wavelengths=wl,
        raw=raw,
        linear_baseline=lin_base,
        linear_corrected=lin_corr,
        als_baseline=als_base,
        als_corrected=als_corr,
        flags=flags + lin_flags,
    )


def band_descriptors(
    ps: ProcessedSpectrum,
    comp: pd.Series,
    method: str = "linear_edge",
    als_corrected: np.ndarray | None = None,
) -> dict[str, object]:
    corrected = ps.linear_corrected if method == "linear_edge" else als_corrected
    if corrected is None:
        corrected = ps.linear_corrected
    det = detect_and_classify(ps.wavelengths, corrected)
    valid = (ps.wavelengths >= VALID_METRIC_RANGE[0]) & (ps.wavelengths <= VALID_METRIC_RANGE[1]) & np.isfinite(corrected)
    total_area = float(np.trapezoid(corrected[valid], ps.wavelengths[valid])) if valid.sum() > 1 else np.nan
    row: dict[str, object] = {
        "well_id": ps.well_id,
        "read_geometry": ps.geometry,
        "acquisition_number": ps.timepoint,
        "read_number": ps.read_number,
        "elapsed_time_min": ps.elapsed_time_min,
        "baseline_method": method,
        "FAPbI3_pct": comp["FAPbI3_pct"],
        "BDAPbI4_pct": comp["BDAPbI4_pct"],
        "PEA2PbI4_pct": comp["PEA2PbI4_pct"],
        "total_area_500_850": total_area,
        "centroid_nm": np.nan,
        "dominant_peak_wavelength_nm": det["dominant_peak_wavelength_nm"],
        "dominant_peak_intensity": det["dominant_peak_intensity"],
        "secondary_peak_wavelength_nm": det["secondary_peak_wavelength_nm"],
        "secondary_peak_intensity": det["secondary_peak_intensity"],
        "num_peaks_or_shoulders": det["num_peaks"],
        "spectrum_class": det["spectrum_class"],
        "fwhm_nm": det["fwhm_nm"],
        "fwhm_defensible": det["fwhm_defensible"],
        "asymmetry_factor": det["asymmetry_factor"],
        "qc_flags": ";".join(sorted(set(ps.flags + list(det["flags"])))),
    }
    if valid.sum() > 1 and total_area > 0:
        row["centroid_nm"] = float(np.trapezoid(ps.wavelengths[valid] * corrected[valid], ps.wavelengths[valid]) / total_area)
    for name, (lo, hi) in REGIONS.items():
        mask = (ps.wavelengths >= lo) & (ps.wavelengths < hi) & np.isfinite(corrected)
        area = float(np.trapezoid(corrected[mask], ps.wavelengths[mask])) if mask.sum() > 1 else np.nan
        row[f"{name}_area"] = area
        row[f"{name}_fraction"] = area / total_area if np.isfinite(area) and np.isfinite(total_area) and total_area > 0 else np.nan
    row["long_short_area_ratio"] = (
        row["long_area"] / row["short_area"]
        if np.isfinite(row["long_area"]) and np.isfinite(row["short_area"]) and row["short_area"] > 0
        else np.nan
    )
    return row


def plot_spectral_validation(processed: list[ProcessedSpectrum], comp_lookup: pd.DataFrame, label: str) -> None:
    if not processed:
        return
    for ps in processed:
        comp = comp_lookup.loc[ps.well_id]
        title = (
            f"{ps.well_id} {ps.geometry}\n"
            f"FA {comp['FAPbI3_pct']:.0f} / BDA {comp['BDAPbI4_pct']:.0f} / PEA {comp['PEA2PbI4_pct']:.0f}%"
        )
        fig, ax = plt.subplots(figsize=(5.8, 4.2))
        ax.plot(ps.wavelengths, ps.raw, color="0.15", lw=1.0, label="raw")
        ax.plot(ps.wavelengths, ps.linear_baseline, color="#d95f02", lw=1.0, label="linear baseline")
        ax.plot(ps.wavelengths, ps.als_baseline, color="#1b9e77", lw=1.0, label="ALS baseline")
        ax.set_title(f"{title}\nRaw spectrum and baselines", fontweight="bold")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Raw PL")
        ax.legend(frameon=False, loc="upper right", fontsize=7)
        savefig(
            fig,
            SPEC_FIG / f"spectral_processing_validation_{label}_{ps.geometry}_raw_baselines",
        )

        fig, ax = plt.subplots(figsize=(5.8, 4.2))
        ax.plot(ps.wavelengths, ps.linear_corrected, color="#4c78a8", lw=1.1, label="baseline-subtracted")
        det = detect_and_classify(ps.wavelengths, ps.linear_corrected)
        for lo, hi in REGIONS.values():
            ax.axvspan(lo, hi, color="0.9", alpha=0.35, zorder=0)
        for pi in det["peak_indices"]:
            if pi is not None:
                ax.plot(ps.wavelengths[pi], ps.linear_corrected[pi], "o", ms=4, color="#f58518")
        if det["dominant_idx"] is not None:
            di = int(det["dominant_idx"])
            ax.plot(ps.wavelengths[di], ps.linear_corrected[di], "*", ms=9, color="#e45756")
        ax.text(
            0.02,
            0.95,
            f"{det['spectrum_class']}\nflags: {det['flags'][:3]}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            bbox=dict(fc="white", ec="0.75", alpha=0.9, boxstyle="round,pad=0.2"),
        )
        ax.set_title(f"{title}\nBaseline-corrected spectrum", fontweight="bold")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Corrected PL")
        ax.legend(frameon=False, loc="upper right", fontsize=7)
        savefig(
            fig,
            SPEC_FIG / f"spectral_processing_validation_{label}_{ps.geometry}_corrected",
        )


def draw_ternary(ax) -> None:
    def xy_pct(fa: float, bda: float, pea: float) -> tuple[float, float]:
        fa_f = fa / 100.0
        pea_f = pea / 100.0
        return pea_f + 0.5 * fa_f, (math.sqrt(3) / 2.0) * fa_f

    vertices = np.array([[0, 0], [1, 0], [0.5, math.sqrt(3) / 2], [0, 0]])
    ax.plot(vertices[:, 0], vertices[:, 1], color="black", lw=1.4)
    for f in np.linspace(0.2, 0.8, 4):
        x1, y1 = xy_pct(100 * f, 100 * (1 - f), 0)
        x2, y2 = xy_pct(100 * f, 0, 100 * (1 - f))
        ax.plot([x1, x2], [y1, y2], color="0.88", lw=0.6, zorder=0)
        x1, y1 = xy_pct(100 * (1 - f), 100 * f, 0)
        x2, y2 = xy_pct(0, 100 * f, 100 * (1 - f))
        ax.plot([x1, x2], [y1, y2], color="0.88", lw=0.6, zorder=0)
        x1, y1 = xy_pct(100 * (1 - f), 0, 100 * f)
        x2, y2 = xy_pct(0, 100 * (1 - f), 100 * f)
        ax.plot([x1, x2], [y1, y2], color="0.88", lw=0.6, zorder=0)
    ax.text(0.50, math.sqrt(3) / 2 + 0.045, "FAPbI3", ha="center", va="bottom", fontsize=13, fontweight="bold")
    ax.text(-0.055, -0.055, "BDAPbI4", ha="right", va="top", fontsize=13, fontweight="bold")
    ax.text(1.055, -0.055, "PEA2PbI4", ha="left", va="top", fontsize=13, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")


def ternary_metric_map(frame: pd.DataFrame, value_col: str, stem: str, label: str, cmap: str = "viridis") -> None:
    data = frame.dropna(subset=[value_col]).copy()
    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    draw_ternary(ax)
    x, y = base.ternary_xy(data)
    xy = np.column_stack([x, y])
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=data[value_col], cmap=cmap, s=54, edgecolor="black", linewidth=0.45)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label)
    savefig(fig, SUMMARY_FIG / stem)


def correlation_summary(x: pd.Series, y: pd.Series) -> dict[str, float]:
    mask = np.isfinite(x.to_numpy(float)) & np.isfinite(y.to_numpy(float))
    n = int(mask.sum())
    out = {"n": n, "pearson_r": np.nan, "pearson_p": np.nan, "spearman_rho": np.nan, "spearman_p": np.nan}
    if n < 3:
        return out
    xv = x.to_numpy(float)[mask]
    yv = y.to_numpy(float)[mask]
    if np.nanstd(xv) <= 0 or np.nanstd(yv) <= 0:
        return out
    pr = pearsonr(xv, yv)
    sr = spearmanr(xv, yv)
    out.update({"pearson_r": float(pr.statistic), "pearson_p": float(pr.pvalue), "spearman_rho": float(sr.statistic), "spearman_p": float(sr.pvalue)})
    return out


def leave_one_out_pearson(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    mask = np.isfinite(x.to_numpy(float)) & np.isfinite(y.to_numpy(float))
    xv = x.to_numpy(float)[mask]
    yv = y.to_numpy(float)[mask]
    vals = []
    if len(xv) < 4:
        return np.nan, np.nan
    for i in range(len(xv)):
        xx = np.delete(xv, i)
        yy = np.delete(yv, i)
        if np.nanstd(xx) > 0 and np.nanstd(yy) > 0:
            vals.append(float(pearsonr(xx, yy).statistic))
    return (float(np.nanmin(vals)), float(np.nanmax(vals))) if vals else (np.nan, np.nan)


def build_processed_table(
    blocks: list[base.SpectrumBlock],
    compositions: pd.DataFrame,
    matched_wells: list[str],
) -> tuple[pd.DataFrame, dict[tuple[str, str, int], ProcessedSpectrum]]:
    comp_lookup = compositions.set_index("well_id")
    rows = []
    spectra: dict[tuple[str, str, int], ProcessedSpectrum] = {}
    for block in blocks:
        if block.geometry not in {"top", "bottom"}:
            continue
        for well in base.WELLS_ROW_MAJOR:
            ps = process_spectrum(block, well)
            spectra[(well, block.geometry, block.timepoint)] = ps
            comp = comp_lookup.loc[well]
            rows.append(band_descriptors(ps, comp, "linear_edge"))


            if block.timepoint == MATCHED_TIMEPOINT and well in matched_wells:
                als_base = als_baseline(ps.raw)
                als_corr = ps.raw - als_base
                als_corr[als_corr < 0] = 0
                ps.als_baseline = als_base
                ps.als_corrected = als_corr
                rows.append(band_descriptors(ps, comp, "als", als_corrected=als_corr))
    return pd.DataFrame(rows), spectra


def representative_wells(compositions: pd.DataFrame, matched_wells: list[str]) -> dict[str, str]:
    comp = compositions.set_index("well_id")
    reps = {
        "PEA_rich": comp["PEA2PbI4_pct"].idxmax(),
        "BDA_rich": comp["BDAPbI4_pct"].idxmax(),
        "FA_rich": comp["FAPbI3_pct"].idxmax(),
    }
    target = np.array([34.0, 33.0, 33.0])
    vals = comp[["FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]].to_numpy(float)
    reps["intermediate"] = comp.index[int(np.nanargmin(np.linalg.norm(vals - target, axis=1)))]
    if matched_wells:
        reps["matched_high_delta_SP"] = matched_wells[0]
    return reps


def make_spectral_validation_figures(spectra: dict[tuple[str, str, int], ProcessedSpectrum], compositions: pd.DataFrame, matched_wells: list[str]) -> None:
    comp_lookup = compositions.set_index("well_id")
    for well in matched_wells:
        procs = [spectra[(well, geom, MATCHED_TIMEPOINT)] for geom in ["top", "bottom"] if (well, geom, MATCHED_TIMEPOINT) in spectra]
        plot_spectral_validation(procs, comp_lookup, f"matched_{well}_timepoint_{MATCHED_TIMEPOINT}")
    reps = representative_wells(compositions, matched_wells)
    for name, well in reps.items():
        procs = [spectra[(well, geom, MATCHED_TIMEPOINT)] for geom in ["top", "bottom"] if (well, geom, MATCHED_TIMEPOINT) in spectra]
        plot_spectral_validation(procs, comp_lookup, f"representative_{name}_{well}_timepoint_{MATCHED_TIMEPOINT}")


def baseline_validation(desc: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    top22 = desc[(desc.read_geometry.eq("top")) & (desc.acquisition_number.eq(MATCHED_TIMEPOINT))]
    lin = top22[top22.baseline_method.eq("linear_edge")].set_index("well_id")
    als = top22[top22.baseline_method.eq("als")].set_index("well_id")
    common = lin.index.intersection(als.index)
    rows = []
    for col in ["total_area_500_850", "short_area", "intermediate_area", "long_area", "long_fraction", "dominant_peak_intensity"]:
        merged = pd.DataFrame({"linear_edge": lin.loc[common, col], "als": als.loc[common, col]})
        merged["abs_rank_shift"] = merged["linear_edge"].rank(ascending=False) - merged["als"].rank(ascending=False)
        rows.append(
            {
                "descriptor": col,
                "n": int(merged.dropna().shape[0]),
                "median_relative_difference_als_vs_linear": float(np.nanmedian((merged["als"] - merged["linear_edge"]) / merged["linear_edge"].replace(0, np.nan))),
                "max_abs_rank_shift": float(np.nanmax(np.abs(merged["abs_rank_shift"]))),
                "baseline_dependent_wells_rank_shift_ge_10": int((np.abs(merged["abs_rank_shift"]) >= 10).sum()),
            }
        )
    comp = lin.reset_index().merge(als.reset_index()[["well_id", "long_fraction"]], on="well_id", suffixes=("_linear", "_als"))
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    ax.scatter(comp["long_fraction_linear"], comp["long_fraction_als"], s=36, edgecolor="black", color="#4c78a8")
    lim = [0, np.nanmax(comp[["long_fraction_linear", "long_fraction_als"]].to_numpy(float)) * 1.05]
    ax.plot(lim, lim, color="0.3", lw=1, ls="--")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Long fraction, linear-edge baseline")
    ax.set_ylabel("Long fraction, ALS baseline")
    savefig(fig, SUMMARY_FIG / "baseline_validation_long_fraction_linear_vs_als")
    return pd.DataFrame(rows)


def photokpfm_over_time(desc: pd.DataFrame, matched: pd.DataFrame, pl_existing: pd.DataFrame) -> pd.DataFrame:
    matched_small = matched[["well_id", "abs_delta_SP"]].copy()
    rows = []
    top = desc[(desc.read_geometry.eq("top")) & (desc.baseline_method.eq("linear_edge"))].copy()
    top = top.merge(matched_small, on="well_id", how="inner")
    existing = pl_existing.rename(columns={"timepoint": "acquisition_number"})
    if "matched_PL_intensity" in existing.columns:
        top = top.merge(existing[["well_id", "acquisition_number", "matched_PL_intensity"]], on=["well_id", "acquisition_number"], how="left")
    descriptors = [
        "matched_PL_intensity",
        "dominant_peak_intensity",
        "total_area_500_850",
        "short_area",
        "intermediate_area",
        "long_area",
        "short_fraction",
        "intermediate_fraction",
        "long_fraction",
        "long_short_area_ratio",
        "dominant_peak_wavelength_nm",
    ]
    for tp, group in top.groupby("acquisition_number"):
        for descriptor in descriptors:
            if descriptor not in group.columns:
                continue
            stats = correlation_summary(group[descriptor], group["abs_delta_SP"])
            lo, hi = leave_one_out_pearson(group[descriptor], group["abs_delta_SP"])
            x = group[descriptor].to_numpy(float)
            y = group["abs_delta_SP"].to_numpy(float)
            valid = np.isfinite(x) & np.isfinite(y)
            rank_disagree = np.nan
            if valid.sum() >= 3:
                rank_disagree = float(np.nanmean(np.abs(pd.Series(x[valid]).rank(ascending=False).to_numpy() - pd.Series(y[valid]).rank(ascending=False).to_numpy())))
            rows.append(
                {
                    "acquisition_number": int(tp),
                    "elapsed_time_min": float((int(tp) - 1) * 9),
                    "descriptor": descriptor,
                    **stats,
                    "leave_one_out_pearson_min": lo,
                    "leave_one_out_pearson_max": hi,
                    "mean_abs_rank_disagreement": rank_disagree,
                }
            )
    stats_df = pd.DataFrame(rows)
    for stat_col, out_name, ylabel in [
        ("pearson_r", "pl_photokpfm_pearson_r_vs_time", "Pearson r with |Delta SP|"),
        ("spearman_rho", "pl_photokpfm_spearman_rho_vs_time", "Spearman rho with |Delta SP|"),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        show = stats_df[stats_df["descriptor"].isin(["matched_PL_intensity", "long_fraction", "long_area", "dominant_peak_wavelength_nm", "total_area_500_850"])]
        for descriptor, group in show.groupby("descriptor"):
            ax.plot(group["elapsed_time_min"], group[stat_col], marker="o", lw=1.2, label=descriptor.replace("_", " "))
        for tp in [1, 20, 21, 22, int(stats_df["acquisition_number"].max())]:
            ax.axvline((tp - 1) * 9, color="0.85", lw=0.8, zorder=0)
        ax.axhline(0, color="0.25", lw=0.9)
        ax.set_xlabel("Elapsed PL time by existing convention (min)")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=7, ncol=2)
        savefig(fig, SUMMARY_FIG / out_name)
    return stats_df


def wavelength_region_summary(desc: pd.DataFrame) -> pd.DataFrame:
    top22 = desc[(desc.baseline_method.eq("linear_edge")) & (desc.acquisition_number.eq(MATCHED_TIMEPOINT))]
    rows = []
    for geom, group in top22.groupby("read_geometry"):
        for region, (lo, hi) in REGIONS.items():
            peak_in_region = group["dominant_peak_wavelength_nm"].between(lo, hi)
            rows.append(
                {
                    "read_geometry": geom,
                    "region": region,
                    "range_nm": f"{lo:.0f}-{hi:.0f}",
                    "n_dominant_peaks_in_region": int(peak_in_region.sum()),
                    "median_fraction": float(np.nanmedian(group[f"{region}_fraction"])),
                    "mean_fraction": float(np.nanmean(group[f"{region}_fraction"])),
                }
            )
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    vals = top22["dominant_peak_wavelength_nm"].dropna()
    ax.hist(vals, bins=np.arange(500, 855, 10), color="#4c78a8", edgecolor="white")
    for name, (lo, hi) in REGIONS.items():
        ax.axvspan(lo, hi, alpha=0.12, label=f"{name}: {lo:.0f}-{hi:.0f} nm")
    ax.set_xlabel("Detected dominant peak wavelength (nm)")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, fontsize=8)
    savefig(fig, SUMMARY_FIG / "dominant_peak_wavelength_distribution_regions")
    return pd.DataFrame(rows)


def late_time_loss_analysis(desc: pd.DataFrame) -> pd.DataFrame:
    lin = desc[desc.baseline_method.eq("linear_edge")]
    rows = []
    for geom in ["top", "bottom"]:
        before = lin[(lin.read_geometry.eq(geom)) & (lin.acquisition_number.eq(LATE_STEP[0]))].set_index("well_id")
        after = lin[(lin.read_geometry.eq(geom)) & (lin.acquisition_number.eq(LATE_STEP[1]))].set_index("well_id")
        common = before.index.intersection(after.index)
        for metric in ["total_area_500_850", "short_area", "intermediate_area", "long_area", "dominant_peak_wavelength_nm"]:
            b = before.loc[common, metric].astype(float)
            a = after.loc[common, metric].astype(float)
            change = (a - b) / b.replace(0, np.nan) if metric != "dominant_peak_wavelength_nm" else a - b
            rows.append(
                {
                    "read_geometry": geom,
                    "from_acquisition": LATE_STEP[0],
                    "to_acquisition": LATE_STEP[1],
                    "metric": metric,
                    "median_change": float(np.nanmedian(change)),
                    "mean_change": float(np.nanmean(change)),
                    "fraction_wells_decreased_gt_25pct": float(np.nanmean(change < -0.25)) if metric != "dominant_peak_wavelength_nm" else np.nan,
                }
            )
    frame = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    area_rows = frame[frame["metric"].isin(["short_area", "intermediate_area", "long_area"])]
    for geom, group in area_rows.groupby("read_geometry"):
        ax.plot(group["metric"].str.replace("_area", ""), group["median_change"], marker="o", label=geom)
    ax.axhline(0, color="0.25", lw=0.9)
    ax.set_ylabel(f"Median fractional change, acquisition {LATE_STEP[0]} to {LATE_STEP[1]}")
    ax.set_xlabel("Wavelength region")
    ax.legend(frameon=False)
    savefig(fig, SUMMARY_FIG / "late_time_loss_band_resolved_summary")
    return frame


def top_bottom_summary(desc: pd.DataFrame, matched_wells: list[str]) -> pd.DataFrame:
    lin22 = desc[
        (desc.baseline_method.eq("linear_edge"))
        & (desc.acquisition_number.eq(MATCHED_TIMEPOINT))
        & (desc.well_id.isin(matched_wells))
    ]
    top = lin22[lin22.read_geometry.eq("top")].set_index("well_id")
    bottom = lin22[lin22.read_geometry.eq("bottom")].set_index("well_id")
    common = top.index.intersection(bottom.index)
    rows = []
    for metric in ["dominant_peak_wavelength_nm", "total_area_500_850", "short_fraction", "intermediate_fraction", "long_fraction"]:
        diff = bottom.loc[common, metric].astype(float) - top.loc[common, metric].astype(float)
        rows.append(
            {
                "metric": metric,
                "n": int(diff.notna().sum()),
                "median_bottom_minus_top": float(np.nanmedian(diff)),
                "mean_bottom_minus_top": float(np.nanmean(diff)),
            }
        )
    return pd.DataFrame(rows)


def write_status_reports(
    desc: pd.DataFrame,
    baseline: pd.DataFrame,
    corr: pd.DataFrame,
    regions: pd.DataFrame,
    late: pd.DataFrame,
    top_bottom: pd.DataFrame,
    matched_wells: list[str],
) -> None:
    lin22_matched = desc[
        (desc.baseline_method.eq("linear_edge"))
        & (desc.acquisition_number.eq(MATCHED_TIMEPOINT))
        & (desc.well_id.isin(matched_wells))
    ]
    fwhm = lin22_matched.groupby("read_geometry")["fwhm_defensible"].agg(["sum", "count"]).reset_index()
    class_counts = (
        lin22_matched.groupby(["read_geometry", "spectrum_class"]).size().reset_index(name="n").sort_values(["read_geometry", "n"], ascending=[True, False])
    )
    best_corr = corr.dropna(subset=["pearson_r"]).sort_values("pearson_r", ascending=False).head(10)
    best_corr_cols = [
        "acquisition_number",
        "elapsed_time_min",
        "descriptor",
        "n",
        "pearson_r",
        "pearson_p",
        "spearman_rho",
        "spearman_p",
        "leave_one_out_pearson_min",
        "leave_one_out_pearson_max",
        "mean_abs_rank_disagreement",
    ]
    summary = f"""# Additional PL Validation Summary

This folder is a validation package only. It does not replace manuscript figures or the earlier peak-intensity analysis.

## Critical Scope

- PL time-series data and photoKPFM data are treated as separately prepared replicate libraries with nominally matched compositions.
- PL and photoKPFM are not described as simultaneous, spatially colocated, or measured on the same physical droplets.
- PL time evolution is interpreted only for the PL replicate during repeated Cytation acquisition.
- FWHM and full-spectrum area are not used as replacement primary metrics because many spectra are multi-featured, shouldered, asymmetric, or baseline-sensitive.

## Wavelength Regions Used

Neutral optical regions were used for descriptors, not phase assignments:

{base.markdown_table(regions)}

## FWHM Defensibility at Timepoint {MATCHED_TIMEPOINT}

FWHM is retained only for isolated single peaks with valid half-max crossings.

{base.markdown_table(fwhm)}

Spectrum-class counts for matched wells:

{base.markdown_table(class_counts)}

## Baseline Validation

Linear-edge baseline and ALS baseline were compared. Descriptors with large rank shifts should be treated as baseline-dependent.

{base.markdown_table(baseline)}

## PL-photoKPFM Correlations Over PL Acquisition

These statistics compare descriptors from the PL replicate at each top-read acquisition with |Delta SP| from the separate photoKPFM replicate.
They are validation diagnostics, not evidence of colocated time evolution.

Top Pearson-r entries:

{base.markdown_table(best_corr[best_corr_cols])}

## Late-Time PL Loss

Use the neutral wording:

> A broad late-time loss of measured PL occurred in the separately prepared PL replicate during repeated optical acquisition.

Band-resolved summary for acquisition {LATE_STEP[0]} to {LATE_STEP[1]}:

{base.markdown_table(late)}

## Top versus Bottom Reads

These are read-geometry-dependent optical differences only, not proof of vertical phase segregation.

{base.markdown_table(top_bottom)}
"""
    (OUT / "validation_analysis_summary.md").write_text(summary, encoding="utf-8")

    status_rows = [
        ("Yang-style analysis", "Added band-resolved optical descriptors, spectra QC, top/bottom checks, and composition maps; no structural phase assignments from PL."),
        ("Timepoint 22", f"Verified as the 22nd top-read PL acquisition, existing convention approximately {MATCHED_TIME_MIN} min after first top read."),
        ("189 min", "Retained as the existing elapsed-time convention; exact per-read timestamps are not available in exported PL blocks."),
        ("G1/B11 labels", "Matched-well spectral validation plots use well IDs and nominal compositions directly from the composition table."),
        ("Late-time PL loss", "Analyzed spectrally as broad late-time measured PL loss in the separate PL replicate; no claim that photoKPFM samples degraded."),
        ("Whether KPFM samples had degraded", "Not inferable from PL replicate behavior; explicitly excluded."),
        ("Peak intensity vs area", "Compared peak intensity, total area, and band-resolved area/fraction without replacing original peak-intensity metric."),
        ("Peak-position variation", "Dominant and secondary peak wavelengths tabulated and plotted by wavelength region."),
        ("Multiple components", "Classified spectra as single, multi, shouldered/asymmetric, boundary, or too weak; FWHM withheld when not defensible."),
        ("Top/bottom", "Compared as read-geometry-dependent optical differences only."),
        ("n-phase assignment", "Not assigned from PL; regions are called short/intermediate/long wavelength emission."),
        ("Spacer passivation", "Not directly proven by PL descriptor statistics; can be discussed as hypothesis with literature/structural support."),
        ("Single-row figure layout", "No manuscript replacement figures generated in this validation pass."),
        ("Linear fit", "Correlation statistics by acquisition include Pearson/Spearman/leave-one-out; no single best descriptor forced."),
        ("Proposed PEA-stability interpretation", "Supported cautiously as composition-dependent optical behavior; structural mechanism remains outside PL-only validation."),
    ]
    status = "# Spectral Validation Status\n\n" + base.markdown_table(pd.DataFrame(status_rows, columns=["comment_topic", "validation_status"]))
    (OUT / "spectral_validation_status.md").write_text(status, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    compositions, pl_existing, raw_matched = base.build_all_data()
    matched = base.prepare_data(pl_existing, raw_matched)
    matched_wells = matched["well_id"].tolist()
    blocks, start_dt, skipped_empty = base.parse_fluorescence_blocks()

    desc, spectra = build_processed_table(blocks, compositions, matched_wells)
    desc.to_csv(OUT / "band_resolved_pl_descriptors_all_reads.csv", index=False)

    make_spectral_validation_figures(spectra, compositions, matched_wells)

    lin22_top = desc[(desc.baseline_method.eq("linear_edge")) & (desc.read_geometry.eq("top")) & (desc.acquisition_number.eq(MATCHED_TIMEPOINT))]
    for value_col, label, cmap in [
        ("short_fraction", "Short-wavelength fraction", "Blues"),
        ("intermediate_fraction", "Intermediate-wavelength fraction", "viridis"),
        ("long_fraction", "Long-wavelength fraction", "magma"),
        ("dominant_peak_wavelength_nm", "Dominant peak wavelength (nm)", "plasma"),
    ]:
        ternary_metric_map(lin22_top, value_col, f"top_timepoint_{MATCHED_TIMEPOINT}_{value_col}_ternary", label, cmap)

    baseline = baseline_validation(desc, matched)
    baseline.to_csv(OUT / "baseline_method_comparison.csv", index=False)

    regions = wavelength_region_summary(desc)
    regions.to_csv(OUT / "wavelength_region_summary.csv", index=False)

    corr = photokpfm_over_time(desc, matched, pl_existing)
    corr.to_csv(OUT / "pl_photokpfm_by_time_statistics.csv", index=False)

    late = late_time_loss_analysis(desc)
    late.to_csv(OUT / "late_time_spectral_loss_analysis.csv", index=False)

    top_bottom = top_bottom_summary(desc, matched_wells)
    top_bottom.to_csv(OUT / "top_bottom_band_resolved_comparison.csv", index=False)

    write_status_reports(desc, baseline, corr, regions, late, top_bottom, matched_wells)
    reproducibility = {
        "raw_pl_csv": str(base.RAW_PL_CSV),
        "project_root": str(ROOT),
        "output_folder": str(OUT),
        "matched_timepoint": MATCHED_TIMEPOINT,
        "matched_time_min_existing_convention": MATCHED_TIME_MIN,
        "pl_photokpfm_relationship_scope": "nominally matched replicate libraries, not same physical droplets",
        "wavelength_regions_nm": REGIONS,
        "baseline_methods": ["linear_edge", "als"],
        "late_loss_step_checked": LATE_STEP,
        "skipped_fully_empty_pl_blocks": skipped_empty,
        "start_datetime_from_file_header": str(start_dt) if start_dt is not None else "",
    }
    (OUT / "validation_reproducibility.json").write_text(json.dumps(reproducibility, indent=2), encoding="utf-8")
    print(f"Saved validation package to: {OUT}")


if __name__ == "__main__":
    main()
