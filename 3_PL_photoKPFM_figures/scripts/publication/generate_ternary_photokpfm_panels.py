"""Generate standalone ternary photoKPFM manuscript panels from the bundled HDF5 file.

The historical analysis notebook used a BoTorch/GPyTorch exact Gaussian process
with a Matern-3/2 kernel. This release implementation reconstructs the same
model class with NumPy and SciPy so the figures do not require the heavyweight
BoTorch stack. Every saved figure contains one scientific data axis; color bars
are ancillary scale axes. Machine-readable source data and model provenance are
written beside the figure outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import curve_fit, minimize

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5 = ROOT / "data" / "raw" / "photokpfm_measurements.h5"
DEFAULT_OUTPUT = ROOT / "results"
HISTOGRAM_RANGE = (-2.5, -0.5)
HISTOGRAM_BINS = 30
CROP_PIXELS = 10
REPRESENTATIVE_POSITIONS = (1, 3, -2, -1)
SQRT3 = math.sqrt(3.0)
TRIANGLE_HEIGHT = SQRT3 / 2.0
COLORS = {
    "dark": "#243B63",
    "light": "#E0A12B",
    "fit_dark": "#2563A6",
    "fit_light": "#C94C4C",
    "ink": "#252A31",
    "muted": "#66717E",
    "grid": "#E4E7EB",
    "repeat_a": "#3B6FA5",
    "repeat_b": "#E6832A",
}


@dataclass(frozen=True)
class Acquisition:
    acquisition_position: int
    plate_index: int
    well_id: str
    x: float
    y: float
    pea_pct: float
    bda_pct: float
    fa_pct: float
    dark_name: str
    light_name: str
    dark_mean_v: float
    light_mean_v: float
    dark_mean_se_v: float
    light_mean_se_v: float
    delta_sp_v: float
    delta_sp_se_v: float


@dataclass(frozen=True)
class ExactGP:
    x_min: np.ndarray
    x_span: np.ndarray
    y_mean: float
    y_scale: float
    x_train_normalized: np.ndarray
    y_train_standardized: np.ndarray
    lengthscales: np.ndarray
    signal_sd: float
    noise_sd: float
    jitter: float
    cho: tuple[np.ndarray, bool]
    alpha: np.ndarray
    objective: float
    optimization_success: bool
    optimization_message: str
    optimization_starts: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate standalone Figure 3, Figure S1, and Figure S2 photoKPFM panels."
    )
    parser.add_argument("--input-h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "axes.titlesize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.4,
            "svg.fonttype": "none",
            "svg.hashsalt": "photokpfm-ternary-release",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_figure(figure: plt.Figure, stem: Path, dpi: int) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [stem.with_suffix(extension) for extension in (".png", ".svg", ".pdf")]
    figure.savefig(outputs[0], dpi=dpi, bbox_inches="tight", facecolor="white")
    figure.savefig(
        outputs[1],
        bbox_inches="tight",
        facecolor="white",
        metadata={"Date": None},
    )
    figure.savefig(
        outputs[2],
        bbox_inches="tight",
        facecolor="white",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(figure)
    return outputs


def decode_text(value: object) -> str:
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    return str(value)


def validate_h5(handle: h5py.File) -> None:
    expected = {
        "X": (91, 2),
        "X_train": (19, 2),
        "dark_data": (19, 128, 128),
        "light_data": (19, 128, 128),
        "idx": (19,),
        "dark_fn": (19,),
        "light_fn": (19,),
    }
    missing = sorted(set(expected).difference(handle.keys()))
    if missing:
        raise KeyError(f"Required HDF5 datasets are missing: {missing}")
    bad = {
        name: (tuple(handle[name].shape), shape)
        for name, shape in expected.items()
        if tuple(handle[name].shape) != shape
    }
    if bad:
        raise ValueError(f"Unexpected HDF5 dataset shapes: {bad}")


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    safe_sigma = np.copysign(max(abs(float(sigma)), 1e-12), sigma)
    return amplitude * np.exp(-((x - mean) ** 2) / (2.0 * safe_sigma**2))


def fit_histogram(values: np.ndarray) -> dict[str, object]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot fit an empty surface-potential distribution.")
    counts, edges = np.histogram(
        finite,
        bins=HISTOGRAM_BINS,
        range=HISTOGRAM_RANGE,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    if int(counts.max()) <= 0:
        raise ValueError("No surface-potential pixels fall inside the analysis range.")
    initial = [float(counts.max()), float(np.mean(finite)), max(float(np.std(finite)), 1e-6)]
    params, covariance = curve_fit(
        gaussian,
        centers,
        counts,
        p0=initial,
        maxfev=5000,
    )
    mean_se = float(np.sqrt(max(float(covariance[1, 1]), 0.0)))
    if not np.all(np.isfinite(params)) or not np.isfinite(mean_se):
        raise ValueError("The Gaussian surface-potential fit returned non-finite values.")
    return {
        "counts": counts.astype(int),
        "edges": edges,
        "centers": centers,
        "amplitude": float(params[0]),
        "mean": float(params[1]),
        "sigma": abs(float(params[2])),
        "mean_se": mean_se,
    }


def cartesian_to_composition(points: np.ndarray) -> np.ndarray:
    array = np.asarray(points, dtype=float)
    fa = array[:, 1] / TRIANGLE_HEIGHT
    pea = array[:, 0] - 0.5 * fa
    bda = 1.0 - pea - fa
    components = np.column_stack([pea, bda, fa])
    components[np.abs(components) < 1e-12] = 0.0
    return components


def plate_index_to_well(index: int) -> str:
    if not 0 <= index < 96:
        raise ValueError(f"Invalid zero-based plate index: {index}")
    return f"{'ABCDEFGH'[index % 8]}{index // 8 + 1}"


def load_acquisitions(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Acquisition]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as handle:
        validate_h5(handle)
        prediction_points = np.asarray(handle["X"], dtype=float)
        training_points = np.asarray(handle["X_train"], dtype=float)
        plate_indices = np.asarray(handle["idx"], dtype=int)
        dark_images = np.asarray(handle["dark_data"], dtype=float)
        light_images = np.asarray(handle["light_data"], dtype=float)
        dark_names = [decode_text(value) for value in np.asarray(handle["dark_fn"])]
        light_names = [decode_text(value) for value in np.asarray(handle["light_fn"])]
    components = cartesian_to_composition(training_points) * 100.0
    acquisitions: list[Acquisition] = []
    for position in range(len(plate_indices)):
        dark_values = dark_images[position, CROP_PIXELS:-CROP_PIXELS, CROP_PIXELS:-CROP_PIXELS].ravel()
        light_values = light_images[position, CROP_PIXELS:-CROP_PIXELS, CROP_PIXELS:-CROP_PIXELS].ravel()
        dark_fit = fit_histogram(dark_values)
        light_fit = fit_histogram(light_values)
        dark_mean = float(dark_fit["mean"])
        light_mean = float(light_fit["mean"])
        delta_se = math.hypot(float(dark_fit["mean_se"]), float(light_fit["mean_se"]))
        acquisitions.append(
            Acquisition(
                acquisition_position=position,
                plate_index=int(plate_indices[position]),
                well_id=plate_index_to_well(int(plate_indices[position])),
                x=float(training_points[position, 0]),
                y=float(training_points[position, 1]),
                pea_pct=float(components[position, 0]),
                bda_pct=float(components[position, 1]),
                fa_pct=float(components[position, 2]),
                dark_name=dark_names[position],
                light_name=light_names[position],
                dark_mean_v=dark_mean,
                light_mean_v=light_mean,
                dark_mean_se_v=float(dark_fit["mean_se"]),
                light_mean_se_v=float(light_fit["mean_se"]),
                delta_sp_v=light_mean - dark_mean,
                delta_sp_se_v=delta_se,
            )
        )
    y_train = np.asarray([row.delta_sp_v for row in acquisitions], dtype=float)
    return prediction_points, training_points, y_train, acquisitions


def matern32_kernel(
    left: np.ndarray,
    right: np.ndarray,
    lengthscales: np.ndarray,
    signal_sd: float,
) -> np.ndarray:
    scaled = (left[:, None, :] - right[None, :, :]) / lengthscales[None, None, :]
    distance = np.sqrt(np.sum(np.square(scaled), axis=2))
    scaled_distance = SQRT3 * distance
    return signal_sd**2 * (1.0 + scaled_distance) * np.exp(-scaled_distance)


def gp_covariance(
    x_train: np.ndarray,
    lengthscales: np.ndarray,
    signal_sd: float,
    noise_sd: float,
    jitter: float,
) -> np.ndarray:
    covariance = matern32_kernel(x_train, x_train, lengthscales, signal_sd)
    diagonal = noise_sd**2 + jitter
    return covariance + np.eye(len(x_train), dtype=float) * diagonal


def gp_objective(
    theta: np.ndarray,
    x_train: np.ndarray,
    y_train: np.ndarray,
    jitter: float,
) -> float:
    lengthscales = np.exp(theta[:2])
    signal_sd = float(np.exp(theta[2]))
    noise_sd = float(np.exp(theta[3]))
    try:
        covariance = gp_covariance(x_train, lengthscales, signal_sd, noise_sd, jitter)
        factor = cho_factor(covariance, lower=True, check_finite=False)
        alpha = cho_solve(factor, y_train, check_finite=False)
    except (ValueError, np.linalg.LinAlgError):
        return 1e30
    log_determinant = 2.0 * float(np.log(np.diag(factor[0])).sum())
    value = 0.5 * float(y_train @ alpha)
    value += 0.5 * log_determinant
    value += 0.5 * len(y_train) * math.log(2.0 * math.pi)
    return value if np.isfinite(value) else 1e30


def fit_exact_gp(x_train: np.ndarray, y_train: np.ndarray) -> ExactGP:
    x_min = np.min(x_train, axis=0)
    x_span = np.ptp(x_train, axis=0)
    if np.any(x_span <= 0):
        raise ValueError("Training coordinates cannot be normalized because one dimension is constant.")
    x_normalized = (x_train - x_min) / x_span
    y_mean = float(np.mean(y_train))
    y_scale = float(np.std(y_train, ddof=1))
    if not np.isfinite(y_scale) or y_scale <= 0:
        raise ValueError("Training responses cannot be standardized.")
    y_standardized = (y_train - y_mean) / y_scale
    jitter = 1e-9
    lower = np.log(np.asarray([0.03, 0.03, 0.05, 0.002], dtype=float))
    upper = np.log(np.asarray([0.70, 0.70, 5.00, 2.000], dtype=float))
    starts = (
        (0.12, 0.12, 1.0, 0.10),
        (0.25, 0.25, 1.0, 0.25),
        (0.50, 0.20, 1.2, 0.50),
        (0.20, 0.50, 0.8, 0.75),
        (0.60, 0.60, 1.5, 1.00),
    )
    results = []
    for start in starts:
        result = minimize(
            gp_objective,
            np.log(np.asarray(start, dtype=float)),
            args=(x_normalized, y_standardized, jitter),
            method="L-BFGS-B",
            bounds=list(zip(lower, upper)),
            options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
        )
        if np.isfinite(result.fun):
            results.append(result)
    if not results:
        raise RuntimeError("Gaussian-process hyperparameter optimization failed for every start.")
    best = min(results, key=lambda item: float(item.fun))
    lengthscales = np.exp(best.x[:2])
    signal_sd = float(np.exp(best.x[2]))
    noise_sd = float(np.exp(best.x[3]))
    covariance = gp_covariance(x_normalized, lengthscales, signal_sd, noise_sd, jitter)
    factor = cho_factor(covariance, lower=True, check_finite=False)
    alpha = cho_solve(factor, y_standardized, check_finite=False)
    return ExactGP(
        x_min=x_min,
        x_span=x_span,
        y_mean=y_mean,
        y_scale=y_scale,
        x_train_normalized=x_normalized,
        y_train_standardized=y_standardized,
        lengthscales=lengthscales,
        signal_sd=signal_sd,
        noise_sd=noise_sd,
        jitter=jitter,
        cho=factor,
        alpha=alpha,
        objective=float(best.fun),
        optimization_success=bool(best.success),
        optimization_message=str(best.message),
        optimization_starts=len(starts),
    )


def predict_exact_gp(model: ExactGP, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normalized = (np.asarray(points, dtype=float) - model.x_min) / model.x_span
    cross = matern32_kernel(
        model.x_train_normalized,
        normalized,
        model.lengthscales,
        model.signal_sd,
    )
    mean_standardized = cross.T @ model.alpha
    solved = cho_solve(model.cho, cross, check_finite=False)
    latent_variance_standardized = model.signal_sd**2 - np.sum(cross * solved, axis=0)
    latent_variance_standardized = np.maximum(latent_variance_standardized, 0.0)
    mean = model.y_mean + model.y_scale * mean_standardized
    standard_deviation = model.y_scale * np.sqrt(latent_variance_standardized)
    return mean, standard_deviation


def draw_ternary(axis: plt.Axes) -> None:
    axis.plot([0.0, 1.0, 0.5, 0.0], [0.0, 0.0, TRIANGLE_HEIGHT, 0.0], color=COLORS["ink"], lw=1.0)
    for fraction in (0.2, 0.4, 0.6, 0.8):
        axis.plot(
            [0.5 * fraction, 1.0 - 0.5 * fraction],
            [TRIANGLE_HEIGHT * fraction] * 2,
            color=COLORS["grid"],
            lw=0.55,
            zorder=0,
        )
        axis.plot(
            [fraction, 0.5 + 0.5 * fraction],
            [0.0, TRIANGLE_HEIGHT * (1.0 - fraction)],
            color=COLORS["grid"],
            lw=0.55,
            zorder=0,
        )
        axis.plot(
            [1.0 - fraction, 0.5 * (1.0 - fraction)],
            [0.0, TRIANGLE_HEIGHT * (1.0 - fraction)],
            color=COLORS["grid"],
            lw=0.55,
            zorder=0,
        )
    axis.text(-0.015, -0.055, r"BDAPbI$_4$", ha="center", va="top", fontsize=9.0)
    axis.text(1.015, -0.055, r"PEA$_2$PbI$_4$", ha="center", va="top", fontsize=9.0)
    axis.text(0.5, TRIANGLE_HEIGHT + 0.035, r"FAPbI$_3$", ha="center", va="bottom", fontsize=9.0)
    axis.set_xlim(-0.08, 1.08)
    axis.set_ylim(-0.10, TRIANGLE_HEIGHT + 0.09)
    axis.set_aspect("equal")
    axis.axis("off")


def aggregate_measurements(acquisitions: list[Acquisition]) -> list[dict[str, float | int | str]]:
    grouped: dict[int, list[Acquisition]] = {}
    for row in acquisitions:
        grouped.setdefault(row.plate_index, []).append(row)
    output: list[dict[str, float | int | str]] = []
    for plate_index, rows in grouped.items():
        output.append(
            {
                "plate_index": plate_index,
                "well_id": rows[0].well_id,
                "x": rows[0].x,
                "y": rows[0].y,
                "mean_delta_sp_v": float(np.mean([row.delta_sp_v for row in rows])),
                "propagated_fit_se_of_mean_v": float(
                    math.sqrt(sum(row.delta_sp_se_v**2 for row in rows)) / len(rows)
                ),
                "repeat_count": len(rows),
            }
        )
    return output


def plot_gp_map(
    points: np.ndarray,
    values: np.ndarray,
    acquisitions: list[Acquisition],
    stem: Path,
    colorbar_label: str,
    cmap: str,
    dpi: int,
    show_fit_uncertainty: bool,
) -> list[Path]:
    figure, axis = plt.subplots(figsize=(5.3, 4.6))
    draw_ternary(axis)
    lower = (
        float(min(0.0, np.nanmin(values)))
        if show_fit_uncertainty
        else float(np.nanmin(values))
    )
    upper = float(np.nanmax(values))
    if show_fit_uncertainty:
        upper = max(0.15, upper)
    if not np.isfinite(upper) or upper <= lower:
        upper = lower + 1.0
    image = axis.scatter(
        points[:, 0],
        points[:, 1],
        c=values,
        cmap=cmap,
        vmin=lower,
        vmax=upper,
        s=68,
        marker="h",
        linewidths=0,
        zorder=2,
    )
    aggregated = aggregate_measurements(acquisitions)
    measured_x = np.asarray([float(row["x"]) for row in aggregated])
    measured_y = np.asarray([float(row["y"]) for row in aggregated])
    axis.scatter(
        measured_x,
        measured_y,
        marker="x",
        c="black",
        s=30,
        linewidths=1.1,
        zorder=5,
    )
    legend_items = [
        Line2D([], [], marker="x", color="black", linestyle="none", markersize=5.5, label="Measured composition")
    ]
    if show_fit_uncertainty:
        fit_se = np.asarray([float(row["propagated_fit_se_of_mean_v"]) for row in aggregated])
        if float(np.ptp(fit_se)) > 0:
            sizes = 38.0 + 175.0 * (fit_se - fit_se.min()) / np.ptp(fit_se)
        else:
            sizes = np.full_like(fit_se, 90.0)
        axis.scatter(
            measured_x,
            measured_y,
            s=sizes,
            facecolors="none",
            edgecolors="white",
            linewidths=1.1,
            zorder=4,
        )
        axis.scatter(
            measured_x,
            measured_y,
            s=sizes,
            facecolors="none",
            edgecolors="black",
            linewidths=0.45,
            zorder=4,
        )
        legend_items.append(
            Line2D(
                [],
                [],
                marker="o",
                markerfacecolor="none",
                markeredgecolor="black",
                linestyle="none",
                markersize=7.0,
                label="Circle area ∝ propagated fit SE",
            )
        )
        selected_positions = [
            position if position >= 0 else len(acquisitions) + position
            for position in REPRESENTATIVE_POSITIONS
        ]
        annotation_offsets = ((0.025, 0.028), (-0.045, 0.030), (0.022, 0.032), (-0.040, 0.032))
        for panel_number, (position, offset) in enumerate(
            zip(selected_positions, annotation_offsets), start=1
        ):
            row = acquisitions[position]
            axis.text(
                row.x + offset[0],
                row.y + offset[1],
                str(panel_number),
                ha="center",
                va="center",
                fontsize=8.0,
                fontweight="bold",
                color="black",
                bbox={"boxstyle": "circle,pad=0.12", "facecolor": "white", "edgecolor": "black", "linewidth": 0.55},
                zorder=7,
            )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.046, pad=0.035)
    colorbar.set_label(colorbar_label)
    axis.legend(
        handles=legend_items,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.10),
        frameon=False,
        ncol=len(legend_items),
        handletextpad=0.45,
        columnspacing=1.2,
    )
    figure.subplots_adjust(left=0.03, right=0.88, bottom=0.14, top=0.97)
    return save_figure(figure, stem, dpi)


def histogram_panel(
    panel_number: int,
    position: int,
    dark_image: np.ndarray,
    light_image: np.ndarray,
    acquisition: Acquisition,
    stem: Path,
    dpi: int,
) -> tuple[list[Path], list[dict[str, object]]]:
    dark_values = dark_image[CROP_PIXELS:-CROP_PIXELS, CROP_PIXELS:-CROP_PIXELS].ravel()
    light_values = light_image[CROP_PIXELS:-CROP_PIXELS, CROP_PIXELS:-CROP_PIXELS].ravel()
    dark_fit = fit_histogram(dark_values)
    light_fit = fit_histogram(light_values)
    centers = np.asarray(dark_fit["centers"], dtype=float)
    width = float(np.diff(np.asarray(dark_fit["edges"], dtype=float))[0])
    figure, axis = plt.subplots(figsize=(4.25, 3.25))
    axis.bar(
        centers,
        np.asarray(dark_fit["counts"], dtype=float),
        width=width,
        color=COLORS["dark"],
        edgecolor="white",
        linewidth=0.25,
        alpha=0.66,
        label="Dark",
    )
    axis.bar(
        centers,
        np.asarray(light_fit["counts"], dtype=float),
        width=width,
        color=COLORS["light"],
        edgecolor="white",
        linewidth=0.25,
        alpha=0.62,
        label="Illuminated",
    )
    fit_x = np.linspace(HISTOGRAM_RANGE[0], HISTOGRAM_RANGE[1], 400)
    dark_curve = gaussian(
        fit_x,
        float(dark_fit["amplitude"]),
        float(dark_fit["mean"]),
        float(dark_fit["sigma"]),
    )
    light_curve = gaussian(
        fit_x,
        float(light_fit["amplitude"]),
        float(light_fit["mean"]),
        float(light_fit["sigma"]),
    )
    axis.plot(fit_x, dark_curve, color=COLORS["fit_dark"], linestyle="--", label="Dark Gaussian fit")
    axis.plot(
        fit_x,
        light_curve,
        color=COLORS["fit_light"],
        linestyle="--",
        label="Illuminated Gaussian fit",
    )
    axis.set_xlim(*HISTOGRAM_RANGE)
    axis.set_xlabel("Surface potential (V)")
    axis.set_ylabel("Pixel count")
    axis.set_title(
        rf"{panel_number}. {acquisition.well_id}: $\Delta SP$ = {acquisition.delta_sp_v:.3f} $\pm$ {acquisition.delta_sp_se_v:.3f} V"
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, ncol=2, fontsize=7.2, loc="upper left")
    figure.tight_layout()
    paths = save_figure(figure, stem, dpi)
    source_rows: list[dict[str, object]] = []
    dark_counts = np.asarray(dark_fit["counts"], dtype=int)
    light_counts = np.asarray(light_fit["counts"], dtype=int)
    for bin_index, center in enumerate(centers):
        source_rows.append(
            {
                "panel_number": panel_number,
                "selected_position_expression": str(position),
                "acquisition_position_zero_based": acquisition.acquisition_position,
                "plate_index_zero_based": acquisition.plate_index,
                "well_id": acquisition.well_id,
                "PEA2PbI4_pct": acquisition.pea_pct,
                "BDAPbI4_pct": acquisition.bda_pct,
                "FAPbI3_pct": acquisition.fa_pct,
                "histogram_bin_center_v": float(center),
                "dark_pixel_count": int(dark_counts[bin_index]),
                "illuminated_pixel_count": int(light_counts[bin_index]),
                "dark_gaussian_fit_at_bin": float(
                    gaussian(
                        np.asarray([center]),
                        float(dark_fit["amplitude"]),
                        float(dark_fit["mean"]),
                        float(dark_fit["sigma"]),
                    )[0]
                ),
                "illuminated_gaussian_fit_at_bin": float(
                    gaussian(
                        np.asarray([center]),
                        float(light_fit["amplitude"]),
                        float(light_fit["mean"]),
                        float(light_fit["sigma"]),
                    )[0]
                ),
                "SP_dark_v": acquisition.dark_mean_v,
                "SP_dark_fit_se_v": acquisition.dark_mean_se_v,
                "SP_illuminated_v": acquisition.light_mean_v,
                "SP_illuminated_fit_se_v": acquisition.light_mean_se_v,
                "delta_SP_v": acquisition.delta_sp_v,
                "delta_SP_fit_se_v": acquisition.delta_sp_se_v,
            }
        )
    return paths, source_rows


def plot_repeat_variability(
    acquisitions: list[Acquisition],
    stem: Path,
    dpi: int,
) -> list[Path]:
    selected = [row for row in acquisitions if row.plate_index in {0, 78}]
    counts = {index: sum(row.plate_index == index for row in selected) for index in (0, 78)}
    if counts != {0: 4, 78: 4}:
        raise ValueError(f"Expected four repeated acquisitions each for plate indices 0 and 78, found {counts}.")
    figure, axis = plt.subplots(figsize=(4.5, 3.35))
    jitter = np.asarray([-0.12, -0.04, 0.04, 0.12], dtype=float)
    for category_position, plate_index in enumerate((0, 78)):
        rows = [row for row in selected if row.plate_index == plate_index]
        values = np.asarray([row.delta_sp_v for row in rows], dtype=float)
        errors = np.asarray([row.delta_sp_se_v for row in rows], dtype=float)
        color = COLORS["repeat_a"] if plate_index == 0 else COLORS["repeat_b"]
        axis.errorbar(
            category_position + jitter,
            values,
            yerr=errors,
            fmt="o",
            ms=5.5,
            color=color,
            ecolor=color,
            elinewidth=1.0,
            capsize=2.0,
            label=f"Index {plate_index} acquisitions",
            zorder=3,
        )
        mean_value = float(np.mean(values))
        axis.hlines(
            mean_value,
            category_position - 0.20,
            category_position + 0.20,
            color=COLORS["ink"],
            linewidth=1.4,
            zorder=2,
        )
    axis.axhline(0.0, color=COLORS["muted"], linewidth=0.7, linestyle=":", zorder=1)
    axis.set_xticks([0, 1])
    axis.set_xticklabels([f"Index 0\n({plate_index_to_well(0)})", f"Index 78\n({plate_index_to_well(78)})"])
    axis.set_ylabel(r"Fitted $\Delta SP$ (V)")
    axis.set_xlabel("Repeatedly sampled plate location")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, loc="best", fontsize=7.5)
    figure.tight_layout()
    return save_figure(figure, stem, dpi)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = fieldnames or (list(rows[0]) if rows else [])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def relative_label(path: Path, base: Path) -> str:
    try:
        return path.resolve().relative_to(base.resolve()).as_posix()
    except ValueError:
        return path.name


def main() -> None:
    args = parse_args()
    input_h5 = args.input_h5.resolve()
    output_root = args.output_root.resolve()
    if args.dpi < 150:
        raise ValueError("DPI must be at least 150 for manuscript output.")
    main_dir = output_root / "figures" / "main"
    supplementary_dir = output_root / "figures" / "supplementary"
    source_dir = output_root / "source_data" / "ternary_photokpfm"
    for directory in (main_dir, supplementary_dir, source_dir):
        directory.mkdir(parents=True, exist_ok=True)
    configure_plotting()
    prediction_points, training_points, y_train, acquisitions = load_acquisitions(input_h5)
    gp_model = fit_exact_gp(training_points, y_train)
    gp_mean, gp_standard_deviation = predict_exact_gp(gp_model, prediction_points)
    all_outputs: list[Path] = []
    all_outputs.extend(
        plot_gp_map(
            prediction_points,
            gp_mean,
            acquisitions,
            main_dir / "Figure3_GP_mean_map",
            r"GP mean $\Delta SP$ (V)",
            "viridis",
            args.dpi,
            True,
        )
    )
    all_outputs.extend(
        plot_gp_map(
            prediction_points,
            gp_standard_deviation,
            acquisitions,
            supplementary_dir / "FigureS2_GP_uncertainty_map",
            r"GP posterior standard deviation of $\Delta SP$ (V)",
            "magma",
            args.dpi,
            False,
        )
    )
    with h5py.File(input_h5, "r") as handle:
        dark_images = np.asarray(handle["dark_data"], dtype=float)
        light_images = np.asarray(handle["light_data"], dtype=float)
    histogram_source_rows: list[dict[str, object]] = []
    normalized_positions = [position if position >= 0 else len(acquisitions) + position for position in REPRESENTATIVE_POSITIONS]
    for panel_number, (expression, position) in enumerate(
        zip(REPRESENTATIVE_POSITIONS, normalized_positions), start=1
    ):
        panel_paths, source_rows = histogram_panel(
            panel_number,
            expression,
            dark_images[position],
            light_images[position],
            acquisitions[position],
            main_dir / f"Figure3_representative_SP_histogram_{panel_number}",
            args.dpi,
        )
        all_outputs.extend(panel_paths)
        histogram_source_rows.extend(source_rows)
    all_outputs.extend(
        plot_repeat_variability(
            acquisitions,
            supplementary_dir / "FigureS1_repeat_location_variability",
            args.dpi,
        )
    )
    acquisition_rows = [
        {
            "acquisition_position_zero_based": row.acquisition_position,
            "plate_index_zero_based": row.plate_index,
            "well_id": row.well_id,
            "ternary_x": row.x,
            "ternary_y": row.y,
            "PEA2PbI4_pct": row.pea_pct,
            "BDAPbI4_pct": row.bda_pct,
            "FAPbI3_pct": row.fa_pct,
            "dark_acquisition": row.dark_name,
            "illuminated_acquisition": row.light_name,
            "SP_dark_v": row.dark_mean_v,
            "SP_dark_fit_se_v": row.dark_mean_se_v,
            "SP_illuminated_v": row.light_mean_v,
            "SP_illuminated_fit_se_v": row.light_mean_se_v,
            "delta_SP_v": row.delta_sp_v,
            "delta_SP_fit_se_v": row.delta_sp_se_v,
            "repeat_location_for_FigureS1": row.plate_index in {0, 78},
        }
        for row in acquisitions
    ]
    grid_components = cartesian_to_composition(prediction_points) * 100.0
    grid_rows = [
        {
            "grid_position_zero_based": position,
            "ternary_x": float(prediction_points[position, 0]),
            "ternary_y": float(prediction_points[position, 1]),
            "PEA2PbI4_pct": float(grid_components[position, 0]),
            "BDAPbI4_pct": float(grid_components[position, 1]),
            "FAPbI3_pct": float(grid_components[position, 2]),
            "GP_mean_delta_SP_v": float(gp_mean[position]),
            "GP_posterior_standard_deviation_v": float(gp_standard_deviation[position]),
        }
        for position in range(len(prediction_points))
    ]
    figure_triplets = len(all_outputs) // 3
    if len(all_outputs) % 3 != 0:
        raise RuntimeError("Figure output count is not divisible into PNG/SVG/PDF triplets.")
    acquisition_path = source_dir / "photoKPFM_acquisitions_all_19.csv"
    grid_path = source_dir / "Figure3_and_FigureS2_GP_grid_source_data.csv"
    histogram_path = source_dir / "Figure3_representative_histograms_source_data.csv"
    repeat_path = source_dir / "FigureS1_repeat_location_variability_source_data.csv"
    write_csv(acquisition_path, acquisition_rows)
    write_csv(grid_path, grid_rows)
    write_csv(histogram_path, histogram_source_rows)
    write_csv(
        repeat_path,
        [row for row in acquisition_rows if bool(row["repeat_location_for_FigureS1"])],
    )
    all_outputs.extend([acquisition_path, grid_path, histogram_path, repeat_path])
    metadata = {
        "schema_version": "1.0",
        "source_h5": relative_label(input_h5, ROOT),
        "source_h5_sha256": sha256(input_h5),
        "historical_notebook_provenance": {
            "name": "20260323_photoKPFM_ternary.ipynb",
            "use": "Read-only provenance for measured-index selection, histogram analysis, and model class.",
            "code_copied_into_release": False,
            "historical_model": "BoTorch/GPyTorch exact GP, Matern-3/2 ARD kernel, normalized inputs, standardized response.",
        },
        "surface_potential_extraction": {
            "crop_pixels_each_edge": CROP_PIXELS,
            "histogram_bins": HISTOGRAM_BINS,
            "histogram_range_v": list(HISTOGRAM_RANGE),
            "fit_function": "A * exp(-(x - mean)^2 / (2 * sigma^2))",
            "response": "delta_SP = illuminated Gaussian mean - dark Gaussian mean",
            "response_uncertainty": "sqrt(SE_dark_mean^2 + SE_illuminated_mean^2)",
            "fit_uncertainty_source": "square root of Gaussian-fit covariance diagonal for the mean parameter",
        },
        "gaussian_process": {
            "implementation": "Transparent exact GP implemented with NumPy and SciPy",
            "kernel": "ARD Matern-3/2",
            "input_transform": "per-coordinate min-max normalization from the 19 training acquisitions",
            "outcome_transform": "sample-mean centering and sample-standard-deviation scaling",
            "hyperparameter_fit": "deterministic five-start L-BFGS-B maximum log marginal likelihood",
            "lengthscale_bounds_normalized_units": [0.03, 0.70],
            "signal_sd_bounds_standardized_units": [0.05, 5.0],
            "noise_sd_bounds_standardized_units": [0.002, 2.0],
            "lengthscales_normalized_units": gp_model.lengthscales.tolist(),
            "signal_sd_standardized_units": gp_model.signal_sd,
            "noise_sd_standardized_units": gp_model.noise_sd,
            "training_response_mean_v": gp_model.y_mean,
            "training_response_sample_sd_v": gp_model.y_scale,
            "negative_log_marginal_likelihood": gp_model.objective,
            "optimizer_success": gp_model.optimization_success,
            "optimizer_message": gp_model.optimization_message,
            "optimizer_starts": gp_model.optimization_starts,
            "jitter_standardized_variance": gp_model.jitter,
            "reported_uncertainty": "latent GP posterior standard deviation, excluding observation noise",
            "compatibility_caveat": "This SciPy reconstruction matches the historical GP model class but is not expected to be bitwise identical to BoTorch/GPyTorch hyperparameter fitting.",
        },
        "representative_acquisition_position_expressions": list(REPRESENTATIVE_POSITIONS),
        "representative_acquisition_positions_zero_based": normalized_positions,
        "repeat_location_plate_indices_zero_based": [0, 78],
        "panel_policy": "Every figure triplet has exactly one scientific data axis; color bars are ancillary scale axes.",
        "outputs": [relative_label(path, ROOT) for path in all_outputs],
    }
    metadata_path = source_dir / "ternary_GP_model_and_figure_provenance.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    all_outputs.append(metadata_path)
    manifest_rows = [
        {
            "path": relative_label(path, ROOT),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(all_outputs, key=lambda item: item.as_posix())
    ]
    manifest_path = source_dir / "ternary_photokpfm_output_manifest.csv"
    write_csv(manifest_path, manifest_rows)
    print(f"Input HDF5: {relative_label(input_h5, ROOT)}")
    print(f"Acquisitions: {len(acquisitions)}; unique plate locations: {len(set(row.plate_index for row in acquisitions))}")
    print(f"GP lengthscales: {gp_model.lengthscales.tolist()}")
    print(f"GP signal SD: {gp_model.signal_sd:.8g}; noise SD: {gp_model.noise_sd:.8g}")
    print(f"Figure triplets: {figure_triplets}")
    print(f"Output manifest: {relative_label(manifest_path, ROOT)}")


if __name__ == "__main__":
    main()
