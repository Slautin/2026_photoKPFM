#!/usr/bin/env python3
"""Validate PL/photo-KPFM relationships against the matched experimental plate.

The phase-fitting CSV was created in Colab but is not always downloaded. This
script can recover the displayed initial-timepoint table from the saved notebook
without rerunning the fit. It then reports raw and composition-adjusted
associations. The output is deliberately associative, not causal.
"""

from __future__ import annotations

import argparse
import json
import math
from io import StringIO
from pathlib import Path
from typing import Iterable

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]

KPFM_METRICS = [
    "photovoltage_signed_V",
    "photovoltage_abs_V",
    "mu_dark_V",
    "mu_light_V",
]

PHASE_COLUMNS = [
    "frac_n=1-like / 2D",
    "frac_n=2-like",
    "frac_n=3-like",
    "frac_n=4-like",
    "frac_n>=5 / high-n",
    "frac_3D / 3D-like",
]

PL_METRICS = [
    "Top_PL_peak_intensity",
    "dominant_peak_nm",
    *PHASE_COLUMNS,
    "frac_low_n_1_to_4",
]

DISPLAY_LABELS = {
    "photovoltage_signed_V": "signed dV",
    "photovoltage_abs_V": "|dV|",
    "mu_dark_V": "dark SP",
    "mu_light_V": "light SP",
    "Top_PL_peak_intensity": "PL peak intensity",
    "dominant_peak_nm": "dominant peak",
    "frac_n=1-like / 2D": "n=1 / 2D",
    "frac_n=2-like": "n=2",
    "frac_n=3-like": "n=3",
    "frac_n=4-like": "n=4",
    "frac_n>=5 / high-n": "high-n",
    "frac_3D / 3D-like": "3D-like",
    "frac_low_n_1_to_4": "n=1-4 sum",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze matched PL phase and photo-KPFM measurements."
    )
    parser.add_argument("--matched-csv", required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--phase-csv", type=Path)
    source.add_argument("--phase-notebook", type=Path)
    parser.add_argument("--literature-relationships", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260816)
    return parser.parse_args()


def _as_text(value: object) -> str:
    if isinstance(value, list):
        return "".join(str(item) for item in value)
    return "" if value is None else str(value)


def _read_html_table(html: str) -> pd.DataFrame:
    tables = pd.read_html(StringIO(html))
    if not tables:
        raise ValueError("Notebook display output did not contain a table.")
    frame = tables[0]
    frame = frame.loc[:, ~frame.columns.astype(str).str.startswith("Unnamed:")]
    return frame


def recover_initial_phase_table(notebook_path: Path) -> pd.DataFrame:
    """Recover the first initial-timepoint phase/KPFM display table."""
    with notebook_path.open("r", encoding="utf-8") as handle:
        notebook = json.load(handle)

    waiting_for_initial_table = False
    candidates: list[pd.DataFrame] = []

    for cell in notebook.get("cells", []):
        for output in cell.get("outputs", []):
            stream_text = _as_text(output.get("text"))
            if "Final plotting table for initial timepoint" in stream_text:
                waiting_for_initial_table = True
                continue

            data = output.get("data", {})
            html = _as_text(data.get("text/html"))
            if not waiting_for_initial_table or "<table" not in html:
                continue

            frame = _read_html_table(html)
            waiting_for_initial_table = False
            if {"Well", "dominant_peak_nm", *PHASE_COLUMNS}.issubset(frame.columns):
                candidates.append(frame)

    for frame in candidates:
        if "photovoltage_abs_V" in frame.columns:
            return frame.copy()

    if candidates:
        return candidates[0].copy()
    raise ValueError(
        "Could not find an initial-timepoint phase table in the notebook outputs."
    )


def load_phase_table(args: argparse.Namespace) -> tuple[pd.DataFrame, str]:
    if args.phase_csv:
        return pd.read_csv(args.phase_csv), portable_path(args.phase_csv)
    return recover_initial_phase_table(args.phase_notebook), portable_path(
        args.phase_notebook
    )


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.name


def prepare_experimental_table(
    matched_csv: Path, phase_table: pd.DataFrame
) -> pd.DataFrame:
    matched = pd.read_csv(matched_csv)
    if "Well" not in matched and "well_id" in matched:
        matched = matched.rename(columns={"well_id": "Well"})
    if "Well" not in phase_table and "well_id" in phase_table:
        phase_table = phase_table.rename(columns={"well_id": "Well"})
    if "Well" not in matched or "Well" not in phase_table:
        raise ValueError("Both inputs must contain a Well or well_id column.")

    phase_keep = [
        column
        for column in [
            "Well",
            *KPFM_METRICS,
            "Top_PL_peak_intensity",
            "PEA_pct",
            "BDA_pct",
            "FA_pct",
            "plot_order",
            "num_peaks",
            "dominant_peak_nm",
            "dominant_phase",
            *PHASE_COLUMNS,
        ]
        if column in phase_table.columns
    ]
    phase = phase_table[phase_keep].copy()
    merged = matched.merge(
        phase,
        on="Well",
        how="inner",
        validate="one_to_one",
        suffixes=("", "_phase"),
    )
    for column in phase_keep:
        phase_column = f"{column}_phase"
        if column == "Well" or phase_column not in merged:
            continue
        if column not in merged:
            merged[column] = merged[phase_column]
        else:
            merged[column] = merged[column].where(merged[column].notna(), merged[phase_column])
        merged = merged.drop(columns=phase_column)

    numeric_columns = [
        *KPFM_METRICS,
        "Top_PL_peak_intensity",
        "PEA_pct",
        "BDA_pct",
        "FA_pct",
        "num_peaks",
        "dominant_peak_nm",
        *PHASE_COLUMNS,
    ]
    for column in numeric_columns:
        if column in merged:
            merged[column] = pd.to_numeric(merged[column], errors="coerce")

    missing = [column for column in KPFM_METRICS + PHASE_COLUMNS if column not in merged]
    if missing:
        raise ValueError(f"Required experimental columns are missing: {missing}")

    merged["frac_low_n_1_to_4"] = merged[PHASE_COLUMNS[:4]].sum(axis=1)
    merged["frac_high_n_plus_3D"] = merged[PHASE_COLUMNS[4:]].sum(axis=1)
    merged["phase_fraction_sum"] = merged[PHASE_COLUMNS].sum(axis=1)
    merged["pl_near_upper_limit"] = merged["Top_PL_peak_intensity"] >= 98000
    return merged.sort_values("plot_order" if "plot_order" in merged else "Well")


def _valid_pair(frame: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    pair = frame[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    return pair


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    if len(x) < 4 or np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return math.nan, math.nan
    if method == "pearson":
        result = stats.pearsonr(x, y)
    else:
        result = stats.spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def bootstrap_correlation_ci(
    x: np.ndarray,
    y: np.ndarray,
    method: str,
    samples: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    values: list[float] = []
    n = len(x)
    for _ in range(samples):
        index = rng.integers(0, n, n)
        bx = x[index]
        by = y[index]
        if np.nanstd(bx) == 0 or np.nanstd(by) == 0:
            continue
        coefficient, _ = _safe_corr(bx, by, method)
        if np.isfinite(coefficient):
            values.append(coefficient)
    if len(values) < max(100, samples // 10):
        return math.nan, math.nan, len(values)
    return (
        float(np.percentile(values, 2.5)),
        float(np.percentile(values, 97.5)),
        len(values),
    )


def _residualize(values: np.ndarray, controls: np.ndarray) -> tuple[np.ndarray, int]:
    design = np.column_stack([np.ones(len(values)), controls])
    rank = int(np.linalg.matrix_rank(design))
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    return values - fitted, rank - 1


def partial_correlation(
    frame: pd.DataFrame,
    x_col: str,
    y_col: str,
    controls: Iterable[str],
    method: str,
) -> tuple[float, float, int, int]:
    columns = [x_col, y_col, *controls]
    subset = frame[columns].replace([np.inf, -np.inf], np.nan).dropna()
    n = len(subset)
    if n < 6:
        return math.nan, math.nan, n, 0

    x = subset[x_col].to_numpy(dtype=float)
    y = subset[y_col].to_numpy(dtype=float)
    control_values = subset[list(controls)].to_numpy(dtype=float)
    if method == "spearman":
        x = stats.rankdata(x)
        y = stats.rankdata(y)
        control_values = np.column_stack(
            [stats.rankdata(control_values[:, i]) for i in range(control_values.shape[1])]
        )

    x_residual, effective_controls = _residualize(x, control_values)
    y_residual, _ = _residualize(y, control_values)
    coefficient, _ = _safe_corr(x_residual, y_residual, "pearson")
    degrees_freedom = n - effective_controls - 2
    if not np.isfinite(coefficient) or degrees_freedom <= 0 or abs(coefficient) >= 1:
        p_value = 0.0 if abs(coefficient) == 1 else math.nan
    else:
        statistic = coefficient * math.sqrt(degrees_freedom / (1 - coefficient**2))
        p_value = float(2 * stats.t.sf(abs(statistic), degrees_freedom))
    return coefficient, p_value, n, effective_controls


def leave_one_out_range(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    values = []
    for index in range(len(x)):
        keep = np.arange(len(x)) != index
        coefficient, _ = _safe_corr(x[keep], y[keep], method)
        if np.isfinite(coefficient):
            values.append(coefficient)
    if not values:
        return math.nan, math.nan
    return float(min(values)), float(max(values))


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    finite = values.dropna().sort_values()
    if finite.empty:
        return result
    count = len(finite)
    adjusted = finite.to_numpy() * count / np.arange(1, count + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result.loc[finite.index] = np.clip(adjusted, 0, 1)
    return result


def calculate_associations(
    frame: pd.DataFrame, bootstrap_samples: int, seed: int
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(seed)
    controls = [column for column in ["PEA_pct", "BDA_pct"] if column in frame]

    for kpfm_metric in KPFM_METRICS:
        for pl_metric in PL_METRICS:
            if kpfm_metric not in frame or pl_metric not in frame:
                continue
            pair = _valid_pair(frame, pl_metric, kpfm_metric)
            x = pair[pl_metric].to_numpy(dtype=float)
            y = pair[kpfm_metric].to_numpy(dtype=float)
            pearson_r, pearson_p = _safe_corr(x, y, "pearson")
            spearman_rho, spearman_p = _safe_corr(x, y, "spearman")
            pearson_low, pearson_high, pearson_boot_n = bootstrap_correlation_ci(
                x, y, "pearson", bootstrap_samples, rng
            )
            spearman_low, spearman_high, spearman_boot_n = bootstrap_correlation_ci(
                x, y, "spearman", bootstrap_samples, rng
            )
            partial_r, partial_p, partial_n, effective_controls = partial_correlation(
                frame, pl_metric, kpfm_metric, controls, "pearson"
            )
            partial_rho, partial_spearman_p, _, _ = partial_correlation(
                frame, pl_metric, kpfm_metric, controls, "spearman"
            )
            loo_low, loo_high = leave_one_out_range(x, y, "spearman")
            rows.append(
                {
                    "kpfm_metric": kpfm_metric,
                    "pl_metric": pl_metric,
                    "n": len(pair),
                    "pearson_r": pearson_r,
                    "pearson_p": pearson_p,
                    "pearson_ci_low": pearson_low,
                    "pearson_ci_high": pearson_high,
                    "pearson_bootstrap_valid": pearson_boot_n,
                    "spearman_rho": spearman_rho,
                    "spearman_p": spearman_p,
                    "spearman_ci_low": spearman_low,
                    "spearman_ci_high": spearman_high,
                    "spearman_bootstrap_valid": spearman_boot_n,
                    "spearman_loo_low": loo_low,
                    "spearman_loo_high": loo_high,
                    "composition_controls": ";".join(controls),
                    "effective_control_count": effective_controls,
                    "partial_n": partial_n,
                    "partial_pearson_r": partial_r,
                    "partial_pearson_p": partial_p,
                    "partial_spearman_rho": partial_rho,
                    "partial_spearman_p": partial_spearman_p,
                }
            )

    results = pd.DataFrame(rows)
    results["spearman_q_bh"] = benjamini_hochberg(results["spearman_p"])
    results["partial_spearman_q_bh"] = benjamini_hochberg(
        results["partial_spearman_p"]
    )
    results["association_is_exploratory"] = True
    results["causal_claim_supported"] = False
    return results


def plot_association_matrix(results: pd.DataFrame, output_base: Path) -> None:
    raw = results.pivot(index="kpfm_metric", columns="pl_metric", values="spearman_rho")
    adjusted = results.pivot(
        index="kpfm_metric", columns="pl_metric", values="partial_spearman_rho"
    )
    raw = raw.reindex(index=KPFM_METRICS, columns=PL_METRICS)
    adjusted = adjusted.reindex(index=KPFM_METRICS, columns=PL_METRICS)

    for matrix, suffix, title in zip(
        [raw, adjusted],
        ["raw", "composition_adjusted"],
        ["Raw rank associations", "Composition-adjusted rank associations"],
    ):
        fig, axis = plt.subplots(figsize=(10.8, 4.5), constrained_layout=True)
        image = axis.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1)
        axis.set_xticks(range(len(matrix.columns)))
        axis.set_xticklabels(
            [DISPLAY_LABELS.get(column, column) for column in matrix.columns],
            rotation=35,
            ha="right",
        )
        axis.set_yticks(range(len(matrix.index)))
        axis.set_yticklabels([DISPLAY_LABELS.get(row, row) for row in matrix.index])
        axis.set_title(title, loc="left", fontsize=12, fontweight="bold")
        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix.iat[row, column]
                if np.isfinite(value):
                    color = "white" if abs(value) >= 0.55 else "black"
                    axis.text(column, row, f"{value:.2f}", ha="center", va="center", color=color, fontsize=8)
        axis.set_xlabel("PL-derived measurement")
        axis.set_ylabel("photo-KPFM measurement")

        colorbar = fig.colorbar(image, ax=axis, shrink=0.82, pad=0.02)
        colorbar.set_label("Spearman coefficient")
        stem = output_base.parent / f"{output_base.name}_{suffix}"
        fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
        plt.close(fig)


def _add_linear_fit(axis: plt.Axes, x: np.ndarray, y: np.ndarray) -> None:
    if len(x) < 3 or np.std(x) == 0:
        return
    coefficients = np.polyfit(x, y, 1)
    x_fit = np.linspace(np.min(x), np.max(x), 100)
    axis.plot(x_fit, np.polyval(coefficients, x_fit), color="#B53A2D", lw=1.6)


def plot_selected_relationships(frame: pd.DataFrame, output_base: Path) -> None:
    panels = [
        ("Top_PL_peak_intensity", "photovoltage_signed_V"),
        ("dominant_peak_nm", "photovoltage_signed_V"),
        ("frac_n=1-like / 2D", "photovoltage_signed_V"),
        ("frac_3D / 3D-like", "photovoltage_signed_V"),
    ]
    label_offsets = [(5, 8), (5, -10), (-16, 0), (5, 0), (-14, 8), (-14, -10)]
    for index, (x_col, y_col) in enumerate(panels, start=1):
        fig, axis = plt.subplots(figsize=(5.4, 4.4), constrained_layout=True)
        subset = frame[["Well", x_col, y_col]].dropna()
        x = subset[x_col].to_numpy(dtype=float)
        y = subset[y_col].to_numpy(dtype=float)
        axis.scatter(x, y, s=45, color="#277DA1", edgecolor="white", linewidth=0.7, zorder=3)
        _add_linear_fit(axis, x, y)
        for row_number, (_, row) in enumerate(subset.iterrows()):
            offset = label_offsets[row_number % len(label_offsets)]
            axis.annotate(
                row["Well"],
                (row[x_col], row[y_col]),
                xytext=offset,
                textcoords="offset points",
                fontsize=7,
            )
        coefficient, p_value = _safe_corr(x, y, "spearman")
        axis.set_title(
            f"Spearman rho = {coefficient:.2f}; "
            f"p = {p_value:.3f}; n = {len(subset)}",
            loc="left",
            fontsize=10.5,
            fontweight="bold",
        )
        axis.set_xlabel(DISPLAY_LABELS.get(x_col, x_col))
        axis.set_ylabel("signed light-dark potential (V)")
        axis.margins(x=0.08, y=0.10)
        axis.grid(alpha=0.2)
        stem = output_base.parent / f"{output_base.name}_{index}"
        fig.savefig(stem.with_suffix(".png"), dpi=350, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
        plt.close(fig)


def classify_kpfm_observable(value: str) -> str:
    text = str(value).lower()
    if "photovoltage" in text or "light" in text and "dark" in text:
        return "Photovoltage / light-dark potential"
    if "work function" in text:
        return "Work function"
    if "contact potential" in text or "cpd" in text:
        return "Contact potential difference"
    if "surface potential" in text or "potential" in text:
        return "Surface potential"
    return "Other KPFM observable"


def classify_pl_observable(value: str) -> str:
    text = str(value).lower()
    if "lifetime" in text or "trpl" in text or "decay" in text:
        return "Carrier lifetime / TRPL"
    if "plqy" in text or "intensity" in text or "photoluminescence" in text:
        return "PL intensity / PLQY"
    if "peak" in text or "wavelength" in text or "phase" in text:
        return "PL peak / phase signature"
    return "Other PL observable"


def make_alignment_matrix(relationships_path: Path | None) -> pd.DataFrame:
    columns = [
        "literature_kpfm_class",
        "literature_pl_class",
        "relationship_count",
        "experiment_kpfm_variables",
        "experiment_pl_variables",
        "alignment_level",
        "interpretation_limit",
    ]
    if relationships_path is None or not relationships_path.exists():
        return pd.DataFrame(columns=columns)

    relationships = pd.read_csv(relationships_path)
    kpfm_column = (
        "kpfm_observable"
        if "kpfm_observable" in relationships
        else "kpfm_observable_and_direct_result"
    )
    pl_column = (
        "pl_observable"
        if "pl_observable" in relationships
        else "pl_family_observable_and_direct_result"
    )
    if kpfm_column not in relationships or pl_column not in relationships:
        raise ValueError("The relationship table does not contain recognized KPFM and PL observable columns.")
    relationships["literature_kpfm_class"] = relationships[kpfm_column].map(
        classify_kpfm_observable
    )
    relationships["literature_pl_class"] = relationships[pl_column].map(
        classify_pl_observable
    )
    counts = (
        relationships.groupby(["literature_kpfm_class", "literature_pl_class"])
        .size()
        .reset_index(name="relationship_count")
    )

    kpfm_map = {
        "Photovoltage / light-dark potential": (
            "photovoltage_signed_V; photovoltage_abs_V; mu_dark_V; mu_light_V",
            "direct",
        ),
        "Contact potential difference": (
            "mu_dark_V; mu_light_V; photovoltage_signed_V",
            "related",
        ),
        "Surface potential": ("mu_dark_V; mu_light_V", "direct"),
        "Work function": ("mu_dark_V; mu_light_V", "proxy_only"),
        "Other KPFM observable": ("not directly mapped", "unresolved"),
    }
    pl_map = {
        "PL intensity / PLQY": ("Top_PL_peak_intensity", "related_not_equivalent"),
        "Carrier lifetime / TRPL": ("not measured in matched table", "not_measured"),
        "PL peak / phase signature": (
            "dominant_peak_nm; PL-derived phase fractions",
            "direct",
        ),
        "Other PL observable": ("not directly mapped", "unresolved"),
    }

    rows = []
    for _, row in counts.iterrows():
        kpfm_variables, kpfm_level = kpfm_map[row["literature_kpfm_class"]]
        pl_variables, pl_level = pl_map[row["literature_pl_class"]]
        if "not_measured" in (kpfm_level, pl_level):
            alignment = "partial"
        elif "unresolved" in (kpfm_level, pl_level):
            alignment = "unresolved"
        elif "proxy_only" in (kpfm_level, pl_level):
            alignment = "proxy"
        elif "related_not_equivalent" in (kpfm_level, pl_level):
            alignment = "related"
        else:
            alignment = "direct"
        rows.append(
            {
                **row.to_dict(),
                "experiment_kpfm_variables": kpfm_variables,
                "experiment_pl_variables": pl_variables,
                "alignment_level": alignment,
                "interpretation_limit": (
                    "Co-measurement supports complementary interpretation only; "
                    "it does not establish a causal KPFM-to-PL pathway."
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def write_interpretation(
    output_path: Path,
    frame: pd.DataFrame,
    results: pd.DataFrame,
    source_paths: dict[str, str],
) -> None:
    strongest_raw = results.iloc[results["spearman_rho"].abs().argmax()]
    valid_partial = results.dropna(subset=["partial_spearman_rho"])
    strongest_partial = valid_partial.iloc[valid_partial["partial_spearman_rho"].abs().argmax()]
    significant_raw = int((results["spearman_q_bh"] < 0.05).sum())
    significant_partial = int((results["partial_spearman_q_bh"] < 0.05).sum())
    near_limit = int(frame["pl_near_upper_limit"].sum())

    text = f"""# PL/photo-KPFM Experimental Validation

## Scope

- Matched wells: **{len(frame)}**
- PL values near the chosen upper-limit flag (>=98,000 a.u.): **{near_limit}**
- Raw rank associations surviving BH q < 0.05: **{significant_raw}**
- Composition-adjusted rank associations surviving BH q < 0.05: **{significant_partial}**

This is a small, composition-designed plate subset. The statistics are
exploratory and do not establish causality. Dark potential, light potential,
and light-dark photovoltage are related quantities rather than independent
replicates. PL-derived phase fractions are also compositional and sum to one.

## Strongest Observed Associations

- Raw: `{strongest_raw['pl_metric']}` vs `{strongest_raw['kpfm_metric']}`,
  Spearman rho = **{strongest_raw['spearman_rho']:.3f}**,
  p = **{strongest_raw['spearman_p']:.3g}**,
  BH q = **{strongest_raw['spearman_q_bh']:.3g}**.
- After adjustment for PEA and BDA percentages (FA is the dependent remainder):
  `{strongest_partial['pl_metric']}` vs `{strongest_partial['kpfm_metric']}`,
  partial Spearman rho = **{strongest_partial['partial_spearman_rho']:.3f}**,
  p = **{strongest_partial['partial_spearman_p']:.3g}**,
  BH q = **{strongest_partial['partial_spearman_q_bh']:.3g}**.

These are the largest coefficients among the tested pairs, not confirmed
effects. Use the confidence intervals, adjusted q values, and leave-one-out
ranges in `experimental_kpfm_pl_associations.csv` when judging robustness.

## What the Literature Graph Validates

The literature pilot supports treating PL and KPFM as complementary views of
the same material/device state. PL reports emissive-state populations,
recombination, and phase signatures. KPFM reports local surface potential,
contact-potential/work-function contrast, and light-induced potential changes.
The evidence graph therefore supports joint interpretation, while preserving
alternative explanations such as composition, morphology, interfaces, and
measurement condition.

## Provenance

- Matched experimental CSV: `{source_paths['matched_csv']}`
- Phase source: `{source_paths['phase_source']}`
- Literature relationship table: `{source_paths.get('literature_relationships', 'not supplied')}`
"""
    output_path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    phase_table, phase_source = load_phase_table(args)
    experimental = prepare_experimental_table(args.matched_csv, phase_table)
    associations = calculate_associations(experimental, args.bootstrap, args.seed)
    alignment = make_alignment_matrix(args.literature_relationships)

    experimental_path = args.out_dir / "recovered_initial_phase_kpfm_table.csv"
    associations_path = args.out_dir / "experimental_kpfm_pl_associations.csv"
    alignment_path = args.out_dir / "literature_experiment_alignment_matrix.csv"
    experimental.to_csv(experimental_path, index=False)
    associations.to_csv(associations_path, index=False)
    alignment.to_csv(alignment_path, index=False)

    plot_association_matrix(
        associations, args.out_dir / "experimental_kpfm_pl_association_matrix"
    )
    plot_selected_relationships(
        experimental, args.out_dir / "experimental_selected_kpfm_pl_relationships"
    )

    source_paths = {
        "matched_csv": portable_path(args.matched_csv),
        "phase_source": phase_source,
    }
    if args.literature_relationships:
        source_paths["literature_relationships"] = portable_path(
            args.literature_relationships
        )
    write_interpretation(
        args.out_dir / "EXPERIMENTAL_VALIDATION_INTERPRETATION.md",
        experimental,
        associations,
        source_paths,
    )

    phase_error = float((experimental["phase_fraction_sum"] - 1).abs().max())
    summary = {
        "matched_wells": int(len(experimental)),
        "association_tests": int(len(associations)),
        "phase_fraction_max_sum_error": phase_error,
        "pl_values_ge_98000": int(experimental["pl_near_upper_limit"].sum()),
        "raw_spearman_bh_q_lt_0_05": int((associations["spearman_q_bh"] < 0.05).sum()),
        "partial_spearman_bh_q_lt_0_05": int(
            (associations["partial_spearman_q_bh"] < 0.05).sum()
        ),
        "causal_claim_supported": False,
        "source_paths": source_paths,
        "outputs": {
            "experimental_table": portable_path(experimental_path),
            "associations": portable_path(associations_path),
            "alignment_matrix": portable_path(alignment_path),
        },
    }
    (args.out_dir / "experimental_validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("Experimental PL/photo-KPFM validation complete")
    print(f"Matched wells: {len(experimental)}")
    print(f"Association tests: {len(associations)}")
    print(f"Outputs: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
