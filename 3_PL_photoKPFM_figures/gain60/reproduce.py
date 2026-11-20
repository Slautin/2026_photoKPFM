"""Recompute the gain-60 PL tables, matched statistics and six individual panels."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import types

import h5py
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
BOOTSTRAPS = 10000
SEED = 20260816
COMPOSITIONS = ["nominal_FA_pct", "nominal_BDA_pct", "nominal_PEA_pct"]
PANELS = ["4a", "4b", "S4", "5a", "5b", "S5"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def load_sources():
    provenance = json.loads((ROOT / "provenance.json").read_text(encoding="utf-8"))
    for source in provenance["sources"]:
        path = (ROOT / source["path"]).resolve()
        if not path.is_relative_to(PROJECT.resolve()) or sha(path) != source["sha256"]:
            raise ValueError("Source identity differs: " + source["path"])
    with np.load(ROOT / "data/top_spectra_raw.npz", allow_pickle=False) as source:
        raw = {key: source[key].copy() for key in source.files}
    mapping = pd.read_csv(ROOT / "data/well_map.csv", keep_default_na=False)
    metadata = json.loads((ROOT / "data/measurement_metadata.json").read_text(encoding="utf-8"))
    legacy = pd.read_csv(PROJECT / "data/processed/matched_pl_photokpfm_metrics.csv").set_index("well_id")
    phase = pd.read_csv(PROJECT / "data/processed/initial_phase_kpfm_table.csv").set_index("Well")
    return raw, mapping, metadata, legacy, phase


def peak_arrays(raw):
    values = raw["numeric_values"]
    waves = raw["wavelengths_nm"]
    mask = (waves >= 455) & (waves <= 850)
    valid_cells = np.isfinite(values[:, :, mask]) & (raw["cell_type_codes"][:, :, mask] < 3)
    valid = valid_cells.all(axis=2)
    if not valid.all():
        raise ValueError("This release requires every source band channel to be valid; invalid spectra are never imputed or silently omitted")
    peaks = values[:, :, mask].max(axis=2)
    denominators = peaks.max(axis=0)
    if np.any(denominators == 0):
        raise ValueError("A zero series maximum gives an undefined normalization")
    normalized = peaks / denominators[None, :]
    peak_wavelengths = waves[mask][values[:, :, mask].argmax(axis=2)]
    return peaks, normalized, peak_wavelengths


def gaussian_checks(legacy, phase):
    path = PROJECT / "scripts/experimental/data_processing.py"
    archived = types.ModuleType("archived_gain100_measurement_helpers")
    archived.__file__ = str(path)
    exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), archived.__dict__)
    archived.KPFM_H5 = PROJECT / "data/raw/photokpfm_measurements.h5"
    acquisitions = archived.load_kpfm_acquisitions()
    grouped = acquisitions.groupby("well_id")
    figure5 = grouped.delta_SP.mean().abs().reindex(legacy.index)
    s5 = grouped.abs_delta_SP.mean().round(6).reindex(phase.index)
    errors = (grouped.delta_SP_se.agg(lambda x: np.sqrt(np.sum(np.square(x)))) / grouped.size()).reindex(legacy.index)
    difference = float(np.max(np.abs(figure5 - legacy.abs_delta_SP)))
    error_difference = float(np.max(np.abs(errors - legacy.delta_SP_se)))
    if difference >= 1e-12 or error_difference >= 1e-12 or not np.array_equal(s5, phase.photovoltage_abs_V):
        raise ValueError("Archived Gaussian values or fit-center errors do not reproduce from the raw scan pairs")
    return {"raw_scan_pairs": len(acquisitions), "nominal_compositions": len(grouped), "figure5_max_abs_difference_V": difference,
            "s5_exact_after_six_decimal_rounding": True, "s5_fit_error_max_abs_difference_V": error_difference}


def validate_sources():
    """Validate public inputs without creating output files or rendering figures."""
    raw, mapping, metadata, legacy, phase = load_sources()
    values = raw["numeric_values"]
    wells = list(raw["wells"].astype(str))
    expected_wells = [f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)]
    if values.shape != (22, 96, 416) or wells != expected_wells:
        raise ValueError("Unexpected acquisition/well/channel array")
    if raw["raw_cells"].shape != values.shape or raw["cell_type_codes"].shape != values.shape:
        raise ValueError("Raw cell and type arrays must retain every measured channel")
    if not np.array_equal(raw["wavelengths_nm"], np.arange(455, 871)) or not np.array_equal(raw["raw_read_numbers"], np.arange(1, 65, 3)):
        raise ValueError("Observed wavelength or top-read indices differ")
    if not np.array_equal(raw["acquisition_numbers"], np.arange(1, 23)) or not np.array_equal(raw["source_well_excel_columns"], np.arange(3, 99)):
        raise ValueError("Source acquisition/column indices differ")
    if list(mapping.well_id) != wells or mapping.well_id.nunique() != 96 or not np.allclose(mapping[COMPOSITIONS].sum(axis=1), 100):
        raise ValueError("Invalid per-well composition map")
    if len(mapping[COMPOSITIONS].drop_duplicates()) != 91:
        raise ValueError("Expected 91 nominal recipes across 96 unpooled wells")
    if str(raw["input_xlsx_sha256"]) != metadata["source_workbook_sha256"] or metadata["gain"] != 60:
        raise ValueError("Measurement provenance differs")
    numeric = raw["cell_type_codes"] < 3
    if not np.array_equal(raw["raw_cells"][numeric].astype(float), values[numeric]):
        raise ValueError("Literal numeric cells differ from their source values")
    peaks, normalized, _ = peak_arrays(raw)
    if not np.array_equal(peaks, raw["peak_max_counts"]):
        raise ValueError("Stored peak cache differs from a direct observed maximum")
    if not np.array_equal(normalized, raw["normalized_peak_series"]):
        raise ValueError("Stored normalization cache differs from the complete series")
    if not np.array_equal(raw["peak_valid"], np.ones((22, 96), dtype=bool)):
        raise ValueError("Cached validity flags differ")
    with h5py.File(PROJECT / "data/raw/photokpfm_measurements.h5", "r") as handle:
        ids = np.asarray(handle["idx"], dtype=int).ravel()
        if handle["dark_data"].shape != (19, 128, 128) or handle["light_data"].shape != (19, 128, 128):
            raise ValueError("Unexpected raw photoKPFM dimensions")
    afm_wells = [f'{"ABCDEFGH"[i % 8]}{i // 8 + 1}' for i in ids]
    if len(ids) != 19 or set(afm_wells) != set(legacy.index) or set(legacy.index) != set(phase.index) or len(legacy) != 13:
        raise ValueError("Nominal PL/photoKPFM crosswalk differs")
    for well in legacy.index:
        source = mapping.set_index("well_id").loc[well, COMPOSITIONS].to_numpy(float)
        historical = phase.loc[well, ["FA_pct", "BDA_pct", "PEA_pct"]].to_numpy(float)
        if not np.allclose(source, historical, rtol=0, atol=1e-9):
            raise ValueError("Matched nominal composition differs: " + well)
    gaussian = gaussian_checks(legacy, phase)
    return {"status": "passed", "top_acquisitions": 22, "wells": 96, "distinct_nominal_compositions": 91,
            "raw_channels": 416, "analysis_channels": 396, "valid_peak_measurements": int(peaks.size),
            "normalization_eligible_wells": 96, "comparison_acquisition": 22, "comparison_raw_read": 64,
            "photoKPFM_nominal_matches": 13, "photoKPFM_scan_pairs": 19, "gaussian_reference_check": gaussian,
            "literal_source_cells_match": True, "source_hashes_match": True, "outputs_created": False}


def partial_spearman(x, y, controls):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    ranked = np.column_stack([stats.rankdata(controls[:, i]) for i in range(controls.shape[1])])
    design = np.column_stack([np.ones(len(x)), ranked])
    k = int(np.linalg.matrix_rank(design)) - 1
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    coefficient = float(stats.pearsonr(ex, ey).statistic)
    df = len(x) - k - 2
    statistic = coefficient * math.sqrt(df / (1-coefficient**2)) if abs(coefficient) < 1 else math.inf
    return coefficient, float(2*stats.t.sf(abs(statistic), df)), k, df


def association(x, y, controls, indices):
    pearson, spearman = stats.pearsonr(x, y), stats.spearmanr(x, y)
    coefficient, probability, k, df = partial_spearman(x, y, controls)
    bx, by = x[indices], y[indices]
    dx, dy = bx-bx.mean(axis=1, keepdims=True), by-by.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(dx*dx, axis=1)*np.sum(dy*dy, axis=1))
    eligible = denominator > 0
    correlations = np.sum(dx[eligible]*dy[eligible], axis=1)/denominator[eligible]
    correlations = correlations[np.isfinite(correlations)]
    low, high = np.percentile(correlations, [2.5, 97.5])
    return {"n": len(x), "pearson_r": float(pearson.statistic), "pearson_p": float(pearson.pvalue),
            "pearson_ci_low": float(low), "pearson_ci_high": float(high), "pearson_bootstrap_valid": len(correlations),
            "spearman_rho": float(spearman.statistic), "spearman_p": float(spearman.pvalue),
            "partial_spearman_rho": coefficient, "partial_spearman_p": probability,
            "effective_control_count": k, "partial_degrees_of_freedom": df, "new_q_values_reported": False}


def derive_tables():
    raw, mapping, metadata, legacy, phase = load_sources()
    peaks, normalized, wavelengths = peak_arrays(raw)
    final = mapping.copy()
    final["acquisition_number"], final["read_block_number"], final["gain"] = 22, 64, 60
    final["acquisition22_status"] = "valid_observed_band_max"
    final["acquisition22_valid_observed_max_counts"] = peaks[-1]
    final["acquisition22_peak_wavelength_nm_first"] = wavelengths[-1]
    final["all22_maxima_valid"], final["valid_acquisition_count"], final["retention_eligible"] = True, 22, True
    final["retention_denominator_max_counts"] = peaks.max(axis=0)
    final["retention_I22_over_series_max"] = normalized[-1]
    series = []
    for wi, row in mapping.iterrows():
        for ai in range(22):
            series.append({"well_id": row.well_id, "acquisition_number": ai+1, "read_block_number": int(raw["raw_read_numbers"][ai]),
                           "gain": 60, "peak_max_counts": float(peaks[ai, wi]), "normalized_peak_over_all22_max": float(normalized[ai, wi])})
    indexed = final.set_index("well_id")
    with h5py.File(PROJECT / "data/raw/photokpfm_measurements.h5", "r") as handle:
        ids = np.asarray(handle["idx"], dtype=int).ravel()
        y = np.asarray(handle["y_train"], dtype=float).ravel()
    wells = np.asarray([f'{"ABCDEFGH"[i % 8]}{i // 8 + 1}' for i in ids])
    matched = []
    for well in legacy.index:
        source = indexed.loc[well]
        matched.append({"well_id": well, **{key: float(source[key]) for key in COMPOSITIONS},
                        "PL_raw_max_counts": float(source.acquisition22_valid_observed_max_counts),
                        "AFM_Figure5_abs_mean_Gaussian_V": float(legacy.loc[well].abs_delta_SP),
                        "AFM_S5_mean_abs_Gaussian_V": float(phase.loc[well].photovoltage_abs_V),
                        "AFM_archived_delta_SP_se_V": float(legacy.loc[well].delta_SP_se),
                        "AFM_retained_y_abs_mean_signed_V": abs(float(np.mean(y[wells == well]))),
                        "n_AFM_acquisitions": int(np.sum(wells == well))})
    matched = pd.DataFrame(matched)
    matched["PL_norm_to_matched13_max"] = matched.PL_raw_max_counts / matched.PL_raw_max_counts.max()
    rank_method = "average" if matched.PL_raw_max_counts.duplicated().any() or matched.AFM_Figure5_abs_mean_Gaussian_V.duplicated().any() else "first"
    matched["PL_rank"] = matched.PL_raw_max_counts.rank(ascending=False, method=rank_method)
    matched["AFM_rank"] = matched.AFM_Figure5_abs_mean_Gaussian_V.rank(ascending=False, method=rank_method)
    controls = matched[["nominal_PEA_pct", "nominal_BDA_pct"]].to_numpy(float)
    indices = np.random.default_rng(SEED).integers(0, len(matched), size=(BOOTSTRAPS, len(matched)))
    x = matched.PL_raw_max_counts.to_numpy(float)
    estimates = {name: association(x, matched[column].to_numpy(float), controls, indices) for name, column in {
        "S5_mean_absolute_Gaussian_rounded6": "AFM_S5_mean_abs_Gaussian_V",
        "Figure5_absolute_mean_signed_Gaussian": "AFM_Figure5_abs_mean_Gaussian_V",
        "retained_y_sensitivity_only": "AFM_retained_y_abs_mean_signed_V"}.items()}
    statistics = {"primary": estimates["S5_mean_absolute_Gaussian_rounded6"], "descriptor_sensitivities": estimates,
                  "bootstrap": {"resamples": BOOTSTRAPS, "seed": SEED, "rng": "NumPy default_rng PCG64", "method": "paired nominal-composition-row percentile bootstrap for Pearson r", "interval_percentiles": [2.5, 97.5]},
                  "partial_spearman": "Rank x, y, PEA and BDA; residualize on intercept and ranked controls; correlate residuals; two-sided t test with n-k-2 degrees of freedom",
                  "rank_method": rank_method, "rank_category_rules": "Top four means rank <= 4; rank lead means difference >= 5. Display categories are not significance tests.",
                  "cross_library_scope": "13 nominal-composition matches across independently prepared PL and photoKPFM libraries; 19 photoKPFM scan pairs include within-droplet repeated locations",
                  "multiple_testing": "No new 36-test family or q values: new spectral-component fits were not produced. Historical gain-100 q values are not gain-60 statistics.",
                  "measurement_uncertainty": "FigureS5 preserves archived delta_SP_se from the Gaussian fit-center variance aggregation; no new error bars are inferred from PL counts.",
                  "causal_claim_supported": False}
    return final, pd.DataFrame(series), matched, statistics


def reproduce(output):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("Choose a fresh output directory")
    validation = validate_sources()
    final, series, matched, statistics = derive_tables()
    (output / "data").mkdir(parents=True)
    final.to_csv(output / "data/final_acquisition_and_retention_by_well.csv", index=False)
    series.to_csv(output / "data/normalized_peak_series.csv", index=False)
    matched.to_csv(output / "data/matched13_panel_source.csv", index=False)
    write_json(output / "data/statistics.json", statistics)
    source = ROOT / "render.py"
    renderer = types.ModuleType("gain60_render")
    renderer.__file__ = str(source)
    exec(compile(source.read_text(encoding="utf-8"), str(source), "exec"), renderer.__dict__)
    figures = renderer.render_all(final, series, matched, output / "figures")
    validation["outputs_created"] = True
    validation["figure_files"] = [path.relative_to(output).as_posix() for path in figures]
    validation["figure_count"] = len(figures)
    write_json(output / "validation.json", validation)
    print(json.dumps({"status": "passed", "panels": 6, "figure_files": len(figures), "valid_peak_measurements": 2112, "nominal_matches": 13}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.validate_only:
        print(json.dumps(validate_sources(), indent=2))
    elif args.output is None:
        parser.error("--output <fresh-directory> is required unless --validate-only is selected")
    else:
        reproduce(args.output)


if __name__ == "__main__":
    main()
