from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
import math
import os
import re
import runpy
import tempfile
import tokenize
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import h5py
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_NAME = "MANIFEST.csv"
DUPLICATE_AUDIT_NAME = "DUPLICATE_AUDIT.csv"
REPORT_NAME = "VALIDATION_REPORT.json"
EXCLUDED_DIRECTORIES = {".git", ".venv", "__pycache__", "reproduced_results"}
TRACKING_FILES = {MANIFEST_NAME, DUPLICATE_AUDIT_NAME, REPORT_NAME}
CANONICAL_TABLE_PATH = "results/tables/supplementary_tables_S1-S4.docx"
CANONICAL_TABLE_PROVENANCE_PATH = (
    "results/tables/supplementary_tables_S1-S4.provenance.json"
)
CANONICAL_TABLE_GENERATOR = "scripts/publication/generate_supplementary_tables.py"
LEGACY_TABLE_PATH = "results/tables/supplementary_tables_scientifically_corrected.docx"
TERNARY_GENERATOR = "scripts/publication/generate_ternary_photokpfm_panels.py"
TERNARY_SOURCE_DIRECTORY = "results/source_data/ternary_photokpfm"
TERNARY_SOURCE_FILES = (
    f"{TERNARY_SOURCE_DIRECTORY}/photoKPFM_acquisitions_all_19.csv",
    f"{TERNARY_SOURCE_DIRECTORY}/Figure3_and_FigureS2_GP_grid_source_data.csv",
    f"{TERNARY_SOURCE_DIRECTORY}/Figure3_representative_histograms_source_data.csv",
    f"{TERNARY_SOURCE_DIRECTORY}/FigureS1_repeat_location_variability_source_data.csv",
    f"{TERNARY_SOURCE_DIRECTORY}/ternary_GP_model_and_figure_provenance.json",
    f"{TERNARY_SOURCE_DIRECTORY}/ternary_photokpfm_output_manifest.csv",
)
AUTHOR_SUPPLIED_ASSETS = {
    "figure1": "results/figures/author_supplied/Figure1_workflow_author_supplied.png",
    "figure2": "results/figures/author_supplied/Figure2_binary_photokpfm_author_supplied.jpeg",
}
STALE_FIGURE_STEMS = (
    "results/figures/main/Figure4A_PL_189min_ternary",
    "results/figures/main/Figure4B_PL_retention_ternary",
    "results/figures/main/Figure5A_matched_screening_space",
    "results/figures/main/Figure5B_candidate_rank_comparison",
    "results/figures/supplementary/FigureS4_representative_PL_time_traces",
    "results/figures/main/Figure4C_representative_PL_time_traces",
    "results/figures/main/figure_photokpfm_knowledge_graph_corrected",
    "results/figures/supplementary/FigureS_all_13_PL_time_traces",
    "results/figures/supplementary/FigureS_PL_photoKPFM_association",
    "results/figures/supplementary/FigureS_plate_composition_map",
    "results/figures/supplementary/figure_literature_measurement_pairings_corrected",
)
STALE_RELEASE_PATHS = (
    LEGACY_TABLE_PATH,
    "results/tables/supplementary_tables_current_manuscript_S1-S4.docx",
    "results/tables/supplementary_tables_current_manuscript_S1-S4.provenance.json",
    "_publication_refresh",
)
EXPECTED_WELLS = [f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)]
EXPECTED_GRAPH_NODE_TYPES = {
    "EvidenceRecord": 15,
    "KPFMObservation": 15,
    "Paper": 7,
    "ProposedLiteratureInterpretation": 15,
    "SampleComparison": 15,
    "ScientificClaim": 15,
    "SteadyStatePLObservation": 8,
    "TRPLContextObservation": 7,
}
EXPECTED_GRAPH_EDGE_TYPES = {
    "COMPLEMENTS": 15,
    "HAS_KPFM_OBSERVATION": 15,
    "HAS_OPTICAL_OBSERVATION": 15,
    "INTERPRETED_AS": 15,
    "REPORTS_COMPARISON": 15,
    "SUPPORTS": 15,
    "SUPPORTS_CLAIM": 30,
}
REQUIRED_MANUSCRIPT_IDS = {
    "figure1",
    "figure2",
    "figure3",
    "figure4a",
    "figure4b",
    "figure5a",
    "figure5b",
    "figures1",
    "figures2",
    "figures3",
    "figures4",
    "figures5",
    "figures6",
    "figures7",
    "tables1",
    "tables2",
    "tables3",
    "tables4",
}
ALLOWED_ASSET_STATUSES = {
    "reproducible", "author-supplied-with-explicit-scope", "historical-gain100"
}
GAIN60_GENERATOR = "gain60/reproduce.py"
GAIN60_S7_DIRECTORY = "gain60/figure_s7"
GAIN60_S7_STEM = f"{GAIN60_S7_DIRECTORY}/FigureS7_gain60"
GAIN60_S7_GENERATOR = f"{GAIN60_S7_DIRECTORY}/FigureS7_plot.py"
GAIN60_S7_SOURCES = (
    "gain60/data/top_spectra_raw.npz",
    "gain60/data/well_map.csv",
    f"{GAIN60_S7_DIRECTORY}/FigureS7_source_spectra.csv",
)
GAIN60_ASSETS = {
    'figure4a': (
        'gain60/figures/Figure4a_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'gain60/data/final_acquisition_and_retention_by_well.csv',
            'gain60/data/matched13_panel_source.csv',
        ),
    ),
    'figure4b': (
        'gain60/figures/Figure4b_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'gain60/data/final_acquisition_and_retention_by_well.csv',
            'gain60/data/matched13_panel_source.csv',
        ),
    ),
    'figure5a': (
        'gain60/figures/Figure5a_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'data/raw/photokpfm_measurements.h5',
            'data/processed/initial_phase_kpfm_table.csv',
            'data/processed/matched_pl_photokpfm_metrics.csv',
            'gain60/data/matched13_panel_source.csv',
        ),
    ),
    'figure5b': (
        'gain60/figures/Figure5b_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'data/raw/photokpfm_measurements.h5',
            'data/processed/initial_phase_kpfm_table.csv',
            'data/processed/matched_pl_photokpfm_metrics.csv',
            'gain60/data/matched13_panel_source.csv',
        ),
    ),
    'figures4': (
        'gain60/figures/FigureS4_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'gain60/data/normalized_peak_series.csv',
        ),
    ),
    'figures5': (
        'gain60/figures/FigureS5_gain60',
        (
            'gain60/data/top_spectra_raw.npz',
            'gain60/data/well_map.csv',
            'data/raw/photokpfm_measurements.h5',
            'data/processed/initial_phase_kpfm_table.csv',
            'data/processed/matched_pl_photokpfm_metrics.csv',
            'gain60/data/matched13_panel_source.csv',
            'gain60/data/statistics.json',
        ),
    ),
}
RELEASE_FIGURES = {
    'results/figures/main/Figure3_GP_mean_map': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/main/Figure3_representative_SP_histogram_1': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/main/Figure3_representative_SP_histogram_2': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/main/Figure3_representative_SP_histogram_3': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/main/Figure3_representative_SP_histogram_4': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/supplementary/FigureS1_repeat_location_variability': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/supplementary/FigureS2_GP_uncertainty_map': ((TERNARY_GENERATOR,), ('data/raw/photokpfm_measurements.h5',)),
    'results/figures/supplementary/FigureS3_plate_composition_map': (('scripts/publication/generate_publication_figures.py',), ('gain60/data/well_map.csv',)),
    'results/figures/supplementary/FigureS6_evidence_graph': (('scripts/literature/build_evidence_graph.py', 'scripts/literature/plot_evidence_graph.py'), ('data/literature/evidence_graph_nodes.csv', 'data/literature/evidence_graph_edges.csv', 'data/literature/table_s4_retained_relationships.csv')),
}
RELEASE_FIGURES.update({
    stem: ((GAIN60_GENERATOR,), sources) for stem, sources in GAIN60_ASSETS.values()
})
RELEASE_FIGURES[GAIN60_S7_STEM] = ((GAIN60_S7_GENERATOR,), GAIN60_S7_SOURCES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate the PL/photoKPFM release package.")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--refresh-manifest", action="store_true")
    parser.add_argument("--write-report", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def included_path(path: Path, root: Path) -> bool:
    relative = path.relative_to(root)
    return not any(part in EXCLUDED_DIRECTORIES for part in relative.parts)


def actual_comments(source: str) -> list[tuple[int, str]]:
    comments = []
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            comments.append((token.start[0], token.string))
    return comments


def allowed_comment(line_number: int, comment: str) -> bool:
    lowered = comment.lower()
    return (line_number == 1 and comment.startswith("#!")) or any(
        marker in lowered for marker in ("# noqa", "# type: ignore", "# pragma", "# nosec")
    )


def role_for(relative_path: str) -> str:
    if relative_path.startswith(f"{GAIN60_S7_DIRECTORY}/"):
        suffix = Path(relative_path).suffix.lower()
        if suffix in {".png", ".svg", ".pdf"}:
            return "figure"
        if suffix == ".csv":
            return "figure_source_data"
        if suffix == ".txt":
            return "validation_data" if relative_path.endswith("_validation.txt") else "documentation"
    if relative_path == "gain60/data/top_spectra_raw.npz":
        return "raw_data"
    if relative_path.startswith("gain60/data/"):
        return "processed_data"
    if relative_path.startswith("gain60/figures/"):
        return "figure"
    if relative_path.startswith("gain60/") and relative_path.endswith(".py"):
        return "analysis_code"
    if relative_path.startswith("data/raw/"):
        return "raw_data"
    if relative_path.startswith("data/processed/"):
        return "processed_data"
    if relative_path.startswith("data/validation/"):
        return "validation_data"
    if relative_path.startswith("data/literature/"):
        return "literature_evidence"
    if relative_path.startswith("scripts/"):
        return "analysis_code"
    if relative_path.startswith("notebooks/"):
        return "analysis_notebook"
    if relative_path.startswith("results/figures/"):
        return "figure"
    if relative_path.startswith("results/source_data/"):
        return "figure_source_data"
    if relative_path.startswith("results/tables/"):
        return "supplementary_table"
    if relative_path.startswith("docs/"):
        return "documentation"
    return "repository_metadata"


def validate_required_files(root: Path, refresh_manifest: bool) -> dict[str, object]:
    required = [
        root / "README.md",
        root / "requirements.txt",
        root / "docs" / "REPRODUCIBILITY_GUIDE.md",
        root / "docs" / "FIGURE_PROVENANCE.md",
        root / "docs" / "MANUSCRIPT_ASSET_INDEX.csv",
        root / "notebooks" / "README.md",
        root / "scripts" / "reproduce_all.py",
        root / GAIN60_GENERATOR,
        root / "scripts" / "publication" / "generate_publication_figures.py",
        root / TERNARY_GENERATOR,
        root / "scripts" / "literature" / "build_evidence_graph.py",
        root / "scripts" / "literature" / "plot_evidence_graph.py",
        root / "data" / "raw" / "pl_plate_reader_export.csv",
        root / "data" / "raw" / "photokpfm_measurements.h5",
        root / "data" / "literature" / "table_s3_relationship_summary.csv",
        root / "data" / "literature" / "table_s4_retained_relationships.csv",
        root / CANONICAL_TABLE_PATH,
        root / CANONICAL_TABLE_PROVENANCE_PATH,
        root / CANONICAL_TABLE_GENERATOR,
        *(root / path for path in AUTHOR_SUPPLIED_ASSETS.values()),
        *(root / path for path in TERNARY_SOURCE_FILES),
    ]
    if not refresh_manifest:
        required.append(root / MANIFEST_NAME)
    missing = [path.relative_to(root).as_posix() for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Required release files are missing: {missing}")
    return {"required_files": len(required), "missing": 0}


def validate_relationships(root: Path) -> dict[str, object]:
    relationships = pd.read_csv(root / "data" / "literature" / "table_s4_retained_relationships.csv")
    summary = pd.read_csv(root / "data" / "literature" / "table_s3_relationship_summary.csv")
    if len(relationships) != 15:
        raise ValueError(f"Expected 15 corrected relationships, found {len(relationships)}.")
    if relationships["relationship_id"].duplicated().any():
        raise ValueError("Duplicate relationship identifiers found in Table S4 export.")
    if "R08" in set(relationships["relationship_id"].astype(str)):
        raise ValueError("Removed relationship R08 is still present.")
    if int(pd.to_numeric(summary["retained_relationships"], errors="raise").sum()) != 15:
        raise ValueError("Table S3 retained counts do not sum to 15.")
    statuses = sorted(relationships["causal_status"].dropna().astype(str).unique().tolist())
    if statuses != ["complementary_only"]:
        raise ValueError(f"Unexpected causal status values: {statuses}")
    return {
        "retained_relationships": int(len(relationships)),
        "unique_papers": int(relationships["doi"].nunique()),
        "causal_status": statuses,
    }


def numeric_signal(value: str) -> bool:
    text = str(value).strip()
    if not text or text.upper() in {"OVRFLW", "OVERFLOW", "OVER", "NAN"}:
        return False
    try:
        return math.isfinite(float(text))
    except ValueError:
        return False


def validate_raw_data(root: Path) -> dict[str, object]:
    pl_path = root / "data" / "raw" / "pl_plate_reader_export.csv"
    with pl_path.open(newline="", errors="replace") as handle:
        rows = list(csv.reader(handle))
    starts = []
    for index, row in enumerate(rows):
        first = row[0].strip() if row else ""
        match = re.match(r"Read\s+(\d+):(.*)", first)
        if match:
            starts.append((index, int(match.group(1))))
    if not starts:
        raise ValueError("No Cytation5 read blocks found in the raw PL CSV.")
    valid_top_reads = []
    for block_index, (start, read_number) in enumerate(starts):
        if (read_number - 1) % 3:
            continue
        end = starts[block_index + 1][0] if block_index + 1 < len(starts) else len(rows)
        header_index = next(
            (
                line_number
                for line_number in range(start, end)
                if len(rows[line_number]) > 1 and rows[line_number][1].strip() == "Wavelength"
            ),
            None,
        )
        if header_index is None:
            raise ValueError(f"Missing wavelength header in top-read block {read_number}.")
        wells = [cell.strip() for cell in rows[header_index][2:] if cell.strip()]
        if wells != EXPECTED_WELLS:
            raise ValueError(f"Unexpected well header in top-read block {read_number}.")
        observed = [False] * len(EXPECTED_WELLS)
        for row in rows[header_index + 1 : end]:
            if len(row) < 2 or not numeric_signal(row[1]):
                continue
            wavelength = float(row[1])
            if not 450.0 <= wavelength <= 850.0:
                continue
            for column, value in enumerate(row[2:98]):
                if numeric_signal(value):
                    observed[column] = True
        if any(observed):
            if not all(observed):
                missing_wells = [well for well, present in zip(EXPECTED_WELLS, observed) if not present]
                raise ValueError(f"PL top-read block {read_number} is incomplete: {missing_wells}")
            valid_top_reads.append(read_number)
    expected_read_numbers = list(range(1, 65, 3))
    if valid_top_reads != expected_read_numbers:
        raise ValueError(
            f"Expected 22 complete PL top reads {expected_read_numbers}, found {valid_top_reads}."
        )
    h5_path = root / "data" / "raw" / "photokpfm_measurements.h5"
    with h5py.File(h5_path, "r") as handle:
        required = {"X", "X_train", "coord_array", "dark_data", "dark_fn", "idx", "light_data", "light_fn", "y_train"}
        missing = sorted(required.difference(handle.keys()))
        if missing:
            raise ValueError(f"Raw photoKPFM HDF5 datasets are missing: {missing}")
        idx = [int(value) for value in handle["idx"][:]]
        shapes = {name: tuple(int(value) for value in handle[name].shape) for name in required}
    expected_shapes = {
        "X": (91, 2),
        "X_train": (19, 2),
        "coord_array": (96, 2),
        "dark_data": (19, 128, 128),
        "dark_fn": (19,),
        "idx": (19,),
        "light_data": (19, 128, 128),
        "light_fn": (19,),
        "y_train": (19,),
    }
    if shapes != expected_shapes:
        raise ValueError(f"Unexpected raw photoKPFM dataset shapes: {shapes}")
    if len(idx) != 19 or len(set(idx)) != 13 or not all(0 <= value < 96 for value in idx):
        raise ValueError("Expected 19 valid photoKPFM acquisition pairs spanning 13 unique wells.")
    return {
        "pl_dataset_scope": "historical_gain100",
        "pl_wells": 96,
        "pl_top_timepoints": len(valid_top_reads),
        "pl_parsed_rows": 96 * len(valid_top_reads),
        "photokpfm_acquisition_pairs": len(idx),
        "photokpfm_unique_wells": len(set(idx)),
        "pl_sha256": sha256(pl_path),
        "photokpfm_sha256": sha256(h5_path),
    }


def validate_processed_data(root: Path) -> dict[str, object]:
    processed = root / "data" / "processed"
    specifications = {
        "pl_retention_metrics.csv": (96, "well_id", 96),
        "plate_composition_map.csv": (96, "Well", 96),
        "matched_pl_photokpfm_metrics.csv": (13, "well_id", 13),
        "candidate_rank_comparison.csv": (13, "well_id", 13),
        "representative_pl_trace_selection.csv": (13, "well_id", 13),
        "representative_pl_traceability.csv": (5, "well_id", 5),
        "representative_pl_traces.csv": (110, "well_id", 5),
        "all_matched_pl_traces.csv": (286, "well_id", 13),
        "initial_phase_kpfm_table.csv": (13, "Well", 13),
        "exploratory_association_statistics.csv": (36, None, None),
        "photokpfm_acquisition_log.csv": (19, "Well ID", 13),
    }
    frames = {}
    inventory = {}
    for name, (expected_rows, identifier, expected_unique) in specifications.items():
        path = processed / name
        if not path.is_file():
            raise FileNotFoundError(f"Required processed dataset is missing: data/processed/{name}")
        frame = pd.read_csv(path)
        frames[name] = frame
        if len(frame) != expected_rows:
            raise ValueError(f"{name} has {len(frame)} rows; expected {expected_rows}.")
        unique = None
        if identifier:
            if identifier not in frame.columns:
                raise ValueError(f"{name} is missing identifier column {identifier}.")
            unique = int(frame[identifier].astype(str).nunique())
            if unique != expected_unique:
                raise ValueError(f"{name} has {unique} unique {identifier} values; expected {expected_unique}.")
        inventory[name] = {"rows": len(frame), "unique_identifiers": unique, "sha256": sha256(path)}
    composition_columns = ["FAPbI3_pct", "BDAPbI4_pct", "PEA2PbI4_pct"]
    retention = frames["pl_retention_metrics.csv"]
    if retention[composition_columns].drop_duplicates().shape[0] != 91:
        raise ValueError("Expected 91 unique nominal composition triplets across 96 PL wells.")
    selection = frames["representative_pl_trace_selection.csv"]
    selected_mask = selection["selected_for_representative_panel"].astype(str).str.lower().eq("true")
    selected_wells = set(selection.loc[selected_mask, "well_id"].astype(str))
    if len(selected_wells) != 5:
        raise ValueError(f"Expected five representative PL wells, found {sorted(selected_wells)}.")
    traceability_wells = set(frames["representative_pl_traceability.csv"]["well_id"].astype(str))
    representative_wells = set(frames["representative_pl_traces.csv"]["well_id"].astype(str))
    if selected_wells != traceability_wells or selected_wells != representative_wells:
        raise ValueError("Representative PL selection, traceability, and trace tables disagree.")
    representative_counts = frames["representative_pl_traces.csv"].groupby("well_id").size()
    all_trace_counts = frames["all_matched_pl_traces.csv"].groupby("well_id").size()
    if not representative_counts.eq(22).all() or not all_trace_counts.eq(22).all():
        raise ValueError("Each supplied PL trace must contain exactly 22 timepoints.")
    with h5py.File(root / "data" / "raw" / "photokpfm_measurements.h5", "r") as handle:
        h5_indices = [int(index) for index in handle["idx"][:]]
        h5_wells = {
            f"{'ABCDEFGH'[index % 8]}{index // 8 + 1}" for index in h5_indices
        }
    matched_wells = set(frames["matched_pl_photokpfm_metrics.csv"]["well_id"].astype(str))
    candidate_wells = set(frames["candidate_rank_comparison.csv"]["well_id"].astype(str))
    phase_wells = set(frames["initial_phase_kpfm_table.csv"]["Well"].astype(str))
    if matched_wells != h5_wells or candidate_wells != h5_wells or phase_wells != h5_wells:
        raise ValueError("Processed matched, rank, phase, and raw photoKPFM well sets disagree.")
    acquisition_log = frames["photokpfm_acquisition_log.csv"]
    acquisition_indices = pd.to_numeric(
        acquisition_log["Plate index (zero-based)"], errors="raise"
    ).astype(int).tolist()
    acquisition_iterations = pd.to_numeric(
        acquisition_log["Acquisition iteration"], errors="raise"
    ).astype(int).tolist()
    if acquisition_indices != h5_indices or acquisition_iterations != list(range(1, 20)):
        raise ValueError("Processed photoKPFM acquisition log does not preserve HDF5 acquisition order.")
    ranks = frames["candidate_rank_comparison.csv"]
    expected_ranks = set(range(1, 14))
    if set(pd.to_numeric(ranks["PL_rank"], errors="raise").astype(int)) != expected_ranks:
        raise ValueError("Candidate PL ranks are not a complete 1-13 permutation.")
    if set(pd.to_numeric(ranks["photoKPFM_rank"], errors="raise").astype(int)) != expected_ranks:
        raise ValueError("Candidate photoKPFM ranks are not a complete 1-13 permutation.")
    return {
        "pl_dataset_scope": "historical_gain100",
        "datasets": inventory,
        "pl_wells": 96,
        "unique_nominal_compositions": 91,
        "matched_wells": 13,
        "representative_wells": sorted(selected_wells),
    }


def validate_gain60(root: Path) -> dict[str, object]:
    namespace = runpy.run_path(str(root / GAIN60_GENERATOR))
    source_checks = namespace["validate_sources"]()
    if source_checks.get("status") != "passed":
        raise ValueError(f"Gain-60 source validation did not pass: {source_checks}")
    json.dumps(source_checks, allow_nan=False)
    from PIL import Image

    expected_files = {
        Path(f"{stem}{suffix}").name
        for stem, _sources in GAIN60_ASSETS.values()
        for suffix in (".png", ".svg", ".pdf")
    }
    actual_files = {
        path.name for path in (root / "gain60" / "figures").iterdir()
        if path.is_file() and path.suffix.lower() in {".png", ".svg", ".pdf"}
    }
    if actual_files != expected_files:
        raise ValueError(
            f"Gain-60 figure inventory differs: missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}"
        )
    forbidden_annotations = (
        "comparison timepoint", "Normalized PL peak-intensity trajectories",
        "PL and photoKPFM screening", "Candidate reprioritization",
        "PL peak intensity versus photoKPFM response", "Pearson", "bootstrap CI",
    )
    figures = []
    for stem, _sources in GAIN60_ASSETS.values():
        png = root / f"{stem}.png"
        svg = root / f"{stem}.svg"
        pdf = root / f"{stem}.pdf"
        for path in (png, svg, pdf):
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing or empty gain-60 figure: {path.relative_to(root)}")
        with Image.open(png) as image:
            dpi = image.info.get("dpi", ())
            if len(dpi) != 2 or any(abs(float(value) - 600.0) > 0.1 for value in dpi):
                raise ValueError(f"Gain-60 PNG is not 600 dpi: {png.name}, {dpi}")
        svg_root = ET.parse(svg).getroot()
        if any(element.tag.rsplit("}", 1)[-1] == "image" for element in svg_root.iter()):
            raise ValueError(f"Gain-60 SVG contains an embedded raster image: {svg.name}")
        text_nodes = [
            "".join(element.itertext()).strip()
            for element in svg_root.iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
        ]
        if not text_nodes:
            raise ValueError(f"Gain-60 SVG lacks editable scientific labels: {svg.name}")
        visible = " ".join(text_nodes)
        if any(term.casefold() in visible.casefold() for term in forbidden_annotations):
            raise ValueError(f"Gain-60 plot contains a title or statistics annotation: {svg.name}")
        if any(re.fullmatch(r"\(?[aAbB]\)?[.)]?", text) for text in text_nodes):
            raise ValueError(f"Gain-60 plot contains a panel letter: {svg.name}")
        figures.append({"stem": stem, "formats": 3, "png_dpi": 600, "svg_embedded_images": 0})
    return {
        "dataset_scope": "active_gain60",
        "source_checks": source_checks,
        "individual_figures": len(figures),
        "figure_files": len(figures) * 3,
        "figures": figures,
    }


def validate_gain60_figure_s7(root: Path) -> dict[str, object]:
    import numpy as np
    from PIL import Image

    directory = root / GAIN60_S7_DIRECTORY
    expected_files = {
        "FigureS7_gain60.svg", "FigureS7_gain60.pdf", "FigureS7_gain60.png",
        "FigureS7_source_spectra.csv", "FigureS7_plot.py",
        "FigureS7_caption.txt", "FigureS7_validation.txt",
    }
    if not directory.is_dir():
        raise FileNotFoundError(f"Missing Figure S7 directory: {GAIN60_S7_DIRECTORY}")
    actual_files = {path.name for path in directory.iterdir() if path.is_file()}
    if actual_files != expected_files:
        raise ValueError(
            f"Figure S7 inventory differs: missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}"
        )
    if any(not (directory / name).stat().st_size for name in expected_files):
        raise ValueError("Figure S7 contains an empty required file.")
    columns = [
        "well_id", "nominal_FA_pct", "nominal_BDA_pct", "nominal_PEA_pct",
        "acquisition_number", "read_block_number", "gain", "wavelength_nm",
        "raw_intensity_counts", "in_display_band", "display_normalized_PL",
    ]
    frame = pd.read_csv(directory / "FigureS7_source_spectra.csv", keep_default_na=False)
    selected_wells = ("A1", "F2", "G6", "D3", "C11")
    if frame.columns.tolist() != columns or len(frame) != 5 * 416:
        raise ValueError("Figure S7 must contain the declared 11 columns and 2080 spectral rows.")
    if set(frame["well_id"]) != set(selected_wells):
        raise ValueError("Figure S7 well selection differs from A1/F2/G6/D3/C11.")
    for name, value in (("acquisition_number", 22), ("read_block_number", 64), ("gain", 60)):
        if not pd.to_numeric(frame[name], errors="raise").eq(value).all():
            raise ValueError(f"Figure S7 has unexpected {name} values.")
    mask_tokens = frame["in_display_band"].astype(str).str.casefold()
    if not mask_tokens.isin({"true", "false"}).all():
        raise ValueError("Figure S7 display-band flags must be True or False.")
    frame["in_display_band"] = mask_tokens.eq("true")
    for name in ("wavelength_nm", "raw_intensity_counts", "nominal_FA_pct", "nominal_BDA_pct", "nominal_PEA_pct"):
        frame[name] = pd.to_numeric(frame[name], errors="raise")
        if not np.isfinite(frame[name].to_numpy(dtype=float)).all():
            raise ValueError(f"Figure S7 contains nonfinite {name} values.")
    if frame.duplicated(["well_id", "wavelength_nm"]).any():
        raise ValueError("Figure S7 contains duplicate well/wavelength rows.")
    mapping = pd.read_csv(root / "gain60/data/well_map.csv").set_index("well_id")
    with np.load(root / "gain60/data/top_spectra_raw.npz", allow_pickle=False) as source:
        wavelengths = source["wavelengths_nm"].copy()
        source_wells = source["wells"].astype(str).tolist()
        acquisitions = source["acquisition_numbers"].tolist()
        acquisition_index = acquisitions.index(22)
        if int(source["raw_read_numbers"][acquisition_index]) != 64:
            raise ValueError("Figure S7 acquisition 22 does not resolve to source read 64.")
        spectra = source["numeric_values"][acquisition_index].copy()
    if not np.array_equal(wavelengths, np.arange(455.0, 871.0)):
        raise ValueError("Figure S7 source must preserve all 416 observed 455–870 nm channels.")
    maxima = {}
    for well in selected_wells:
        rows = frame.loc[frame["well_id"].eq(well)].sort_values("wavelength_nm")
        if len(rows) != 416 or not np.array_equal(rows["wavelength_nm"].to_numpy(), wavelengths):
            raise ValueError(f"Figure S7 wavelength channels differ from source for {well}.")
        raw = rows["raw_intensity_counts"].to_numpy(dtype=float)
        if not np.array_equal(raw, spectra[source_wells.index(well)]):
            raise ValueError(f"Figure S7 raw intensities differ from source for {well}.")
        for name in ("nominal_FA_pct", "nominal_BDA_pct", "nominal_PEA_pct"):
            if not rows[name].eq(float(mapping.loc[well, name])).all():
                raise ValueError(f"Figure S7 composition differs from source for {well}: {name}")
        display = (wavelengths >= 455.0) & (wavelengths <= 850.0)
        if not np.array_equal(rows["in_display_band"].to_numpy(), display):
            raise ValueError(f"Figure S7 display-band flags differ for {well}.")
        if not rows.loc[~display, "display_normalized_PL"].astype(str).str.strip().eq("").all():
            raise ValueError(f"Figure S7 normalized values must be blank outside 455–850 nm for {well}.")
        maximum = float(np.max(raw[display]))
        if maximum <= 0:
            raise ValueError(f"Figure S7 has no positive display-band maximum for {well}.")
        normalized = pd.to_numeric(rows.loc[display, "display_normalized_PL"], errors="raise").to_numpy(dtype=float)
        if not np.allclose(normalized, raw[display] / maximum, rtol=1e-12, atol=1e-12):
            raise ValueError(f"Figure S7 display normalization differs from source for {well}.")
        maxima[well] = maximum
    with Image.open(directory / "FigureS7_gain60.png") as image:
        dpi = image.info.get("dpi", ())
        if len(dpi) != 2 or any(abs(float(value) - 600.0) > 0.1 for value in dpi):
            raise ValueError("Figure S7 PNG is not 600 dpi.")
    svg_root = ET.parse(directory / "FigureS7_gain60.svg").getroot()
    tags = {element.tag.rsplit("}", 1)[-1] for element in svg_root.iter()}
    if "image" in tags or "text" not in tags:
        raise ValueError("Figure S7 SVG must contain editable labels without embedded raster images.")
    if not (directory / "FigureS7_gain60.pdf").read_bytes().startswith(b"%PDF-"):
        raise ValueError("Figure S7 PDF has an invalid file header.")
    return {
        "dataset_scope": "active_gain60", "files": len(expected_files),
        "selected_wells": list(selected_wells), "source_rows": len(frame),
        "observed_channels_per_well": 416, "display_channels_per_well": 396,
        "comparison_acquisition": 22, "comparison_raw_read": 64, "gain": 60,
        "raw_values_match_source": True, "outside_display_normalization_blank": True,
        "display_band_maxima_counts": maxima, "png_dpi": 600,
    }


def validate_graph(root: Path) -> dict[str, object]:
    literature = root / "data" / "literature"
    nodes = pd.read_csv(literature / "evidence_graph_nodes.csv")
    edges = pd.read_csv(literature / "evidence_graph_edges.csv")
    if len(nodes) != 97 or len(edges) != 120:
        raise ValueError(f"Expected graph size 97 nodes/120 edges, found {len(nodes)}/{len(edges)}.")
    if nodes["node_id"].duplicated().any() or edges["edge_id"].duplicated().any():
        raise ValueError("Duplicate graph node or edge identifiers found.")
    node_ids = set(nodes["node_id"].astype(str))
    dangling = edges.loc[
        ~edges["source_id"].astype(str).isin(node_ids) | ~edges["target_id"].astype(str).isin(node_ids)
    ]
    if not dangling.empty:
        raise ValueError(f"Graph contains {len(dangling)} dangling edges.")
    node_types = nodes.groupby("node_type").size().astype(int).to_dict()
    edge_types = edges.groupby("relationship").size().astype(int).to_dict()
    if node_types != EXPECTED_GRAPH_NODE_TYPES:
        raise ValueError(f"Unexpected graph node-type counts: {node_types}")
    if edge_types != EXPECTED_GRAPH_EDGE_TYPES:
        raise ValueError(f"Unexpected graph edge-type counts: {edge_types}")
    summary_path = literature / "evidence_graph_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    expected_summary = {
        "papers": 7,
        "retained_relationships": 15,
        "nodes": 97,
        "edges": 120,
        "causal_status": ["complementary_only"],
        "steady_state_relationships": 8,
        "trpl_context_relationships": 7,
    }
    if summary != expected_summary:
        raise ValueError(f"Graph summary does not match canonical graph: {summary}")
    stale_counts = {}
    node_count_path = literature / "knowledge_graph_node_counts_corrected.csv"
    edge_count_path = literature / "knowledge_graph_edge_counts_corrected.csv"
    if node_count_path.is_file():
        value = int(pd.to_numeric(pd.read_csv(node_count_path)["count"], errors="raise").sum())
        stale_counts[node_count_path.name] = value
        if value != 97:
            raise ValueError(f"{node_count_path.name} sums to {value}, not 97.")
    if edge_count_path.is_file():
        value = int(pd.to_numeric(pd.read_csv(edge_count_path)["count"], errors="raise").sum())
        stale_counts[edge_count_path.name] = value
        if value != 120:
            raise ValueError(f"{edge_count_path.name} sums to {value}, not 120.")
    return {
        "nodes": len(nodes),
        "edges": len(edges),
        "node_type_counts": node_types,
        "edge_type_counts": edge_types,
        "legacy_count_tables": stale_counts,
    }


def release_figure_files(stems: set[str] | None = None) -> set[str]:
    selected = stems if stems is not None else set(RELEASE_FIGURES)
    return {
        f"{stem}{suffix}"
        for stem in selected
        for suffix in (".png", ".svg", ".pdf")
    }


def validate_ternary_source_data(root: Path) -> dict[str, object]:
    source_paths = {Path(path).name: root / path for path in TERNARY_SOURCE_FILES}
    missing = [path.relative_to(root).as_posix() for path in source_paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Ternary source-data files are missing: {missing}")
    acquisitions = pd.read_csv(source_paths["photoKPFM_acquisitions_all_19.csv"])
    grid = pd.read_csv(source_paths["Figure3_and_FigureS2_GP_grid_source_data.csv"])
    histograms = pd.read_csv(source_paths["Figure3_representative_histograms_source_data.csv"])
    repeats = pd.read_csv(source_paths["FigureS1_repeat_location_variability_source_data.csv"])
    required_acquisition_columns = {
        "acquisition_position_zero_based",
        "plate_index_zero_based",
        "well_id",
        "PEA2PbI4_pct",
        "BDAPbI4_pct",
        "FAPbI3_pct",
        "dark_acquisition",
        "illuminated_acquisition",
        "SP_dark_v",
        "SP_illuminated_v",
        "delta_SP_v",
        "repeat_location_for_FigureS1",
    }
    if not required_acquisition_columns.issubset(acquisitions.columns):
        missing_columns = sorted(required_acquisition_columns.difference(acquisitions.columns))
        raise ValueError(f"Ternary acquisition source is missing columns: {missing_columns}")
    acquisition_positions = pd.to_numeric(
        acquisitions["acquisition_position_zero_based"], errors="raise"
    ).astype(int)
    plate_indices = pd.to_numeric(acquisitions["plate_index_zero_based"], errors="raise").astype(int)
    if len(acquisitions) != 19 or acquisition_positions.tolist() != list(range(19)):
        raise ValueError("Ternary acquisition source must preserve all 19 HDF5 acquisition positions.")
    with h5py.File(root / "data" / "raw" / "photokpfm_measurements.h5", "r") as handle:
        h5_indices = [int(value) for value in handle["idx"][:]]
    if plate_indices.tolist() != h5_indices or acquisitions["well_id"].astype(str).nunique() != 13:
        raise ValueError("Ternary acquisition source does not preserve the HDF5 plate-index sequence.")
    composition_columns = ["PEA2PbI4_pct", "BDAPbI4_pct", "FAPbI3_pct"]
    acquisition_sums = acquisitions[composition_columns].apply(pd.to_numeric, errors="raise").sum(axis=1)
    if acquisition_sums.sub(100.0).abs().gt(1e-6).any():
        raise ValueError("Ternary acquisition compositions do not sum to 100 percent.")
    required_grid_columns = {
        "grid_position_zero_based",
        *composition_columns,
        "GP_mean_delta_SP_v",
        "GP_posterior_standard_deviation_v",
    }
    if not required_grid_columns.issubset(grid.columns) or len(grid) != 91:
        raise ValueError("Figure 3/Figure S2 GP-grid source must contain 91 complete rows.")
    grid_positions = pd.to_numeric(grid["grid_position_zero_based"], errors="raise").astype(int)
    if grid_positions.tolist() != list(range(91)):
        raise ValueError("Figure 3/Figure S2 GP-grid positions must be the ordered range 0-90.")
    grid_sums = grid[composition_columns].apply(pd.to_numeric, errors="raise").sum(axis=1)
    if grid_sums.sub(100.0).abs().gt(1e-6).any():
        raise ValueError("Figure 3/Figure S2 GP-grid compositions do not sum to 100 percent.")
    required_histogram_columns = {
        "panel_number",
        "acquisition_position_zero_based",
        "histogram_bin_center_v",
        "dark_pixel_count",
        "illuminated_pixel_count",
        "dark_gaussian_fit_at_bin",
        "illuminated_gaussian_fit_at_bin",
    }
    if not required_histogram_columns.issubset(histograms.columns) or len(histograms) != 120:
        raise ValueError("Figure 3 histogram source must contain 120 complete rows.")
    histogram_counts = (
        pd.to_numeric(histograms["panel_number"], errors="raise").astype(int).value_counts().to_dict()
    )
    if histogram_counts != {1: 30, 2: 30, 3: 30, 4: 30}:
        raise ValueError(f"Figure 3 histogram panels must each contain 30 bins: {histogram_counts}")
    repeat_mask = acquisitions["repeat_location_for_FigureS1"].astype(str).str.lower().eq("true")
    repeat_positions = set(acquisition_positions[repeat_mask].tolist())
    published_repeat_positions = set(
        pd.to_numeric(repeats["acquisition_position_zero_based"], errors="raise").astype(int).tolist()
    )
    published_repeat_indices = set(
        pd.to_numeric(repeats["plate_index_zero_based"], errors="raise").astype(int).tolist()
    )
    if len(repeats) != 8 or published_repeat_positions != repeat_positions or published_repeat_indices != {0, 78}:
        raise ValueError("Figure S1 source must contain the eight acquisitions at repeated indices 0 and 78.")
    provenance_path = source_paths["ternary_GP_model_and_figure_provenance.json"]
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    ternary_stems = {
        stem
        for stem, (generators, _) in RELEASE_FIGURES.items()
        if TERNARY_GENERATOR in generators
    }
    generated_csvs = set(TERNARY_SOURCE_FILES[:4])
    expected_metadata_outputs = release_figure_files(ternary_stems) | generated_csvs
    metadata_outputs = {str(value).replace("\\", "/") for value in provenance.get("outputs", [])}
    source_h5 = root / "data" / "raw" / "photokpfm_measurements.h5"
    if provenance.get("schema_version") != "1.0":
        raise ValueError("Unexpected ternary model-provenance schema version.")
    if provenance.get("source_h5") != "data/raw/photokpfm_measurements.h5":
        raise ValueError("Ternary model provenance identifies the wrong source HDF5 file.")
    if provenance.get("source_h5_sha256") != sha256(source_h5):
        raise ValueError("Ternary model provenance source-HDF5 hash does not match.")
    if provenance.get("representative_acquisition_positions_zero_based") != [1, 3, 17, 18]:
        raise ValueError("Ternary model provenance has unexpected representative acquisitions.")
    if provenance.get("repeat_location_plate_indices_zero_based") != [0, 78]:
        raise ValueError("Ternary model provenance has unexpected repeated plate indices.")
    if provenance.get("gaussian_process", {}).get("optimizer_success") is not True:
        raise ValueError("Ternary GP provenance does not record a successful optimizer fit.")
    if metadata_outputs != expected_metadata_outputs:
        missing_outputs = sorted(expected_metadata_outputs.difference(metadata_outputs))
        extra_outputs = sorted(metadata_outputs.difference(expected_metadata_outputs))
        raise ValueError(
            f"Ternary provenance output crosswalk differs: missing={missing_outputs}, extra={extra_outputs}"
        )
    manifest_path = source_paths["ternary_photokpfm_output_manifest.csv"]
    manifest = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    if list(manifest.columns) != ["path", "size_bytes", "sha256"]:
        raise ValueError("Ternary output manifest must contain path, size_bytes, and sha256 columns.")
    if manifest["path"].duplicated().any():
        raise ValueError("Ternary output manifest contains duplicate paths.")
    expected_manifest_paths = expected_metadata_outputs | {
        f"{TERNARY_SOURCE_DIRECTORY}/ternary_GP_model_and_figure_provenance.json"
    }
    manifest_paths = set(manifest["path"].str.replace("\\", "/", regex=False))
    if manifest_paths != expected_manifest_paths:
        missing_manifest = sorted(expected_manifest_paths.difference(manifest_paths))
        extra_manifest = sorted(manifest_paths.difference(expected_manifest_paths))
        raise ValueError(
            f"Ternary output manifest crosswalk differs: missing={missing_manifest}, extra={extra_manifest}"
        )
    for row in manifest.to_dict("records"):
        relative = str(row["path"]).replace("\\", "/")
        path = safe_relative(root, relative)
        if not path.is_file():
            raise FileNotFoundError(f"Ternary manifested output is missing: {relative}")
        if int(row["size_bytes"]) != path.stat().st_size or row["sha256"] != sha256(path):
            raise ValueError(f"Ternary manifested output has drifted: {relative}")
    return {
        "acquisitions": len(acquisitions),
        "unique_wells": int(acquisitions["well_id"].astype(str).nunique()),
        "gp_grid_rows": len(grid),
        "histogram_rows": len(histograms),
        "repeat_rows": len(repeats),
        "manifested_outputs": len(manifest),
    }


def safe_relative(root: Path, value: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        raise ValueError(f"Asset index path must be release-relative: {value}")
    resolved = (root / candidate).resolve()
    root_resolved = root.resolve()
    if resolved != root_resolved and root_resolved not in resolved.parents:
        raise ValueError(f"Asset index path escapes the release root: {value}")
    return resolved


def split_paths(value: str) -> list[str]:
    return [item.strip().replace("\\", "/") for item in str(value).split(";") if item.strip()]


def normalized_asset_id(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def expected_asset_index() -> dict[str, dict[str, object]]:
    publication_generator = "scripts/publication/generate_publication_figures.py"
    table_output = (CANONICAL_TABLE_PATH,)
    table_generator = (CANONICAL_TABLE_GENERATOR,)
    figure3_stems = {
        "results/figures/main/Figure3_GP_mean_map",
        "results/figures/main/Figure3_representative_SP_histogram_1",
        "results/figures/main/Figure3_representative_SP_histogram_2",
        "results/figures/main/Figure3_representative_SP_histogram_3",
        "results/figures/main/Figure3_representative_SP_histogram_4",
    }
    assets = {
        'figure1': {'asset_type': 'figure', 'status': 'author-supplied-with-explicit-scope', 'outputs': (AUTHOR_SUPPLIED_ASSETS['figure1'],), 'generators': (), 'sources': ()},
        'figure2': {'asset_type': 'figure', 'status': 'author-supplied-with-explicit-scope', 'outputs': (AUTHOR_SUPPLIED_ASSETS['figure2'],), 'generators': (), 'sources': ()},
        'figure3': {'asset_type': 'figure', 'status': 'reproducible', 'outputs': tuple(sorted(release_figure_files(figure3_stems))), 'generators': (TERNARY_GENERATOR,), 'sources': ('data/raw/photokpfm_measurements.h5',)},
        'figures1': {'asset_type': 'figure', 'status': 'reproducible', 'outputs': tuple(sorted(release_figure_files({'results/figures/supplementary/FigureS1_repeat_location_variability'}))), 'generators': (TERNARY_GENERATOR,), 'sources': ('data/raw/photokpfm_measurements.h5',)},
        'figures2': {'asset_type': 'figure', 'status': 'reproducible', 'outputs': tuple(sorted(release_figure_files({'results/figures/supplementary/FigureS2_GP_uncertainty_map'}))), 'generators': (TERNARY_GENERATOR,), 'sources': ('data/raw/photokpfm_measurements.h5',)},
        'figures3': {'asset_type': 'figure', 'status': 'reproducible', 'outputs': tuple(sorted(release_figure_files({'results/figures/supplementary/FigureS3_plate_composition_map'}))), 'generators': (publication_generator,), 'sources': ('gain60/data/well_map.csv',)},
        'figures6': {'asset_type': 'figure', 'status': 'reproducible', 'outputs': tuple(sorted(release_figure_files({'results/figures/supplementary/FigureS6_evidence_graph'}))), 'generators': ('scripts/literature/plot_evidence_graph.py',), 'sources': ('data/literature/evidence_graph_nodes.csv', 'data/literature/evidence_graph_edges.csv', 'data/literature/table_s4_retained_relationships.csv')},
        'tables1': {'asset_type': 'table', 'status': 'reproducible', 'outputs': table_output, 'generators': table_generator, 'sources': ('data/raw/photokpfm_measurements.h5', 'data/processed/plate_composition_map.csv')},
        'tables2': {'asset_type': 'table', 'status': 'reproducible', 'outputs': table_output, 'generators': table_generator, 'sources': ('data/literature/table_s3_relationship_summary.csv',)},
        'tables3': {'asset_type': 'table', 'status': 'reproducible', 'outputs': table_output, 'generators': table_generator, 'sources': ('data/literature/table_s4_retained_relationships.csv',)},
        'tables4': {'asset_type': 'table', 'status': 'historical-gain100', 'outputs': table_output, 'generators': table_generator, 'sources': ('data/processed/exploratory_association_statistics.csv',)},
    }
    for item_id, (stem, sources) in GAIN60_ASSETS.items():
        assets[item_id] = {
            "asset_type": "figure",
            "status": "reproducible",
            "outputs": tuple(sorted(release_figure_files({stem}))),
            "generators": (GAIN60_GENERATOR,),
            "sources": sources,
        }
    assets["figures7"] = {
        "asset_type": "figure", "status": "reproducible",
        "outputs": tuple(sorted(release_figure_files({GAIN60_S7_STEM}))),
        "generators": (GAIN60_S7_GENERATOR,), "sources": GAIN60_S7_SOURCES,
    }
    return assets


def validate_asset_index(root: Path) -> dict[str, object]:
    index_path = root / "docs" / "MANUSCRIPT_ASSET_INDEX.csv"
    if not index_path.is_file():
        raise FileNotFoundError("Required manuscript asset index is missing: docs/MANUSCRIPT_ASSET_INDEX.csv")
    with index_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = reader.fieldnames or []
    required_columns = [
        "manuscript_id",
        "asset_type",
        "output_paths",
        "generator_script",
        "source_paths",
        "status",
        "scope_note",
    ]
    if columns != required_columns:
        raise ValueError(
            f"MANUSCRIPT_ASSET_INDEX.csv columns must be {required_columns}, found {columns}."
        )
    indexed = {}
    for row_number, row in enumerate(rows, start=2):
        item_id = normalized_asset_id(row["manuscript_id"])
        if not item_id or item_id in indexed:
            raise ValueError(f"Missing or duplicate manuscript_id at asset-index row {row_number}.")
        asset_type = row["asset_type"].strip().lower()
        if asset_type not in {"figure", "table"}:
            raise ValueError(f"Invalid asset_type at asset-index row {row_number}: {asset_type}")
        status = row["status"].strip().lower()
        if status not in ALLOWED_ASSET_STATUSES:
            raise ValueError(f"Invalid status at asset-index row {row_number}: {status}")
        output_paths = split_paths(row["output_paths"])
        if not output_paths:
            raise ValueError(f"No output_paths at asset-index row {row_number}.")
        for value in output_paths:
            path = safe_relative(root, value)
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing or empty manuscript asset: {value}")
        source_paths = split_paths(row["source_paths"])
        generators = split_paths(row["generator_script"])
        if status in {"reproducible", "historical-gain100"}:
            if not generators:
                raise ValueError(f"Missing generator_script at asset-index row {row_number}.")
            for generator in generators:
                if not safe_relative(root, generator).is_file():
                    raise ValueError(f"Missing indexed generator script: {generator}")
            if not source_paths:
                raise ValueError(f"No source_paths for reproducible asset at row {row_number}.")
            for value in source_paths:
                if not safe_relative(root, value).is_file():
                    raise ValueError(f"Missing indexed source file: {value}")
            if status == "historical-gain100" and not re.search(
                r"historical.*gain.?100|gain.?100.*historical", row["scope_note"], re.IGNORECASE
            ):
                raise ValueError(f"Historical gain-100 asset lacks dataset scope at row {row_number}.")
        else:
            if generators or source_paths:
                raise ValueError(f"Author-supplied asset claims a generator or source at row {row_number}.")
            if not row["scope_note"].strip():
                raise ValueError(f"Author-supplied asset lacks explicit scope at row {row_number}.")
        indexed[item_id] = {
            "asset_type": asset_type,
            "status": status,
            "outputs": tuple(sorted(output_paths)),
            "generators": tuple(sorted(generators)),
            "sources": tuple(sorted(source_paths)),
        }
    missing_ids = sorted(REQUIRED_MANUSCRIPT_IDS.difference(indexed))
    extra_ids = sorted(set(indexed).difference(REQUIRED_MANUSCRIPT_IDS))
    if missing_ids or extra_ids:
        raise ValueError(
            f"Manuscript asset-index identifiers differ: missing={missing_ids}, extra={extra_ids}"
        )
    expected = expected_asset_index()
    mismatches = {}
    for item_id in sorted(expected):
        expected_record = expected[item_id]
        actual_record = indexed[item_id]
        differences = {}
        for field in ("asset_type", "status", "outputs", "generators", "sources"):
            expected_value = expected_record[field]
            if isinstance(expected_value, tuple):
                expected_value = tuple(sorted(expected_value))
            if actual_record[field] != expected_value:
                differences[field] = {
                    "expected": expected_value,
                    "actual": actual_record[field],
                }
        if differences:
            mismatches[item_id] = differences
    if mismatches:
        raise ValueError(f"Manuscript asset-index mappings differ from the release crosswalk: {mismatches}")
    status_counts = defaultdict(int)
    for row in indexed.values():
        status_counts[row["status"]] += 1
    return {
        "indexed_assets": len(rows),
        "required_manuscript_assets": len(REQUIRED_MANUSCRIPT_IDS),
        "status_counts": dict(sorted(status_counts.items())),
    }


def validate_release_assets(root: Path) -> dict[str, object]:
    missing = []
    issues = []
    for stem, (generators, sources) in RELEASE_FIGURES.items():
        for suffix in (".png", ".svg", ".pdf"):
            path = root / f"{stem}{suffix}"
            if not path.is_file() or not path.stat().st_size:
                missing.append(path.relative_to(root).as_posix())
        for generator in generators:
            generator_path = root / generator
            if not generator_path.is_file():
                missing.append(generator)
        for source in sources:
            if not (root / source).is_file():
                missing.append(source)
    for asset in AUTHOR_SUPPLIED_ASSETS.values():
        path = root / asset
        if not path.is_file() or not path.stat().st_size:
            missing.append(asset)
    for source in TERNARY_SOURCE_FILES:
        path = root / source
        if not path.is_file() or not path.stat().st_size:
            missing.append(source)
    table = root / CANONICAL_TABLE_PATH
    table_provenance_path = root / CANONICAL_TABLE_PROVENANCE_PATH
    table_generator = root / CANONICAL_TABLE_GENERATOR
    for path in (table, table_provenance_path, table_generator):
        if not path.is_file() or not path.stat().st_size:
            missing.append(path.relative_to(root).as_posix())
    if missing:
        issues.append(f"release crosswalk missing files: {sorted(set(missing))}")
    stale_artifacts = [
        f"{stem}{suffix}"
        for stem in STALE_FIGURE_STEMS
        for suffix in (".png", ".svg", ".pdf")
        if (root / f"{stem}{suffix}").exists()
    ]
    stale_artifacts.extend(path for path in STALE_RELEASE_PATHS if (root / path).exists())
    if stale_artifacts:
        issues.append(f"stale release artifacts remain: {sorted(stale_artifacts)}")
    expected_table_rows = {"S1": 19, "S2": 6, "S3": 15, "S4": 36}
    if table.is_file() and table_provenance_path.is_file():
        try:
            table_provenance = json.loads(table_provenance_path.read_text(encoding="utf-8"))
            actual_table_rows = {
                label: int(record["rows"])
                for label, record in table_provenance.get("tables", {}).items()
            }
            if table_provenance.get("status") != "complete" or actual_table_rows != expected_table_rows:
                raise ValueError(f"incomplete row record {actual_table_rows}")
            if (
                table_provenance.get("output") != table.name
                or table_provenance.get("output_sha256") != sha256(table)
            ):
                raise ValueError("output hash mismatch")
            for relative, digest in table_provenance.get("source_sha256", {}).items():
                source = safe_relative(root, relative)
                if not source.is_file() or sha256(source) != digest:
                    raise ValueError(f"source hash mismatch: {relative}")
        except Exception as error:
            issues.append(f"supplementary-table provenance invalid: {type(error).__name__}: {error}")
    provenance_path = root / "docs" / "FIGURE_PROVENANCE.md"
    if provenance_path.is_file():
        provenance = provenance_path.read_text(encoding="utf-8")
        provenance_gaps = []
        gain60_stems = {stem for stem, _sources in GAIN60_ASSETS.values()}
        gain60_stems.add(GAIN60_S7_STEM)
        for stem, (generators, sources) in RELEASE_FIGURES.items():
            required_tokens = [
                Path(stem).name,
                *(Path(generator).name for generator in generators),
                *(Path(source).name for source in sources if stem not in gain60_stems),
            ]
            absent = [token for token in required_tokens if token not in provenance]
            if absent:
                provenance_gaps.append({"stem": Path(stem).name, "missing_tokens": absent})
        for item_id, asset in AUTHOR_SUPPLIED_ASSETS.items():
            if Path(asset).name not in provenance:
                provenance_gaps.append({"asset": item_id, "missing_tokens": [Path(asset).name]})
        if provenance_gaps:
            issues.append(f"FIGURE_PROVENANCE.md mappings incomplete: {provenance_gaps}")
    try:
        indexed = validate_asset_index(root)
    except Exception as error:
        indexed = {}
        issues.append(f"manuscript asset index invalid: {type(error).__name__}: {error}")
    if issues:
        raise ValueError("; ".join(issues))
    return {
        "figure_stems": len(RELEASE_FIGURES),
        "figure_files": len(RELEASE_FIGURES) * 3,
        "author_supplied_figures": len(AUTHOR_SUPPLIED_ASSETS),
        "ternary_source_files": len(TERNARY_SOURCE_FILES),
        "table_files": 2,
        "supplementary_table_rows": expected_table_rows,
        "manuscript_index": indexed,
    }


def validate_python(root: Path) -> dict[str, object]:
    scripts = sorted(path for path in root.rglob("*.py") if included_path(path, root))
    forbidden_comments = []
    remaining_comments = []
    absolute_paths = []
    attribution_terms = []
    for path in scripts:
        source = path.read_text(encoding="utf-8")
        compile(source, str(path), "exec")
        relative = path.relative_to(root).as_posix()
        for line_number, comment in actual_comments(source):
            remaining_comments.append({"path": relative, "line": line_number, "comment": comment})
            if not allowed_comment(line_number, comment):
                forbidden_comments.append({"path": relative, "line": line_number, "comment": comment})
        for line_number, line in enumerate(source.splitlines(), start=1):
            if re.search(r"(?<![A-Za-z0-9_])[A-Za-z]:\\(?=[A-Za-z_])", line):
                absolute_paths.append({"path": relative, "line": line_number, "text": line.strip()})
            attribution_pattern = r"\b(?:" + "|".join(("bo" + "ris", "chat" + "gpt", "clau" + "de")) + r")\b"
            if re.search(attribution_pattern, line, flags=re.IGNORECASE):
                attribution_terms.append({"path": relative, "line": line_number, "text": line.strip()})
    if forbidden_comments:
        raise ValueError(f"Nonessential Python comments remain: {forbidden_comments[:5]}")
    if absolute_paths:
        raise ValueError(f"Absolute Windows paths remain in scripts: {absolute_paths[:5]}")
    if attribution_terms:
        raise ValueError(f"Personal or tool-specific terms remain in scripts: {attribution_terms[:5]}")
    return {
        "python_files": len(scripts),
        "syntax_valid": True,
        "remaining_comments": remaining_comments,
        "nonessential_comments": 0,
        "absolute_windows_paths": 0,
    }


def integer_value(node: ast.AST | None) -> int | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return int(node.value)
    return None


def call_dimensions(call: ast.Call) -> tuple[ast.AST | None, ast.AST | None]:
    rows = call.args[0] if call.args else None
    columns = call.args[1] if len(call.args) > 1 else None
    for keyword in call.keywords:
        if keyword.arg == "nrows":
            rows = keyword.value
        elif keyword.arg == "ncols":
            columns = keyword.value
    return rows, columns


def has_main_entrypoint(tree: ast.AST) -> bool:
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        try:
            expression = ast.unparse(node.test)
        except AttributeError:
            expression = ""
        if "__name__" in expression and "__main__" in expression:
            return True
    return False


def plotting_call_violations(
    tree: ast.AST,
    relative: str,
    cell_index: int | None = None,
) -> list[dict[str, object]]:
    violations = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = ""
        if isinstance(node.func, ast.Attribute):
            name = node.func.attr
        elif isinstance(node.func, ast.Name):
            name = node.func.id
        location: dict[str, object] = {"path": relative, "line": node.lineno, "call": name}
        if cell_index is not None:
            location["cell"] = cell_index
        if name == "subplot_mosaic":
            violations.append(location)
            continue
        if name in {"subplots", "add_gridspec", "GridSpec"}:
            row_node, column_node = call_dimensions(node)
            rows = integer_value(row_node)
            columns = integer_value(column_node)
            unsafe_rows = row_node is not None and rows != 1
            unsafe_columns = column_node is not None and columns != 1
            if unsafe_rows or unsafe_columns:
                location.update({"rows": rows, "columns": columns})
                violations.append(location)
            continue
        if name not in {"subplot", "add_subplot"}:
            continue
        rows = None
        columns = None
        if len(node.args) >= 3:
            rows = integer_value(node.args[0])
            columns = integer_value(node.args[1])
            if rows != 1 or columns != 1:
                location.update({"rows": rows, "columns": columns})
                violations.append(location)
        elif node.args:
            encoded = integer_value(node.args[0])
            if encoded is not None and 100 <= encoded <= 999:
                rows = encoded // 100
                columns = (encoded // 10) % 10
                if rows != 1 or columns != 1:
                    location.update({"rows": rows, "columns": columns})
                    violations.append(location)
    return violations


def notebook_python_source(source: str) -> str:
    lines = source.splitlines(keepends=True)
    first_code = next((line.lstrip() for line in lines if line.strip()), "")
    if first_code.startswith("%%"):
        return "".join("\n" if line.endswith(("\n", "\r")) else "" for line in lines)
    return "".join(
        "\n" if line.lstrip().startswith(("!", "%")) and line.endswith(("\n", "\r")) else
        "" if line.lstrip().startswith(("!", "%")) else line
        for line in lines
    )


def validate_single_panel_sources(root: Path) -> dict[str, object]:
    active = []
    violations = []
    plotting_sources = [*(root / "scripts").rglob("*.py"), *(root / "gain60").rglob("*.py")]
    for path in sorted(plotting_sources):
        source = path.read_text(encoding="utf-8")
        if "matplotlib" not in source:
            continue
        tree = ast.parse(source, filename=str(path))
        if not has_main_entrypoint(tree):
            continue
        relative = path.relative_to(root).as_posix()
        active.append(relative)
        violations.extend(plotting_call_violations(tree, relative))
    notebook_cells = 0
    for path in sorted((root / "notebooks").rglob("*.ipynb")):
        relative = path.relative_to(root).as_posix()
        notebook = json.loads(path.read_text(encoding="utf-8"))
        for cell_index, cell in enumerate(notebook.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue
            source = "".join(cell.get("source", []))
            if "matplotlib" not in source and "plt." not in source:
                continue
            notebook_cells += 1
            prepared = notebook_python_source(source)
            try:
                tree = ast.parse(prepared, filename=f"{path}#cell-{cell_index}")
            except SyntaxError as error:
                violations.append(
                    {
                        "path": relative,
                        "cell": cell_index,
                        "line": error.lineno,
                        "call": "unparseable_notebook_code",
                        "message": error.msg,
                    }
                )
                continue
            violations.extend(plotting_call_violations(tree, relative, cell_index))
    if violations:
        raise ValueError(f"Plotting sources create multipanel figures: {violations}")
    return {
        "active_plotting_entrypoints": active,
        "notebook_plotting_cells": notebook_cells,
        "multipanel_calls": 0,
    }


def validate_notebook(root: Path) -> dict[str, object]:
    notebook_path = root / "notebooks" / "pl_phase_fitting_colab.ipynb"
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    output_count = 0
    comment_count = 0
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        output_count += len(cell.get("outputs", []))
        source = "".join(cell.get("source", []))
        try:
            comment_count += sum(
                1
                for line_number, comment in actual_comments(source)
                if not allowed_comment(line_number, comment)
            )
        except (IndentationError, tokenize.TokenError):
            pass
    if output_count:
        raise ValueError(f"Notebook still contains {output_count} saved outputs.")
    if comment_count:
        raise ValueError(f"Notebook still contains {comment_count} nonessential code comments.")
    return {"saved_outputs": output_count, "nonessential_code_comments": comment_count}


def validate_names(root: Path) -> dict[str, object]:
    attribution_terms = "|".join(("bo" + "ris", "chat" + "gpt", "clau" + "de"))
    forbidden = re.compile(
        rf"(?:{attribution_terms}|\bold\b|\bcopy\b|\(1\)|superseded|changelog|\bdraft\b|\brevised\b)",
        re.IGNORECASE,
    )
    bad = [
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if included_path(path, root) and forbidden.search(path.name)
    ]
    if bad:
        raise ValueError(f"Unclear or personal filenames remain: {bad[:10]}")
    return {"forbidden_names": 0}


def inventory_rows(root: Path) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    files = sorted(
        path
        for path in root.rglob("*")
        if path.is_file() and path.name not in TRACKING_FILES and included_path(path, root)
    )
    rows = []
    hashes: dict[str, list[str]] = defaultdict(list)
    for path in files:
        relative = path.relative_to(root).as_posix()
        digest = sha256(path)
        hashes[digest].append(relative)
        rows.append(
            {
                "relative_path": relative,
                "bytes": path.stat().st_size,
                "sha256": digest,
                "role": role_for(relative),
            }
        )
    duplicate_rows = []
    for digest, paths in sorted(hashes.items()):
        if len(paths) < 2:
            continue
        group_id = digest[:12]
        for path in paths:
            duplicate_rows.append({"duplicate_group": group_id, "sha256": digest, "relative_path": path})
    return rows, duplicate_rows


def read_manifest(root: Path) -> list[dict[str, str]]:
    path = root / MANIFEST_NAME
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    expected_columns = {"relative_path", "bytes", "sha256", "role"}
    columns = set(rows[0]) if rows else set()
    if columns != expected_columns:
        raise ValueError(f"Manifest columns must be {sorted(expected_columns)}, found {sorted(columns)}.")
    return rows


def compare_manifest(root: Path, rows: list[dict[str, object]]) -> dict[str, object]:
    expected_rows = read_manifest(root)
    expected = {row["relative_path"]: row for row in expected_rows}
    actual = {str(row["relative_path"]): row for row in rows}
    missing = sorted(set(expected).difference(actual))
    extra = sorted(set(actual).difference(expected))
    mismatches = []
    for relative in sorted(set(expected).intersection(actual)):
        old = expected[relative]
        new = actual[relative]
        if (
            str(old["bytes"]) != str(new["bytes"])
            or str(old["sha256"]).lower() != str(new["sha256"]).lower()
            or str(old["role"]) != str(new["role"])
        ):
            mismatches.append(
                {
                    "relative_path": relative,
                    "manifest_bytes": old["bytes"],
                    "actual_bytes": new["bytes"],
                    "manifest_sha256": old["sha256"],
                    "actual_sha256": new["sha256"],
                    "manifest_role": old["role"],
                    "actual_role": new["role"],
                }
            )
    if missing or extra or mismatches:
        raise ValueError(
            f"Manifest mismatch: missing={missing}, extra={extra}, changed={mismatches}"
        )
    return {"manifest_rows": len(expected_rows), "missing": 0, "extra": 0, "changed": 0}


def write_inventory(root: Path, rows: list[dict[str, object]], duplicate_rows: list[dict[str, str]]) -> None:
    with tempfile.TemporaryDirectory(prefix=".release_inventory_", dir=root) as directory:
        temporary = Path(directory)
        manifest = temporary / MANIFEST_NAME
        duplicates = temporary / DUPLICATE_AUDIT_NAME
        with manifest.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["relative_path", "bytes", "sha256", "role"])
            writer.writeheader()
            writer.writerows(rows)
        with duplicates.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["duplicate_group", "sha256", "relative_path"])
            writer.writeheader()
            writer.writerows(duplicate_rows)
        os.replace(manifest, root / MANIFEST_NAME)
        os.replace(duplicates, root / DUPLICATE_AUDIT_NAME)


def write_report(root: Path, report: dict[str, object]) -> None:
    with tempfile.TemporaryDirectory(prefix=".release_report_", dir=root) as directory:
        temporary = Path(directory) / REPORT_NAME
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, root / REPORT_NAME)


def validate_inventory(root: Path, compare_existing: bool) -> dict[str, object]:
    rows, duplicate_rows = inventory_rows(root)
    if duplicate_rows:
        raise ValueError(f"Exact duplicate files remain: {duplicate_rows[:6]}")
    result = {
        "files": len(rows),
        "bytes": int(sum(int(row["bytes"]) for row in rows)),
        "exact_duplicate_groups": 0,
    }
    if compare_existing:
        result["manifest"] = compare_manifest(root, rows)
    else:
        result["manifest"] = {"mode": "pending_refresh"}
    return result


def run_check(report: dict[str, object], failures: list[dict[str, str]], name: str, function) -> None:
    try:
        report[name] = function()
    except Exception as error:
        failures.append({"check": name, "error": f"{type(error).__name__}: {error}"})


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    report: dict[str, object] = {
        "status": "running", "mode": "read_only",
        "active_PL_dataset": "gain60", "historical_PL_dataset": "gain100",
    }
    failures: list[dict[str, str]] = []
    run_check(report, failures, "required_files", lambda: validate_required_files(root, args.refresh_manifest))
    run_check(report, failures, "names", lambda: validate_names(root))
    run_check(report, failures, "relationships", lambda: validate_relationships(root))
    run_check(report, failures, "raw_data", lambda: validate_raw_data(root))
    run_check(report, failures, "processed_data", lambda: validate_processed_data(root))
    run_check(report, failures, "gain60", lambda: validate_gain60(root))
    run_check(report, failures, "figure_s7", lambda: validate_gain60_figure_s7(root))
    run_check(report, failures, "graph", lambda: validate_graph(root))
    run_check(report, failures, "ternary_source_data", lambda: validate_ternary_source_data(root))
    run_check(report, failures, "assets", lambda: validate_release_assets(root))
    run_check(report, failures, "python", lambda: validate_python(root))
    run_check(report, failures, "single_panel_sources", lambda: validate_single_panel_sources(root))
    run_check(report, failures, "notebook", lambda: validate_notebook(root))
    run_check(report, failures, "inventory", lambda: validate_inventory(root, not args.refresh_manifest))
    if args.refresh_manifest and not failures:
        try:
            rows, duplicate_rows = inventory_rows(root)
            write_inventory(root, rows, duplicate_rows)
            report["inventory"] = validate_inventory(root, True)
            report["mode"] = "manifest_refreshed"
        except Exception as error:
            failures.append({"check": "manifest_refresh", "error": f"{type(error).__name__}: {error}"})
    report["status"] = "failed" if failures else "passed"
    report["failures"] = failures
    if args.write_report:
        write_report(root, report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
