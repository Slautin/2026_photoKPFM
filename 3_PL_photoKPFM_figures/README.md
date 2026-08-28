# PL and photoKPFM analysis of PEA/BDA/FA perovskite libraries

This is the reviewer-facing reproducibility release for the PL/photoKPFM manuscript. It bundles the original plate-reader PL export, the ternary-library photoKPFM HDF5 file, manuscript-level processed tables, source-audited literature relationships, validation code, manuscript-numbered supplementary tables, and single-panel figure exports.

Figure 1 and Figure 2 are included as rendered author-supplied reference images. The editable Figure 1 workflow source and the binary-library raw data and original GP workflow behind Figure 2 were not available for this release, so neither figure is claimed as computationally reproducible. The bundled photoKPFM HDF5 file is the ternary-library dataset used for Figure 3, Figure S1, Figure S2, and Table S1.

## Quick start on Windows PowerShell

Open PowerShell in the extracted directory that contains this README, then run:

```powershell
py -3.11 -m venv .venv
& .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts\validate_release.py
python scripts\reproduce_all.py
```

If the Windows Python launcher is unavailable but `python` already selects Python 3.11, use `python -m venv .venv` for the first command. If local PowerShell policy blocks activation, run `Set-ExecutionPolicy -Scope Process Bypass` in that same window and retry activation. Administrator access, registry edits, and machine-specific paths are not required.

## Quick start on macOS, Linux, or another POSIX shell

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/validate_release.py
python scripts/reproduce_all.py
```

`scripts/reproduce_all.py` works in an isolated copy and writes only to `reproduced_results/`. It refuses to replace an existing output directory unless invoked with `--overwrite`. A successful run writes `reproduced_results/reproduction_summary.json` with `"status": "passed"`.

## Bundled raw measurements

| File | Expected contents |
|---|---|
| `data/raw/pl_plate_reader_export.csv` | Original, nonrectangular Cytation 5 export: 90 labeled instrument read blocks and 43,818 physical CSV rows. The validated PL time course uses 22 complete top-read blocks (read numbers 1, 4, ..., 64), each containing all 96 wells, for 2,112 well-timepoint records. The plate contains 91 unique nominal compositions because five compositions were repeated. |
| `data/raw/photokpfm_measurements.h5` | Ternary-library data with 19 paired dark/illuminated acquisitions across 13 unique plate wells. `dark_data` and `light_data` are each 19×128×128; `X_train` is 19×2; `y_train` and `idx` each contain 19 values; `X` is 91×2; and `coord_array` is 96×2. The manuscript reports 5×5 µm² scans. Exact timestamps and an external absolute-potential calibration are not inferred. |

The validator checks these dimensions and counts before any analysis is reproduced. See `docs/DATA_PROVENANCE.md` for the complete HDF5 inventory and `docs/DATA_DICTIONARY.md` for derived tables.

## Repository layout

- `data/raw`: original PL export and ternary photoKPFM HDF5 measurements.
- `data/processed`: manuscript-level derived tables and exact plotting inputs.
- `data/validation`: spectral-processing, baseline, time-course, and sensitivity checks.
- `data/literature`: source-audited literature relationships and evidence-graph data.
- `scripts/publication`: single-panel publication-figure and supplementary-table generators.
- `scripts/experimental`: raw-data processing and independent analytical checks.
- `scripts/analysis`: matched PL/photoKPFM association analysis.
- `scripts/literature`: evidence-graph construction and plotting.
- `notebooks`: output-free Colab workflow for PL emissive-component fitting.
- `results/figures`: manuscript-numbered figure exports plus explicitly scoped author-supplied references.
- `results/source_data`: machine-readable source data and model provenance emitted for publication panels.
- `results/tables`: the current combined Tables S1–S4 document and provenance record.
- `docs/MANUSCRIPT_ASSET_INDEX.csv`: authoritative manuscript-ID-to-file crosswalk.

## Reproduction routes

The one-command route stages `data/` and `scripts/` in a clean workspace; regenerates the PL, ternary photoKPFM, literature, and table outputs; reruns validation analyses; copies the two explicitly scoped author reference images; and then moves the completed results into `reproduced_results/results/`.

Individual generators may also be run directly with reviewer-selected output locations. Exact commands and overwrite behavior are documented in `docs/REPRODUCIBILITY_GUIDE.md`.

All plotting scripts produce one scientific panel per file. Figure 3 is represented by one GP mean map and four separate histogram files. Composite assembly, panel lettering, resizing, and placement in the manuscript occur outside the plotting scripts and are not separate analyses. PNG, SVG, and PDF files are intentional format variants of a panel.

The interactive Colab notebook follows the same boundary: per-well PL fits and dark/illuminated KPFM distributions are separate files, and the phase-distribution, KPFM-response, and dominant-peak views are standalone. It does not produce the manuscript composite.

## Scientific scope

PL and photoKPFM records are matched by nominal composition, not spatially registered fields of view. Their associations are exploratory and noncausal. TRPL/carrier-lifetime relationships are broader literature context because TRPL was not measured in this study. PL-derived emissive components are spectral-model outputs rather than structurally confirmed phase fractions.

The full LiteratureAgent extraction engine is a separate project. This release contains only the source-audited relationships and graph code required for the manuscript; it excludes model caches, credentials, publisher article files, and raw language-model responses.

See `docs/FIGURE_PROVENANCE.md` for figure-level sources and generators and `docs/MANUSCRIPT_ASSET_INDEX.csv` for the complete manuscript crosswalk.
