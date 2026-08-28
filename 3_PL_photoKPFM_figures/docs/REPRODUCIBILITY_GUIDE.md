# Reproducibility guide

Run every command from the extracted release root: the directory containing `README.md`, `requirements.txt`, `data/`, and `scripts/`. Python 3.11 is recommended. The code resolves files relative to the release and does not require a user-specific drive, home directory, registry change, or administrator access.

## 1. Install on Windows PowerShell

```powershell
py -3.11 -m venv .venv
& .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If the `py` launcher is unavailable but `python` already selects Python 3.11, use `python -m venv .venv` for the first command. If PowerShell reports that script execution is disabled, apply a process-only policy in the current window and activate again:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
& .\.venv\Scripts\Activate.ps1
```

## 2. Install on macOS, Linux, or another POSIX shell

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

If the platform names Python 3.11 explicitly, substitute `python3.11` for `python3` in the first command.

## 3. Validate the supplied release

```text
python scripts/validate_release.py
```

On Windows, the same command can use backslashes:

```powershell
python scripts\validate_release.py
```

Validation is read-only unless `--refresh-manifest` or `--write-report` is explicitly supplied. It checks required files, Python syntax, cleared notebook outputs, raw and processed dimensions, retained relationship counts and causal status, graph structure, figure triplets, supplementary-table provenance, manuscript asset-index paths, personal paths, unclear filenames, manifest hashes, and duplicate content. A validation failure is not a warning: read the named path or invariant before reproducing the analyses.

## 4. Reproduce the complete release in isolation

```text
python scripts/reproduce_all.py
```

This command is identical on Windows PowerShell and POSIX shells; Python accepts the forward-slash path on both.

The command creates an isolated working copy under `reproduced_results/`, then:

1. regenerates Figures 4a, 4b, 5a, 5b, S3, and S4;
2. regenerates the Figure 3 GP map and four histograms plus Figures S1 and S2 from the ternary HDF5 file;
3. regenerates the combined Tables S1–S4 document and Table S1 acquisition log;
4. reruns the PL metric, spectral-processing, and exploratory association checks;
5. rebuilds the literature evidence-graph tables and Figure S6; and
6. copies the explicitly scoped Figure 1 and Figure 2 reference images.

The supplied `data/`, `results/`, and manuscript assets are not modified. On success, inspect:

- `reproduced_results/reproduction_summary.json` for `"status": "passed"` and stage return codes;
- `reproduced_results/results/figures/` for figure panels;
- `reproduced_results/results/source_data/` for panel source data and provenance;
- `reproduced_results/results/tables/` for Tables S1–S4; and
- `reproduced_results/results/validation/` for analysis checks.

A second run refuses to replace the existing directory. Replace it explicitly with:

```text
python scripts/reproduce_all.py --overwrite
```

`--skip-validation-analyses` omits the two longer independent PL/spectral validation stages; it should not be used for a full reviewer verification.

## 5. Use individual generators directly

Use separate output directories when inspecting one workflow. These examples keep generated figures away from the supplied authoritative exports. They run as written in Windows PowerShell and POSIX shells; PowerShell users may substitute backslashes in paths if preferred.

### PL publication panels

```text
python scripts/publication/generate_publication_figures.py --out-dir reproduced_results/direct_pl_panels
```

Add `--overwrite` only when intentionally replacing that complete output directory.

### Ternary photoKPFM panels

```text
python scripts/publication/generate_ternary_photokpfm_panels.py --input-h5 data/raw/photokpfm_measurements.h5 --output-root reproduced_results/direct_ternary
```

This writes one Figure 3 mean map, four separate Figure 3 histogram panels, Figure S1, Figure S2, their PNG/SVG/PDF variants, machine-readable source data, and GP provenance. Choose a new output root for an immutable rerun; this generator refreshes same-named files when a target already exists.

### Evidence graph and Figure S6

```text
python scripts/literature/build_evidence_graph.py --relationships data/literature/table_s4_retained_relationships.csv --nodes-out reproduced_results/direct_graph/evidence_graph_nodes.csv --edges-out reproduced_results/direct_graph/evidence_graph_edges.csv --summary-out reproduced_results/direct_graph/evidence_graph_summary.json
python scripts/literature/plot_evidence_graph.py --nodes reproduced_results/direct_graph/evidence_graph_nodes.csv --edges reproduced_results/direct_graph/evidence_graph_edges.csv --relationships data/literature/table_s4_retained_relationships.csv --output-stem reproduced_results/direct_graph/FigureS6_evidence_graph
```

### Tables S1–S4

```text
python scripts/publication/generate_supplementary_tables.py --output reproduced_results/direct_tables/supplementary_tables_current_manuscript_S1-S4.docx --overwrite
```

The table generator writes the requested DOCX and adjacent provenance JSON and also refreshes `data/processed/photokpfm_acquisition_log.csv`, the machine-readable Table S1. Use the isolated `reproduce_all.py` route when the supplied tree must remain byte-for-byte unchanged.

### Exploratory matched association analysis

```text
python scripts/analysis/analyze_pl_photokpfm_associations.py --matched-csv data/processed/matched_pl_photokpfm_metrics.csv --phase-csv data/processed/initial_phase_kpfm_table.csv --literature-relationships data/literature/table_s4_retained_relationships.csv --out-dir reproduced_results/direct_association
```

## 6. Expected invariants

| Check | Expected value |
|---|---:|
| PL plate wells | 96 |
| Unique nominal PL compositions | 91 |
| Complete PL top-read timepoints | 22 |
| Parsed PL well-timepoint records | 2,112 |
| Ternary photoKPFM acquisition pairs | 19 |
| Unique photoKPFM wells | 13 |
| Dark and illuminated image arrays | 19×128×128 each |
| Retained literature relationships | 15 |
| DOI-linked primary papers | 7 |
| Evidence-graph nodes / edges | 97 / 120 |
| Relationship causal status | `complementary_only` |
| Exploratory association tests passing BH 0.05 | 0 of 36 |
| Tables S1 / S2 / S3 / S4 data rows | 19 / 6 / 15 / 36 |

The PL/photoKPFM associations are exploratory. Successful execution does not convert them into causal evidence.

## 7. Interactive PL emissive-component fitting

`notebooks/pl_phase_fitting_colab.ipynb` is a separate Colab route because it uses upload widgets. Its saved outputs and execution counters are cleared. The checked downstream table used by the local association workflow is supplied as `data/processed/initial_phase_kpfm_table.csv`; see `notebooks/README.md` before rerunning the fit.

The fitted emissive components are spectral-model outputs, not XRD-confirmed structural phase fractions.

The notebook follows the same single-panel policy as the local plotting scripts. Each dark/illuminated KPFM distribution and each per-well PL fit example is saved separately. The phase-distribution, KPFM-response, and dominant-PL-peak views are standalone files. The optional phase-plus-KPFM overlay is one plot with a shared composition axis and secondary y-axis, not a multi-panel grid. The notebook does not generate a manuscript composite or apply panel lettering.

## 8. Figure and manuscript assembly boundary

Each plotting script produces one scientific panel per file. Figure 3's map and four histograms remain separate outputs. Panel letters, composite layout, captions, and placement are applied outside the scripts during manuscript assembly. PNG, SVG, and PDF files with the same stem are format variants of one panel.
