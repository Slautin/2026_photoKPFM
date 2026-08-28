# Data provenance

## Raw PL export

`data/raw/pl_plate_reader_export.csv` is the original Cytation 5 plate-reader export for the 96-well PEA/BDA/FA compositional library. It is an instrument report rather than a rectangular analysis table: the file has 43,818 physical rows, as many as 98 fields per row, and 90 labeled read blocks. The validated longitudinal analysis uses the 22 complete top-read emission-spectrum blocks numbered 1, 4, 7, ..., 64. Each selected block contains all 96 wells, yielding 2,112 well-timepoint observations. The plate represents 91 unique nominal composition triplets because five compositions were repeated.

The raw export is retained unchanged. `data/processed/pl_retention_metrics.csv`, the matched PL tables, and the trace tables are derived release inputs. The representative-trace and all-matched-trace tables preserve plotted numerical coordinates; their endpoints are checked against the independently calculated PL-retention metrics.

## Raw ternary photoKPFM HDF5 file

`data/raw/photokpfm_measurements.h5` contains the ternary-library data used for Figure 3, Figure S1, Figure S2, and Table S1. Its required datasets are:

| Dataset | Expected shape | Interpretation |
|---|---:|---|
| `X` | 91×2 | Unique ternary composition coordinates used for prediction. |
| `X_train` | 19×2 | Acquisition-level GP training coordinates. |
| `coord_array` | 96×2 | Coordinates for the complete 96-well plate layout. |
| `dark_data` | 19×128×128 | Dark surface-potential images. |
| `dark_fn` | 19 | Dark acquisition identifiers. |
| `idx` | 19 | Zero-based plate indices; 13 unique wells are represented. |
| `light_data` | 19×128×128 | Illuminated surface-potential images paired with `dark_data`. |
| `light_fn` | 19 | Illuminated acquisition identifiers. |
| `y_train` | 19 | Training-response values retained in the source HDF5 file; the release generator independently refits the dark/light pixel distributions and recomputes ΔSP. |

The acquisition order, plate/well mapping, nominal compositions, identifiers, and pixel dimensions are also published in `data/processed/photokpfm_acquisition_log.csv` and Table S1. The 5×5 µm² scan area is stated in the current manuscript Methods; exact acquisition timestamps and an external absolute-potential calibration are not present in the HDF5 file and were not inferred.

This HDF5 file is not the binary-library dataset underlying rendered Figure 2. The binary raw data and original Figure 2 GP workflow were unavailable; `results/figures/author_supplied/Figure2_binary_photokpfm_author_supplied.jpeg` is therefore an explicitly scoped reference image rather than a reproducible output.

## Derived and literature data

`data/processed/initial_phase_kpfm_table.csv` is the checked downstream table from the endpoint-calibrated Colab phase-fitting workflow. The output-free notebook and both bundled raw files are supplied, but the fit remains separate from the local one-command route because it uses Colab upload widgets.

The literature data contain 15 retained complementary relationships from seven DOI-linked primary papers after source-level review. Publisher PDFs and article full text are not redistributed. The evidence graph contains 97 nodes and 120 edges; every edge retains `causal_status = complementary_only`.

Working directories, cached language-model responses, noncurrent plots, reviewer markup, and duplicate files are intentionally excluded. File-level SHA-256 values are recorded in `MANIFEST.csv`; generated table and ternary-model provenance records carry their own source hashes.
