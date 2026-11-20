# Reproducibility guide

Run commands from `3_PL_photoKPFM_figures/`. Use Python 3.11 and install
the dependencies in `requirements.txt`. The analysis uses repository-relative
inputs and writes to a fresh output directory by default.

```text
python scripts/validate_release.py
python scripts/reproduce_all.py
```

The reproduction writes the gain-60 results under
`reproduced_results/results/gain60/`. The complete-file manifest records
integrity hashes for the supplied data and outputs. See the
[asset index](MANUSCRIPT_ASSET_INDEX.csv) for the individual source and
output paths. To reproduce only the six PL-related figures:

```text
python gain60/reproduce.py --output reproduced_results/gain60_only
```

## Other figure and table generators

```text
python scripts/publication/generate_publication_figures.py --plate-map-only --out-dir reproduced_results/plate_map
python scripts/publication/generate_ternary_photokpfm_panels.py --output-root reproduced_results/ternary
python scripts/publication/generate_supplementary_tables.py --output reproduced_results/tables/supplementary_tables_S1-S4.docx --overwrite
```

The supplementary table generator writes a 19-row acquisition CSV with
Table S1. Table S4 records the earlier gain-100 association analysis.

## Historical gain-100 workflow

```text
python scripts/reproduce_all.py --historical-gain100 --output-dir reproduced_results_gain100
```

The historical route regenerates its PL plots and the ancillary spectral
checks. Its `--skip-validation-analyses` option omits the two longer spectral
and PL metric checks. The gain-60 method uses acquisition numbers, full
measured spectra and the verified nominal plate map.

PNG, SVG, and PDF variants represent the same scientific plot. Manuscript
titles, lettering, and layout are applied after plotting.
