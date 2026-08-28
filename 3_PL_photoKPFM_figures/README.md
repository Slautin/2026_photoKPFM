# PL/photoKPFM manuscript figure scripts

This folder contains the Python scripts used to generate the PL/photoKPFM experimental and literature-context figures contributed to the manuscript. Generated images and manuscript data are intentionally not duplicated here.

## Scripts

- `scripts/publication/generate_publication_figures.py` generates experimental Figures 4A-4C, Figures 5A-5B, and the associated supplementary PL/photoKPFM figures.
- `scripts/literature/build_evidence_graph.py` constructs the auditable record-level evidence graph from the retained literature relationships.
- `scripts/literature/plot_evidence_graph.py` generates the corrected literature evidence-graph figure.
- `scripts/literature/plot_measurement_pairings.py` generates the corrected literature measurement-pairing matrix.

## Environment

```powershell
python -m venv .venv
& ./.venv/Scripts/Activate.ps1
python -m pip install -r requirements.txt
```

## Experimental figures

The publication generator expects a data directory containing a `processed` subdirectory with the curated manuscript plot tables.

```powershell
python scripts/publication/generate_publication_figures.py `
  --data-dir PATH_TO_MANUSCRIPT_DATA `
  --out-dir reproduced_figures `
  --overwrite
```

The script writes PNG, SVG, and PDF versions of each figure.

## Literature figures

Each literature script exposes explicit input and output arguments:

```powershell
python scripts/literature/build_evidence_graph.py --help
python scripts/literature/plot_evidence_graph.py --help
python scripts/literature/plot_measurement_pairings.py --help
```

The literature relationships describe evidence-linked complementarity and do not establish causality. PL and photoKPFM measurements are matched by nominal composition rather than spatially registered fields of view.
