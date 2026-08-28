# 2026_photoKPFM

This repository contains the supporting workflows and reviewer-facing reproducibility materials for the FAPbI₃–BDAPbI₄–PEA₂PbI₄ combinatorial-library manuscript.

## Repository organization

- `1_photoKPFM_PEA-FAPI.ipynb`: binary-composition photoKPFM workflow.
- `2_photoKPFM_ternary_PEA-BDA-FAPI.ipynb`: ternary-composition photoKPFM workflow.
- `3_PL_photoKPFM_figures/`: complete PL/photoKPFM reproducibility release, including bundled source datasets, analysis and publication scripts, single-panel figure exports, supplementary tables, provenance, and validation records.

## Reproduce the PL/photoKPFM outputs

```powershell
cd 3_PL_photoKPFM_figures
py -3.11 -m venv .venv
& .\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/validate_release.py
python scripts/reproduce_all.py
```

See [`3_PL_photoKPFM_figures/README.md`](3_PL_photoKPFM_figures/README.md) for complete usage instructions, dataset descriptions, figure and table mappings, and scientific-scope limitations.
