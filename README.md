# PL and photoKPFM mapping of a perovskite composition library

This repository contains data, figure source tables, and analysis code for the
FAPbI3-BDAPbI4-PEA2PbI4 combinatorial library.

- [Gain-60 PL measurements and plots](3_PL_photoKPFM_figures/gain60/) contain the
  top-read spectra and the current Figures 4a, 4b, 5a, 5b, S4, and S5.
- [Figure S7](3_PL_photoKPFM_figures/gain60/figure_s7/) contains the five
  representative top-read PL spectra, source data, and plotting script.
- [Figure and table index](3_PL_photoKPFM_figures/docs/MANUSCRIPT_ASSET_INDEX.csv)
  maps each asset to its inputs and generating script.
- [Binary photoKPFM notebook](1_photoKPFM_PEA-FAPI.ipynb) and
  [ternary photoKPFM notebook](2_photoKPFM_ternary_PEA-BDA-FAPI.ipynb) document the
  corresponding experimental workflows.

From `3_PL_photoKPFM_figures/`, install the dependencies and run:

```text
python -m pip install -r requirements.txt
python scripts/validate_release.py
python scripts/reproduce_all.py
```

[Analysis and reproduction guide](3_PL_photoKPFM_figures/README.md)
