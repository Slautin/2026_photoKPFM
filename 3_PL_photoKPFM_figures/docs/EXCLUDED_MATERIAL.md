# Excluded material

The release omits exact duplicate scripts, superseded figures, presentation-only composites, temporary audits, user-specific paths, OAuth credentials, cached model outputs, publisher article full text, raw language-model responses, the full LiteratureAgent implementation, and unrelated PCE/stability-modeling files.

Figure 1 is supplied only as the rendered workflow image extracted from the current manuscript; its editable schematic source was unavailable. Figure 2 is supplied only as the rendered binary-photoKPFM image extracted from the manuscript; its underlying binary-library measurements and original GP generator were unavailable. These limitations are explicit in `FIGURE_PROVENANCE.md` and `MANUSCRIPT_ASSET_INDEX.csv`, and neither asset is claimed as computationally reproducible.

The excluded material does not remove data or code required for the analyses marked `reproducible` in the manuscript asset index. The bundled `data/raw/photokpfm_measurements.h5` is the ternary-library dataset used for Figure 3, Figure S1, Figure S2, and Table S1; it must not be described as the missing Figure 2 binary dataset.
