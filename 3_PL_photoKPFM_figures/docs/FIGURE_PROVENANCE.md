# Figure provenance

The paths below are release-relative. Every reproducible stem has PNG, SVG, and PDF exports. Figure 1 and Figure 2 are author-supplied rendered references and have the narrower scope stated in the table.

| Manuscript ID | Release output stem(s) | Status | Primary source data | Reproduction code |
|---|---|---|---|---|
| Figure 1 | `results/figures/author_supplied/Figure1_workflow_author_supplied.png` | Author-supplied reference | Rendered image extracted without modification from the current manuscript DOCX; editable schematic source unavailable. | None supplied. |
| Figure 2 | `results/figures/author_supplied/Figure2_binary_photokpfm_author_supplied.jpeg` | Author-supplied reference | Rendered image extracted without modification from the current manuscript DOCX; underlying binary-library measurements and original GP workflow unavailable. | None supplied. |
| Figure 3 | `results/figures/main/Figure3_GP_mean_map`; `Figure3_representative_SP_histogram_1`; `Figure3_representative_SP_histogram_2`; `Figure3_representative_SP_histogram_3`; `Figure3_representative_SP_histogram_4` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure 4a | `gain60/figures/Figure4a_gain60` | Reproducible | Observed gain-60 PL peak at top acquisition 22. | `gain60/reproduce.py` |
| Figure 4b | `gain60/figures/Figure4b_gain60` | Reproducible | Within-well gain-60 PL peak retention at acquisition 22. | `gain60/reproduce.py` |
| Figure 5a | `gain60/figures/Figure5a_gain60` | Reproducible | Gain-60 PL peaks and archived Gaussian photoKPFM shifts at 13 nominal matches. | `gain60/reproduce.py` |
| Figure 5b | `gain60/figures/Figure5b_gain60` | Reproducible | Gain-60 PL ranks and archived Gaussian photoKPFM ranks at 13 nominal matches. | `gain60/reproduce.py` |
| Figure S1 | `results/figures/supplementary/FigureS1_repeat_location_variability` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure S2 | `results/figures/supplementary/FigureS2_GP_uncertainty_map` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure S3 | `results/figures/supplementary/FigureS3_plate_composition_map` | Reproducible | `gain60/data/well_map.csv` | `scripts/publication/generate_publication_figures.py --plate-map-only` |
| Figure S4 | `gain60/figures/FigureS4_gain60` | Reproducible | Five gain-60 PL trajectories normalized by their respective maxima. | `gain60/reproduce.py` |
| Figure S5 | `gain60/figures/FigureS5_gain60` | Reproducible | Gain-60 peak intensity versus archived photoKPFM shifts and fit errors. | `gain60/reproduce.py` |
| Figure S6 | `results/figures/supplementary/FigureS6_evidence_graph` | Reproducible | `data/literature/evidence_graph_nodes.csv`; `data/literature/evidence_graph_edges.csv`; `data/literature/table_s4_retained_relationships.csv` | `scripts/literature/build_evidence_graph.py`; `scripts/literature/plot_evidence_graph.py` |
| Figure S7 | `gain60/figure_s7/FigureS7_gain60` | Reproducible | `gain60/figure_s7/FigureS7_source_spectra.csv`, five top-read spectra at acquisition 22 (Read 64) with display normalization over 455–850 nm. | `gain60/figure_s7/FigureS7_plot.py` |

The ternary generator also writes `results/source_data/ternary_photokpfm/photoKPFM_acquisitions_all_19.csv`, `Figure3_and_FigureS2_GP_grid_source_data.csv`, `Figure3_representative_histograms_source_data.csv`, `FigureS1_repeat_location_variability_source_data.csv`, `ternary_GP_model_and_figure_provenance.json`, and `ternary_photokpfm_output_manifest.csv`.

The evidence graph is built from 15 retained relationships. It has 97 nodes and 120 edges, all with `causal_status = complementary_only`. Its seven TRPL/carrier-lifetime relationships are literature context, not measurements from this experiment.

## Single-panel output policy

Plotting scripts emit only standalone scientific panels. A color bar is an ancillary scale axis, not another scientific panel. In particular, Figure 3 is represented by one GP mean map and four separate representative histograms; no plotting script combines them. Panel lettering, alignment, resizing, captions, and final multi-panel manuscript assembly occur in the manuscript-production software outside this repository. The PNG, SVG, and PDF variants of a stem are exports of the same panel, not independent results.

`docs/MANUSCRIPT_ASSET_INDEX.csv` is the machine-readable crosswalk for all manuscript-cited figures and Tables S1–S4.
