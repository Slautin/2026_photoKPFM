# Figure provenance

The paths below are release-relative. Every reproducible stem has PNG, SVG, and PDF exports. Figure 1 and Figure 2 are author-supplied rendered references and have the narrower scope stated in the table.

| Manuscript ID | Release output stem(s) | Status | Primary source data | Reproduction code |
|---|---|---|---|---|
| Figure 1 | `results/figures/author_supplied/Figure1_workflow_author_supplied.png` | Author-supplied reference | Rendered image extracted without modification from the current manuscript DOCX; editable schematic source unavailable. | None supplied. |
| Figure 2 | `results/figures/author_supplied/Figure2_binary_photokpfm_author_supplied.jpeg` | Author-supplied reference | Rendered image extracted without modification from the current manuscript DOCX; underlying binary-library measurements and original GP workflow unavailable. | None supplied. |
| Figure 3 | `results/figures/main/Figure3_GP_mean_map`; `Figure3_representative_SP_histogram_1`; `Figure3_representative_SP_histogram_2`; `Figure3_representative_SP_histogram_3`; `Figure3_representative_SP_histogram_4` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure 4a | `results/figures/main/Figure4A_PL_189min_ternary` | Reproducible | `data/processed/pl_retention_metrics.csv`; `data/processed/matched_pl_photokpfm_metrics.csv`; `data/processed/representative_pl_trace_selection.csv` | `scripts/publication/generate_publication_figures.py` |
| Figure 4b | `results/figures/main/Figure4B_PL_retention_ternary` | Reproducible | `data/processed/pl_retention_metrics.csv`; `data/processed/matched_pl_photokpfm_metrics.csv`; `data/processed/representative_pl_trace_selection.csv` | `scripts/publication/generate_publication_figures.py` |
| Figure 5a | `results/figures/main/Figure5A_matched_screening_space` | Reproducible | `data/processed/matched_pl_photokpfm_metrics.csv`; `data/processed/representative_pl_trace_selection.csv` | `scripts/publication/generate_publication_figures.py` |
| Figure 5b | `results/figures/main/Figure5B_candidate_rank_comparison` | Reproducible | `data/processed/matched_pl_photokpfm_metrics.csv`; `data/processed/candidate_rank_comparison.csv` is the exported rank table regenerated from those metrics. | `scripts/publication/generate_publication_figures.py` |
| Figure S1 | `results/figures/supplementary/FigureS1_repeat_location_variability` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure S2 | `results/figures/supplementary/FigureS2_GP_uncertainty_map` | Reproducible | `data/raw/photokpfm_measurements.h5` | `scripts/publication/generate_ternary_photokpfm_panels.py` |
| Figure S3 | `results/figures/supplementary/FigureS3_plate_composition_map` | Reproducible | `data/processed/plate_composition_map.csv` | `scripts/publication/generate_publication_figures.py` |
| Figure S4 | `results/figures/supplementary/FigureS4_representative_PL_time_traces` | Reproducible | `data/processed/representative_pl_traces.csv`; `data/processed/representative_pl_trace_selection.csv` | `scripts/publication/generate_publication_figures.py` |
| Figure S6 | `results/figures/supplementary/FigureS6_evidence_graph` | Reproducible | `data/literature/evidence_graph_nodes.csv`; `data/literature/evidence_graph_edges.csv`; `data/literature/table_s4_retained_relationships.csv` | `scripts/literature/build_evidence_graph.py`; `scripts/literature/plot_evidence_graph.py` |

The ternary generator also writes `results/source_data/ternary_photokpfm/photoKPFM_acquisitions_all_19.csv`, `Figure3_and_FigureS2_GP_grid_source_data.csv`, `Figure3_representative_histograms_source_data.csv`, `FigureS1_repeat_location_variability_source_data.csv`, `ternary_GP_model_and_figure_provenance.json`, and `ternary_photokpfm_output_manifest.csv`.

The evidence graph is built from 15 retained relationships. It has 97 nodes and 120 edges, all with `causal_status = complementary_only`. Its seven TRPL/carrier-lifetime relationships are literature context, not measurements from this experiment.

## Single-panel output policy

Plotting scripts emit only standalone scientific panels. A color bar is an ancillary scale axis, not another scientific panel. In particular, Figure 3 is represented by one GP mean map and four separate representative histograms; no plotting script combines them. Panel lettering, alignment, resizing, captions, and final multi-panel manuscript assembly occur in the manuscript-production software outside this repository. The PNG, SVG, and PDF variants of a stem are exports of the same panel, not independent results.

`docs/MANUSCRIPT_ASSET_INDEX.csv` is the machine-readable crosswalk for all manuscript-cited figures and Tables S1–S4.
