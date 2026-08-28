# Data dictionary

| Path | Rows or shape | Contents |
|---|---:|---|
| `data/raw/pl_plate_reader_export.csv` | 43,818 physical rows | Original nonrectangular Cytation 5 export. The analysis uses 22 complete top-read blocks × 96 wells. |
| `data/raw/photokpfm_measurements.h5` | 9 required datasets | Original ternary-library dark/illuminated photoKPFM acquisitions; detailed shapes are in `DATA_PROVENANCE.md`. |
| `data/processed/pl_retention_metrics.csv` | 96 | PL peak intensity and retention metrics for the full plate. |
| `data/processed/plate_composition_map.csv` | 96 | Well-to-nominal-composition lookup for the plate layout. |
| `data/processed/matched_pl_photokpfm_metrics.csv` | 13 | Nominally composition-matched PL/photoKPFM records used for exploratory comparisons. |
| `data/processed/initial_phase_kpfm_table.csv` | 13 | Initial PL spectral-component outputs joined to photoKPFM observables. |
| `data/processed/candidate_rank_comparison.csv` | 13 | Composition ranks for top-read PL peak intensity and absolute photoKPFM response. |
| `data/processed/representative_pl_trace_selection.csv` | 13 | Matched-well temporal classes and representative-panel selection flags. |
| `data/processed/representative_pl_traceability.csv` | 5 | Traceability records for the five representative PL wells. |
| `data/processed/representative_pl_traces.csv` | 110 | Five representative PL traces × 22 timepoints. |
| `data/processed/all_matched_pl_traces.csv` | 286 | Thirteen matched PL traces × 22 timepoints. |
| `data/processed/exploratory_association_statistics.csv` | 36 | Prespecified exploratory PL/photoKPFM association tests and multiplicity-adjusted results. |
| `data/processed/photokpfm_acquisition_log.csv` | 19 | Table S1 machine-readable acquisition order, well identity, nominal composition, file IDs, and scan dimensions. |
| `data/validation/` | multiple | Spectral-processing, baseline, time-course, and sensitivity checks. |
| `data/literature/table_s3_relationship_summary.csv` | 6 | Source table for current manuscript Table S2: counts by KPFM and optical-observable class. |
| `data/literature/table_s4_retained_relationships.csv` | 15 | Source table for current manuscript Table S3: source-audited complementary relationships from seven DOI-linked papers. |
| `data/literature/evidence_graph_nodes.csv` | 97 | Canonical evidence-graph nodes. |
| `data/literature/evidence_graph_edges.csv` | 120 | Canonical evidence-graph edges with DOI, evidence strength, and causal status. |
| `results/source_data/ternary_photokpfm/` | multiple | Figure 3/S1/S2 panel source data, GP model provenance, and output manifest. |
| `results/tables/supplementary_tables_current_manuscript_S1-S4.docx` | 19/6/15/36 data rows | Combined, manuscript-numbered Tables S1–S4. |

The historical filenames `table_s3_relationship_summary.csv` and `table_s4_retained_relationships.csv` are retained as stable machine-readable source names; in the current combined supplement their contents appear as Tables S2 and S3, respectively.

`causal_status = complementary_only` means the KPFM and optical observations concern the same treatment or comparison but do not establish a causal pathway. TRPL/carrier-lifetime rows are literature context and were not measured in this study.
