# Reproducibility verification protocol

The former dated build report described a superseded file set and is not evidence for this release. Verification must be performed against the exact extracted directory being reviewed.

## Current execution evidence

On 2026-08-28, a full isolated run of `scripts/reproduce_all.py` completed with exit status 0 in 237.6 seconds. All eight stages returned 0: publication figures, ternary photoKPFM panels, Tables S1-S4, PL-metric validation, spectral-processing validation, exploratory PL/photoKPFM associations, evidence-graph data, and Figure S6. The run reported 353 generated files. SHA-256 comparison found exact matches for all 15 comparable supplied and reproduced PNG files, including the author-supplied Figure 1 reference and every reproducible manuscript panel.

The read-only release validator also passed the raw-data, processed-data, graph, ternary provenance, manuscript asset-index, Python, notebook, and single-panel checks. The release must not be called frozen until the explicitly identified legacy artifacts are removed, `MANIFEST.csv` is refreshed, and the final read-only validation exits 0.

## Package validation

```text
python scripts/validate_release.py
```

Acceptance requires exit status 0 and `"status": "passed"` in the printed JSON. The validator checks raw and processed dimensions, manuscript assets, supplementary-table provenance, graph counts, single-panel source policy, notebook cleanliness, personal paths, manifest integrity, and exact duplicate content.

## Clean reproduction

```text
python scripts/reproduce_all.py
```

Acceptance requires `reproduced_results/reproduction_summary.json` to report `"status": "passed"`, every listed stage to return 0, and the expected outputs to exist under `reproduced_results/results/`. A second test can use `--overwrite` to verify deterministic replacement of only the generated reproduction directory.

## Scientific acceptance invariants

- 96 PL wells, 91 unique nominal compositions, and 22 complete top-read timepoints.
- 19 paired ternary photoKPFM acquisitions spanning 13 unique wells; dark and illuminated arrays are each 19×128×128.
- 15 retained complementary literature relationships from seven DOI-linked papers.
- Evidence graph: 97 nodes and 120 edges, all `complementary_only`.
- Tables S1–S4 contain 19, 6, 15, and 36 data rows, respectively.
- None of the 36 exploratory association tests passes BH-adjusted 0.05.
- Figure 1 and Figure 2 are present only as explicitly scoped author-supplied rendered references.

`MANIFEST.csv` provides file-level hashes for the frozen release. The supplementary-table provenance JSON and ternary GP provenance provide workflow-specific source and output records.
