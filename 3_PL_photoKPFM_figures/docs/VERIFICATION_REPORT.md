# Verification

`python scripts/validate_release.py` checks the dataset identities, figure
exports, table provenance, asset index, and repository manifest.
`python scripts/reproduce_all.py` generates the current outputs in a fresh
directory and records the stage results in `reproduction_summary.json`.

The gain-60 analysis covers 96 wells, 91 distinct nominal compositions, and
22 complete top reads. It uses the 19 dark/illuminated AFM pairs to verify
the archived metrics for 13 nominal matches. Figure 5a has 13 distinct
marker locations and Figure S5 reuses the archived AFM fit errors.

The historical gain-100 Table S4 contains 36 exploratory tests. It belongs
to that run. The current gain-60 results report the peak-intensity estimates
and leave the 36-test multiplicity adjustment uncomputed.
