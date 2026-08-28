# Public release checklist

Do not describe a newly assembled archive as validated until all of the following are true:

- `python scripts/validate_release.py` exits with status 0 against the exact directory to be deposited.
- `MANIFEST.csv` was refreshed only after the final file set was frozen, then a second read-only validation passed.
- `docs/MANUSCRIPT_ASSET_INDEX.csv` resolves every manuscript-cited figure and Table S1–S4 to an existing release-relative file.
- `python scripts/reproduce_all.py` completes in a clean extraction and reports `"status": "passed"` in `reproduced_results/reproduction_summary.json`.
- Figure 1 and Figure 2 remain labeled `author-supplied-with-explicit-scope`; neither is presented as computationally reproduced.
- Generated manuscript assets remain single-panel outputs; composite layout and panel lettering are performed outside plotting scripts.
- The repository owner has selected an appropriate software/data license and added final citation metadata or an archival DOI.

No license, author list, or DOI is inferred by the build or validation scripts.
