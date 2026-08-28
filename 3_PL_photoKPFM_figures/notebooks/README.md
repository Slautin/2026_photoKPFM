# PL phase-fitting notebook

`pl_phase_fitting_colab.ipynb` is the output-free Google Colab workflow used for endpoint-calibrated, residual-aware PL emissive-component fitting.

1. Open the notebook in Colab.
2. Upload `data/raw/pl_plate_reader_export.csv` when prompted.
3. Upload `data/raw/photokpfm_measurements.h5` when prompted for the HDF5 file.
4. Run all cells and retain the generated phase-fit tables and diagnostics.

The notebook contains no credentials and no required private download link. Its output cells and execution counters are cleared in this release. The exact downstream table used by the local association analysis is already included as `data/processed/initial_phase_kpfm_table.csv`.

The fitted emissive components are PL-derived spectral descriptors. They are not XRD-confirmed structural phase fractions. This Colab step is intentionally not hidden inside `scripts/reproduce_all.py`; the separation makes the interactive input boundary explicit.

The notebook saves one scientific panel per file. Each well's PL fit and dark/illuminated KPFM distribution is separate, and the phase-distribution, KPFM-response, and dominant-peak diagnostics are standalone figures. Manuscript composites and panel lettering are not produced in the notebook.
