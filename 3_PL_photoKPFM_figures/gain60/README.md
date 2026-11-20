# Gain-60 PL and nominal photoKPFM comparison

This directory contains the remade gain-60 top-read PL measurements and the six current individual panels: Figures 4a, 4b, S4, 5a, 5b and S5. Each panel is available in PNG, editable SVG and PDF under `figures/`. The plots contain their scientific axes, scales and legends; panel assembly is separate.

From `3_PL_photoKPFM_figures`, with the parent requirements installed:

```bash
python gain60/reproduce.py --validate-only
python gain60/reproduce.py --output reproduced_gain60
```

The destination must not already exist. The second command regenerates the source tables, statistics and all 18 figure files. It reads only the included public measurements and the shared parent inputs; no network, notebook execution or original local workbook path is required. `validate_sources()` is also an importable, read-only validation entrypoint.

## Measurements and sources

`data/top_spectra_raw.npz` retains all 22 top acquisitions, 96 wells and 416 measured emission-center channels from 455 to 870 nm, including the original literal cell strings, numeric/type arrays, acquisition/read indices and workbook-cell indices. Load it with `numpy.load(..., allow_pickle=False)`. Its peak and normalized arrays are checked caches; the public analysis recomputes both from the raw numeric cells. `data/measurement_metadata.json` records gain 60, excitation center 405 nm with 20 nm bandwidth, emission sampling at 1 nm with 20 nm optical bandwidth, and the original workbook SHA-256. The workbook itself is not required or redistributed here.

`data/well_map.csv` is the verified map for the current PL run. It has 96 well records and 91 distinct nominal recipes. The repeated wells remain individual observations. This map differs from the historical display lookup at C9, D9 and H10; those wells are outside the 13-well photoKPFM subset. Historical gain-100 processed tables remain an archive, and their displayed PL values are not substituted into this analysis.

The shared `../data/raw/photokpfm_measurements.h5` supplies 19 dark/light scan pairs across 13 nominal compositions. The shared `../data/processed/matched_pl_photokpfm_metrics.csv`, `../data/processed/initial_phase_kpfm_table.csv` and `../scripts/experimental/data_processing.py` retain the original Gaussian estimators and fitting method. Validation refits all 19 raw pairs and checks the archived panel values and error bars. No shared raw file or historical table is duplicated here. Input hashes and relative paths are recorded in `provenance.json`.

## Quantities and display

The PL quantity is the largest observed value across all 396 recorded channels from 455 through 850 nm, inclusive. A spectrum containing literal `OVRFLW`, a missing cell, an unrecognized cell or a nonfinite value within this band is invalid; no interpolation, smoothing or baseline subtraction is applied. The committed run has 2,112 valid maxima. Acquisition 22 is raw top-read block 64. Retention is `I22 / max(I1, ..., I22)` for a complete valid 22-acquisition well. Tied or edge-channel maxima and low counts are not silently excluded. No ADC ceiling, elapsed time or count conversion between gains is assumed.

Figure 4 retains all 96 solid markers at their 91 recipe positions. At the five shared positions, column 12 is drawn behind column 11. Figure 5a encodes the matched-subset PL intensity by marker area `36 + 180 * normalized_PL` in points squared; its legend uses the same transform.

Figures 5a and 5b retain `abs(mean(signed Gaussian shift))` at the original precision. Figure S5 retains `mean(abs(Gaussian shift))` rounded to six decimal places. These two historical summaries differ, particularly for G10. Figure S5 preserves the original `delta_SP_se` values: Gaussian fit-center variance terms are combined over scan pairs and divided by the number of acquisitions. They are not newly estimated PL errors or between-device replicate variability.

## Statistical scope

`data/matched13_panel_source.csv` and `data/statistics.json` are regenerated with the panels. The S5 Pearson coefficient is 0.696289, with analytic two-sided p = 0.008197 and paired-bootstrap 95% interval [-0.165625, 0.915069]. The bootstrap uses 10,000 nominal-composition-row resamples and NumPy PCG64 seed 20260816. The Spearman coefficient is 0.362637 (p = 0.223316). Partial Spearman, controlling ranked nominal PEA and BDA percentages, is -0.131190 (p = 0.700617; 9 degrees of freedom). The interval and differing tests should be considered together; these are exploratory, noncausal associations.

The PL and photoKPFM libraries were prepared independently and matched by nominal composition. The 19 scan pairs include repeated locations within droplets and are not 19 independent devices. The new PL run does not recover historical overflow values, establish same-specimen or synchronized measurements, or supply new spectral-component fits or structural phase fractions. No new 36-test family or multiple-testing q values are reported; historical gain-100 q values are not gain-60 results.

`validation.json` records the source, array, Gaussian and rendering checks. The six committed PNGs reproduce the approved individual exports exactly in the recorded environment. Plot geometry and numeric inputs are reproducible; font availability and library versions can affect rendering on another platform. All files needed for this route are within this directory and its shared parent data/code tree.
