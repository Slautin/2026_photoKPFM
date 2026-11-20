# PL and photoKPFM analysis

The [gain-60 dataset](gain60/) supplies the current PL plots. It contains 22
top-read acquisitions of 96 wells, each measured from 455 to 870 nm.
Figures 4a, 4b, 5a, 5b, S4, and S5 use the measured peak within 455-850 nm.
Acquisition 22 is the comparison read. The workbook gives acquisition order;
exact elapsed time requires a separate acquisition log.

The [asset index](docs/MANUSCRIPT_ASSET_INDEX.csv) links each figure and table
to its data and script. The [data dictionary](docs/DATA_DICTIONARY.md) and
[data provenance](docs/DATA_PROVENANCE.md) describe the archived measurements.
The paired ternary photoKPFM data are in `data/raw/photokpfm_measurements.h5`.
The comparison uses 13 nominal compositions shared by the separately prepared
PL and photoKPFM libraries.

## Reproduce the outputs

Python 3.11 is recommended. From this directory:

```text
python -m venv .venv
```

Use `.venv/Scripts/python.exe` on Windows or `.venv/bin/python` on macOS/Linux
in place of `python` in these commands:

```text
python -m pip install -r requirements.txt
python scripts/validate_release.py
python scripts/reproduce_all.py
```

The default workflow generates the gain-60 plots, the ternary AFM plots, the
composition map, the literature evidence graph, and the supplementary tables
under `reproduced_results/`. For the gain-60 plots alone, see the
[gain-60 instructions](gain60/README.md). Figures 1 and 2 are rendered reference
images. The editable Figure 1 illustration and the binary data behind the
Figure 2 image are unavailable for computational reproduction.

## Historical gain-100 analysis

The earlier plate-reader export and derived files remain available to
reproduce their original analysis, including the historical 36-test Table S4.
Earlier Figure 4 rendering used a plotting map that misplaced five repeated
wells. The current gain-60 plots use the documented preparation map.
The gain-60 peak-intensity associations are in the gain-60 data directory.
A new 36-test correction would require reproducing the full spectral-component
analysis on the gain-60 measurements.

To rerun the historical analysis, use:

```text
python scripts/reproduce_all.py --historical-gain100 --output-dir reproduced_results_gain100
```

See the [reproducibility guide](docs/REPRODUCIBILITY_GUIDE.md) for the individual
commands and [scientific scope](docs/SCIENTIFIC_SCOPE.md) for interpretation.
