"""Render Figure S7 from the measured spectra in the adjacent CSV."""

import argparse
import csv
from itertools import groupby
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


SOURCE = Path(__file__).with_name("FigureS7_source_spectra.csv")
STEM = "FigureS7_gain60"
DISPLAY_MIN, DISPLAY_MAX = 455, 850
STYLES = {
    "A1": ("#0072B2", "-"),
    "F2": ("#E69F00", "--"),
    "G6": ("#00875F", "-"),
    "D3": ("#D55E00", "-."),
    "C11": ("#6F4B8B", "-"),
}


def read_spectra():
    with SOURCE.open(encoding="utf-8", newline="") as handle:
        rows = sorted(csv.DictReader(handle), key=lambda row: row["well_id"])
    spectra = {well: list(group) for well, group in groupby(rows, key=lambda row: row["well_id"])}
    if set(spectra) != set(STYLES):
        raise ValueError("The source CSV must contain exactly the five Figure S7 wells.")

    for well in STYLES:
        measured = sorted(spectra[well], key=lambda row: int(row["wavelength_nm"]))
        wavelengths = np.array([int(row["wavelength_nm"]) for row in measured])
        if not np.array_equal(wavelengths, np.arange(455, 871)):
            raise ValueError(f"{well}: expected every measured channel from 455 to 870 nm.")
        settings = {(row["acquisition_number"], row["read_block_number"], row["gain"]) for row in measured}
        if settings != {("22", "64", "60")}:
            raise ValueError(f"{well}: expected acquisition 22, Read 64, gain 60.")

        in_band = (wavelengths >= DISPLAY_MIN) & (wavelengths <= DISPLAY_MAX)
        raw = np.array([float(row["raw_intensity_counts"]) for row in measured])[in_band]
        peak = raw.max()
        if not np.isfinite(raw).all() or peak <= 0:
            raise ValueError(f"{well}: expected finite in-band counts with a positive maximum.")
        composition = "/".join(measured[0][f"nominal_{component}_pct"] for component in ("FA", "BDA", "PEA"))
        yield well, composition, wavelengths[in_band], raw / peak, peak


def render(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    available_fonts = set(font_manager.get_font_names())
    font = next(name for name in ("Arial", "Liberation Sans", "DejaVu Sans") if name in available_fonts)
    with matplotlib.rc_context({
        "font.family": font,
        "font.size": 10,
        "axes.labelsize": 11.5,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "svg.hashsalt": STEM,
        "pdf.fonttype": 42,
        "path.simplify": False,
    }):
        figure, axis = plt.subplots(figsize=(7.1, 3.7), facecolor="white")
        figure.subplots_adjust(left=0.107, right=0.70, bottom=0.145, top=0.95)
        for well, composition, wavelengths, normalized, peak in read_spectra():
            color, linestyle = STYLES[well]
            axis.plot(
                wavelengths, normalized, color=color, linestyle=linestyle,
                linewidth=1.65, label=f"{well}  -  {composition}\n{peak:,.0f} counts",
            )
        axis.set(
            xlabel="Emission wavelength (nm)", ylabel="Normalized PL intensity (a.u.)",
            xlim=(DISPLAY_MIN, DISPLAY_MAX), ylim=(-0.022, 1.045),
            xticks=[455, 500, 550, 600, 650, 700, 750, 800, 850],
            yticks=np.linspace(0, 1, 6),
        )
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#E3E3E3", linewidth=0.55)
        axis.legend(
            title="Well  -  FA/BDA/PEA (%)\nRaw peak intensity",
            loc="center left", bbox_to_anchor=(1.04, 0.53), frameon=False,
            fontsize=9.5, title_fontsize=9.5, labelspacing=0.72,
            handlelength=2.3, handletextpad=0.8, borderaxespad=0,
        )
        figure.savefig(output_dir / f"{STEM}.svg", metadata={"Date": None})
        figure.savefig(output_dir / f"{STEM}.pdf", metadata={"CreationDate": None})
        figure.savefig(output_dir / f"{STEM}.png", dpi=600)
        plt.close(figure)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=SOURCE.resolve().parent)
    render(parser.parse_args().output_dir)
