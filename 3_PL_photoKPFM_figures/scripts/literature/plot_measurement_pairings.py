from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "data" / "literature" / "table_s3_relationship_summary.csv"
DEFAULT_OUTPUT = ROOT / "results" / "figures" / "supplementary" / "figure_literature_measurement_pairings_corrected"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot retained KPFM-PL literature pairings.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def save_figure(figure: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    frame = pd.read_csv(args.input)
    required = {
        "kpfm_observable_class",
        "pl_family_observable_class",
        "retained_relationships",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing summary columns: {missing}")
    frame["retained_relationships"] = pd.to_numeric(frame["retained_relationships"], errors="raise").astype(int)
    if int(frame["retained_relationships"].sum()) != 15:
        raise ValueError("The corrected pairing summary must contain 15 retained relationships.")

    kpfm_order = ["Contact potential difference", "Surface potential", "Work function"]
    optical_order = ["PL intensity / PLQY", "Carrier lifetime / TRPL"]
    matrix = (
        frame.pivot_table(
            index="kpfm_observable_class",
            columns="pl_family_observable_class",
            values="retained_relationships",
            aggfunc="sum",
            fill_value=0,
        )
        .reindex(index=kpfm_order, columns=optical_order, fill_value=0)
        .to_numpy(dtype=int)
    )

    figure, axis = plt.subplots(figsize=(5.8, 3.45))
    image = axis.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=max(4, int(matrix.max())), aspect="auto")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = int(matrix[row, column])
            axis.text(
                column,
                row,
                str(value),
                ha="center",
                va="center",
                fontsize=13,
                fontweight="bold",
                color="white" if value >= 3 else "#1D252C",
            )
    axis.set_xticks(np.arange(len(optical_order)), ["Steady-state PL /\nPLQY", "TRPL / carrier lifetime\n(literature context)"])
    axis.set_yticks(np.arange(len(kpfm_order)), kpfm_order)
    axis.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False, length=0, labelsize=9.0, pad=7)
    axis.tick_params(axis="y", length=0, labelsize=9.5, pad=7)
    axis.set_xticks(np.arange(-0.5, len(optical_order), 1), minor=True)
    axis.set_yticks(np.arange(-0.5, len(kpfm_order), 1), minor=True)
    axis.grid(which="minor", color="white", linewidth=1.5)
    axis.tick_params(which="minor", bottom=False, left=False)
    for spine in axis.spines.values():
        spine.set_visible(False)
    colorbar = figure.colorbar(image, ax=axis, fraction=0.052, pad=0.045, ticks=range(0, max(4, int(matrix.max())) + 1))
    colorbar.set_label("Retained relationships", fontsize=9.5)
    colorbar.ax.tick_params(labelsize=8.5)
    figure.subplots_adjust(left=0.30, right=0.90, top=0.74, bottom=0.08)
    save_figure(figure, args.output_stem)


if __name__ == "__main__":
    main()
