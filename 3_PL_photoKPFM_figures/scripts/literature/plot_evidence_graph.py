from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_NODES = ROOT / "data" / "literature" / "knowledge_graph_node_counts_corrected.csv"
DEFAULT_EDGES = ROOT / "data" / "literature" / "knowledge_graph_edge_counts_corrected.csv"
DEFAULT_RELATIONSHIPS = ROOT / "data" / "literature" / "table_s4_retained_relationships.csv"
DEFAULT_OUTPUT = ROOT / "results" / "figures" / "main" / "figure_photokpfm_knowledge_graph_corrected"

COLORS = {"teal": "#2A9D8F", "dark": "#252A31"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot the corrected KPFM-PL literature evidence graph.")
    parser.add_argument("--node-counts", type=Path, default=DEFAULT_NODES)
    parser.add_argument("--edge-counts", type=Path, default=DEFAULT_EDGES)
    parser.add_argument("--relationships", type=Path, default=DEFAULT_RELATIONSHIPS)
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8.0,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.4,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def save_figure(figure: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_knowledge_graph(node_counts: pd.DataFrame, edge_counts: pd.DataFrame, stem: Path) -> None:
    positions = {
        "Paper": (0.065, 0.55),
        "SampleContext": (0.215, 0.55),
        "KPFMObservation": (0.405, 0.73),
        "SteadyStatePLObservation": (0.405, 0.45),
        "TRPLObservation": (0.405, 0.15),
        "ScientificClaim": (0.665, 0.55),
        "MechanismConcept": (0.895, 0.67),
        "Evidence": (0.895, 0.35),
    }
    display = {
        "Paper": "Paper",
        "SampleContext": "Sample context",
        "KPFMObservation": "KPFM observations",
        "SteadyStatePLObservation": "Steady-state PL / PLQY",
        "TRPLObservation": "TRPL / carrier lifetime\nliterature context",
        "ScientificClaim": "Scientific claims",
        "MechanismConcept": "Proposed literature\ninterpretations",
        "Evidence": "Evidence excerpts",
    }
    colors = {
        "Paper": "#4C78A8",
        "SampleContext": "#72B7B2",
        "KPFMObservation": "#E45756",
        "SteadyStatePLObservation": "#F2CF5B",
        "TRPLObservation": "#FFF4CC",
        "ScientificClaim": "#B279A2",
        "MechanismConcept": "#59A14F",
        "Evidence": "#9D9D9D",
    }
    dimensions = {
        "Paper": (0.105, 0.115),
        "SampleContext": (0.135, 0.115),
        "KPFMObservation": (0.155, 0.115),
        "SteadyStatePLObservation": (0.175, 0.125),
        "TRPLObservation": (0.185, 0.135),
        "ScientificClaim": (0.150, 0.115),
        "MechanismConcept": (0.190, 0.135),
        "Evidence": (0.155, 0.115),
    }
    counts = node_counts.set_index("node_type")["count"].astype(int).to_dict()
    lookup = {
        (str(row.source_type), str(row.relationship), str(row.target_type)): int(row["count"])
        for _, row in edge_counts.iterrows()
    }
    figure, axis = plt.subplots(figsize=(7.4, 3.85))

    def count(source: str, relationship: str, target: str) -> int:
        return lookup.get((source, relationship, target), 0)

    def boundary(node_type: str, side: str) -> tuple[float, float]:
        xpos, ypos = positions[node_type]
        width, height = dimensions[node_type]
        offsets = {
            "left": (-width / 2, 0),
            "right": (width / 2, 0),
            "top": (0, height / 2),
            "bottom": (0, -height / 2),
            "upper_left": (-width / 2, height * 0.24),
            "upper_right": (width / 2, height * 0.24),
            "lower_left": (-width / 2, -height * 0.24),
            "lower_right": (width / 2, -height * 0.24),
        }
        dx, dy = offsets[side]
        return xpos + dx, ypos + dy

    def edge(
        source: str,
        target: str,
        relationship: str,
        source_side: str,
        target_side: str,
        label: str | None = None,
        label_xy: tuple[float, float] | None = None,
        dashed: bool = False,
        color: str = "#5C6168",
        curve: float = 0,
        arrowstyle: str = "-|>",
    ) -> None:
        value = count(source, relationship, target)
        if not value:
            return
        axis.add_patch(
            FancyArrowPatch(
                boundary(source, source_side),
                boundary(target, target_side),
                arrowstyle=arrowstyle,
                mutation_scale=8.5,
                linewidth=1.05 if color != COLORS["teal"] else 1.8,
                color=color,
                linestyle=(0, (4, 3)) if dashed else "-",
                connectionstyle=f"arc3,rad={curve}",
                shrinkA=1.5,
                shrinkB=1.5,
                zorder=1,
            )
        )
        if label and label_xy:
            axis.text(
                label_xy[0],
                label_xy[1],
                f"{label} ({value})",
                ha="center",
                va="center",
                fontsize=7.4,
                color="#42474E",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.97, "pad": 0.7},
                zorder=3,
            )

    edge("Paper", "SampleContext", "REPORTS_SAMPLE", "right", "left", "reports", (0.140, 0.635))
    edge("SampleContext", "KPFMObservation", "HAS_KPFM_OBSERVATION", "upper_right", "left", "has KPFM", (0.305, 0.815))
    edge("SampleContext", "SteadyStatePLObservation", "HAS_PL_OBSERVATION", "lower_right", "left", "has PL", (0.300, 0.365))
    edge("SampleContext", "TRPLObservation", "HAS_PL_OBSERVATION", "bottom", "left", "reported TRPL", (0.265, 0.285), dashed=True, curve=-0.12)
    edge("ScientificClaim", "KPFMObservation", "LINKS_KPFM", "upper_left", "right", "claim links", (0.555, 0.690), arrowstyle="<-")
    edge("ScientificClaim", "SteadyStatePLObservation", "LINKS_PL", "lower_left", "right", "claim links", (0.555, 0.475), arrowstyle="<-")
    edge("ScientificClaim", "TRPLObservation", "LINKS_PL", "bottom", "right", "broader context", (0.575, 0.285), dashed=True, curve=0.10, arrowstyle="<-")
    edge("KPFMObservation", "SteadyStatePLObservation", "COMPLEMENTS", "bottom", "top", "complementary links", (0.405, 0.590), color=COLORS["teal"], arrowstyle="<->")
    edge("ScientificClaim", "MechanismConcept", "INTERPRETED_AS", "upper_right", "left", "proposed as", (0.785, 0.765))
    edge("ScientificClaim", "Evidence", "SUPPORTED_BY", "lower_right", "left", "supports", (0.800, 0.430), arrowstyle="<-")

    for node_type, (xpos, ypos) in positions.items():
        width, height = dimensions[node_type]
        context = node_type == "TRPLObservation"
        axis.add_patch(
            Rectangle(
                (xpos - width / 2, ypos - height / 2),
                width,
                height,
                facecolor=colors[node_type],
                edgecolor="#9A7415" if context else "white",
                linewidth=1.15 if context else 0.75,
                linestyle=(0, (4, 2)) if context else "-",
                zorder=4,
            )
        )
        dark_text = node_type in {"SampleContext", "SteadyStatePLObservation", "TRPLObservation"}
        axis.text(
            xpos,
            ypos,
            f"{display[node_type]}\n(n={int(counts.get(node_type, 0))})",
            ha="center",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            color=COLORS["dark"] if dark_text else "white",
            zorder=5,
        )

    solid = Line2D([], [], color="#5C6168", linewidth=1.1, label="Evidence-linked path")
    complementary = Line2D([], [], color=COLORS["teal"], linewidth=1.8, label="Reported complementarity")
    dashed = Line2D([], [], color="#5C6168", linewidth=1.1, linestyle=(0, (4, 3)), label="Literature context only")
    axis.legend(
        handles=[solid, complementary, dashed],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.025),
        ncol=3,
        frameon=False,
        fontsize=7.6,
        handlelength=2.4,
        columnspacing=1.35,
    )
    axis.set_xlim(0.00, 1.00)
    axis.set_ylim(0.035, 0.88)
    axis.axis("off")
    figure.subplots_adjust(left=0.015, right=0.985, bottom=0.14, top=0.985)
    save_figure(figure, stem)


def verify_svg(svg_path: Path) -> None:
    root = ET.parse(svg_path).getroot()
    visible_text = " ".join(text.strip() for text in root.itertext() if text.strip())
    if re.search(r"\b16\b", visible_text):
        raise ValueError("Corrected graph still contains a visible retained-relationship count of 16.")
    for required in ("Scientific claims", "(n=15)", "claim links (15)"):
        if required not in visible_text:
            raise ValueError(f"Corrected graph is missing required visible text: {required}")


def main() -> None:
    args = parse_args()
    relationships = pd.read_csv(args.relationships)
    if len(relationships) != 15:
        raise ValueError("The corrected source table must contain 15 retained relationships.")
    if set(relationships["causal_status"].astype(str)) != {"complementary_only"}:
        raise ValueError("All retained relationships must remain complementary-only.")
    node_counts = pd.read_csv(args.node_counts)
    edge_counts = pd.read_csv(args.edge_counts)
    scientific_claims = int(node_counts.loc[node_counts["node_type"].eq("ScientificClaim"), "count"].iloc[0])
    if scientific_claims != 15:
        raise ValueError("ScientificClaim node count must be 15.")
    stale_edges = edge_counts.loc[pd.to_numeric(edge_counts["count"], errors="coerce").eq(16)]
    if not stale_edges.empty:
        raise ValueError("Corrected edge counts still contain a retained-relationship count of 16.")
    configure_plotting()
    plot_knowledge_graph(node_counts, edge_counts, args.output_stem)
    verify_svg(args.output_stem.with_suffix(".svg"))
    print(args.output_stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
