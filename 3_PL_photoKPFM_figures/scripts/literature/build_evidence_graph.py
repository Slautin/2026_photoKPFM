from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RELATIONSHIPS = ROOT / "data" / "literature" / "table_s4_retained_relationships.csv"
DEFAULT_NODE_OUTPUT = ROOT / "data" / "literature" / "evidence_graph_nodes.csv"
DEFAULT_EDGE_OUTPUT = ROOT / "data" / "literature" / "evidence_graph_edges.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "data" / "literature" / "evidence_graph_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an auditable KPFM-PL evidence graph.")
    parser.add_argument("--relationships", type=Path, default=DEFAULT_RELATIONSHIPS)
    parser.add_argument("--nodes-out", type=Path, default=DEFAULT_NODE_OUTPUT)
    parser.add_argument("--edges-out", type=Path, default=DEFAULT_EDGE_OUTPUT)
    parser.add_argument("--summary-out", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    return parser.parse_args()


def stable_id(prefix: str, *values: object) -> str:
    text = "|".join(str(value).strip() for value in values)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]
    return f"{prefix}:{digest}"


def add_node(
    rows: dict[str, dict[str, object]],
    node_type: str,
    label: str,
    doi: str,
    relationship_id: str,
    attributes: dict[str, object] | None = None,
) -> str:
    node_id = stable_id(node_type.lower(), doi, relationship_id, label)
    rows[node_id] = {
        "node_id": node_id,
        "node_type": node_type,
        "label": label,
        "doi": doi,
        "relationship_id": relationship_id,
        "attributes_json": json.dumps(attributes or {}, ensure_ascii=False, sort_keys=True),
    }
    return node_id


def add_edge(
    rows: list[dict[str, object]],
    source_id: str,
    relationship: str,
    target_id: str,
    source_row: pd.Series,
) -> None:
    relationship_id = str(source_row["relationship_id"])
    doi = str(source_row["doi"])
    edge_id = stable_id("edge", source_id, relationship, target_id, relationship_id)
    rows.append(
        {
            "edge_id": edge_id,
            "source_id": source_id,
            "relationship": relationship,
            "target_id": target_id,
            "relationship_id": relationship_id,
            "doi": doi,
            "evidence_strength": source_row["evidence_strength"],
            "causal_status": source_row["causal_status"],
        }
    )


def optical_type(text: str) -> str:
    lowered = text.lower()
    if any(term in lowered for term in ("trpl", "lifetime", "decay")):
        return "TRPLContextObservation"
    return "SteadyStatePLObservation"


def validate_relationships(frame: pd.DataFrame) -> None:
    required = {
        "relationship_id",
        "paper_title",
        "doi",
        "intervention_or_comparison",
        "kpfm_observable_and_direct_result",
        "pl_family_observable_and_direct_result",
        "evidence_grounded_claim",
        "complementary_interpretation",
        "interpretation_limit",
        "evidence_strength",
        "causal_status",
    }
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Missing relationship columns: {missing}")
    if frame.empty:
        raise ValueError("The corrected relationship table is empty.")
    if frame["relationship_id"].duplicated().any():
        raise ValueError("Relationship identifiers must be unique.")
    if frame["doi"].fillna("").str.strip().eq("").any():
        raise ValueError("Every retained relationship must have a DOI.")
    statuses = set(frame["causal_status"].fillna("").str.strip())
    if statuses != {"complementary_only"}:
        raise ValueError(f"Unexpected causal status values: {sorted(statuses)}")


def build_graph(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes: dict[str, dict[str, object]] = {}
    edges: list[dict[str, object]] = []
    paper_nodes: dict[str, str] = {}

    for _, row in frame.iterrows():
        relationship_id = str(row["relationship_id"])
        doi = str(row["doi"]).strip().lower()
        paper_id = paper_nodes.get(doi)
        if paper_id is None:
            paper_id = stable_id("paper", doi)
            nodes[paper_id] = {
                "node_id": paper_id,
                "node_type": "Paper",
                "label": str(row["paper_title"]),
                "doi": doi,
                "relationship_id": "",
                "attributes_json": json.dumps({"doi": doi}, ensure_ascii=False, sort_keys=True),
            }
            paper_nodes[doi] = paper_id

        comparison_id = add_node(
            nodes,
            "SampleComparison",
            str(row["intervention_or_comparison"]),
            doi,
            relationship_id,
        )
        kpfm_id = add_node(
            nodes,
            "KPFMObservation",
            str(row["kpfm_observable_and_direct_result"]),
            doi,
            relationship_id,
        )
        optical_label = str(row["pl_family_observable_and_direct_result"])
        optical_id = add_node(
            nodes,
            optical_type(optical_label),
            optical_label,
            doi,
            relationship_id,
            {"measured_in_present_study": optical_type(optical_label) == "SteadyStatePLObservation"},
        )
        claim_id = add_node(
            nodes,
            "ScientificClaim",
            str(row["evidence_grounded_claim"]),
            doi,
            relationship_id,
        )
        interpretation_id = add_node(
            nodes,
            "ProposedLiteratureInterpretation",
            str(row["complementary_interpretation"]),
            doi,
            relationship_id,
            {"interpretation_limit": str(row["interpretation_limit"])},
        )
        evidence_id = add_node(
            nodes,
            "EvidenceRecord",
            relationship_id,
            doi,
            relationship_id,
            {"source_table": "Table S4"},
        )

        add_edge(edges, paper_id, "REPORTS_COMPARISON", comparison_id, row)
        add_edge(edges, comparison_id, "HAS_KPFM_OBSERVATION", kpfm_id, row)
        add_edge(edges, comparison_id, "HAS_OPTICAL_OBSERVATION", optical_id, row)
        add_edge(edges, kpfm_id, "COMPLEMENTS", optical_id, row)
        add_edge(edges, kpfm_id, "SUPPORTS_CLAIM", claim_id, row)
        add_edge(edges, optical_id, "SUPPORTS_CLAIM", claim_id, row)
        add_edge(edges, claim_id, "INTERPRETED_AS", interpretation_id, row)
        add_edge(edges, evidence_id, "SUPPORTS", claim_id, row)

    node_frame = pd.DataFrame(nodes.values()).sort_values(["node_type", "node_id"])
    edge_frame = pd.DataFrame(edges).drop_duplicates("edge_id").sort_values("edge_id")
    return node_frame.reset_index(drop=True), edge_frame.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    relationships = pd.read_csv(args.relationships)
    validate_relationships(relationships)
    nodes, edges = build_graph(relationships)
    for destination in (args.nodes_out, args.edges_out, args.summary_out):
        destination.parent.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(args.nodes_out, index=False)
    edges.to_csv(args.edges_out, index=False)
    summary = {
        "papers": int(relationships["doi"].nunique()),
        "retained_relationships": int(len(relationships)),
        "nodes": int(len(nodes)),
        "edges": int(len(edges)),
        "causal_status": sorted(relationships["causal_status"].unique().tolist()),
        "steady_state_relationships": int(
            relationships["pl_family_observable_and_direct_result"].map(optical_type).eq("SteadyStatePLObservation").sum()
        ),
        "trpl_context_relationships": int(
            relationships["pl_family_observable_and_direct_result"].map(optical_type).eq("TRPLContextObservation").sum()
        ),
    }
    args.summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
