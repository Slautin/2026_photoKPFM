from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "reproduced_results"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce the PL/photoKPFM figures and validation outputs."
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-validation-analyses", action="store_true")
    return parser.parse_args()


def remove_output(path: Path, overwrite: bool) -> None:
    if not path.exists():
        return
    if not overwrite:
        raise FileExistsError(f"Output exists: {path}. Use --overwrite to replace it.")
    resolved = path.resolve()
    root = ROOT.resolve()
    if resolved == root or root not in resolved.parents:
        raise RuntimeError(f"Refusing to remove an unsafe output path: {resolved}")
    shutil.rmtree(resolved)


def run_stage(name: str, command: list[str], root: Path) -> dict[str, object]:
    started = time.time()
    print(f"\n[{name}] {' '.join(command)}", flush=True)
    completed = subprocess.run(command, cwd=root, check=False)
    elapsed = round(time.time() - started, 3)
    if completed.returncode:
        raise RuntimeError(f"Stage {name} failed with exit code {completed.returncode}.")
    return {"stage": name, "elapsed_seconds": elapsed, "return_code": completed.returncode}


def csv_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.reader(handle)) - 1


def merge_tree(source: Path, destination: Path) -> None:
    shutil.copytree(source, destination, dirs_exist_ok=True)


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    remove_output(output, args.overwrite)
    workspace = output / "_workspace"
    workspace.mkdir(parents=True)
    shutil.copytree(ROOT / "data", workspace / "data")
    shutil.copytree(
        ROOT / "scripts",
        workspace / "scripts",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    merge_tree(
        ROOT / "results" / "figures" / "author_supplied",
        workspace / "results" / "figures" / "author_supplied",
    )
    python = sys.executable
    logs: list[dict[str, object]] = []
    publication_staging = workspace / "_publication_figures"

    logs.append(
        run_stage(
            "publication_figures",
            [
                python,
                "scripts/publication/generate_publication_figures.py",
                "--out-dir",
                str(publication_staging),
            ],
            workspace,
        )
    )
    merge_tree(publication_staging / "main", workspace / "results" / "figures" / "main")
    merge_tree(
        publication_staging / "supplementary",
        workspace / "results" / "figures" / "supplementary",
    )
    merge_tree(
        publication_staging / "source_data",
        workspace / "results" / "source_data" / "publication_figures",
    )
    merge_tree(
        publication_staging / "audit",
        workspace / "results" / "audit" / "publication_figures",
    )
    shutil.rmtree(publication_staging)
    logs.append(
        run_stage(
            "ternary_photokpfm_panels",
            [
                python,
                "scripts/publication/generate_ternary_photokpfm_panels.py",
                "--output-root",
                str(workspace / "results"),
            ],
            workspace,
        )
    )
    logs.append(
        run_stage(
            "supplementary_tables_S1_S4",
            [
                python,
                "scripts/publication/generate_supplementary_tables.py",
                "--output",
                str(
                    workspace
                    / "results"
                    / "tables"
                    / "supplementary_tables_current_manuscript_S1-S4.docx"
                ),
                "--overwrite",
            ],
            workspace,
        )
    )
    if not args.skip_validation_analyses:
        logs.append(
            run_stage(
                "pl_metric_validation",
                [python, "scripts/experimental/analyze_pl_metrics.py"],
                workspace,
            )
        )
        logs.append(
            run_stage(
                "spectral_processing_validation",
                [python, "scripts/experimental/validate_spectral_processing.py"],
                workspace,
            )
        )
    association_dir = workspace / "results" / "validation" / "association_analysis"
    logs.append(
        run_stage(
            "pl_photokpfm_associations",
            [
                python,
                "scripts/analysis/analyze_pl_photokpfm_associations.py",
                "--matched-csv",
                "data/processed/matched_pl_photokpfm_metrics.csv",
                "--phase-csv",
                "data/processed/initial_phase_kpfm_table.csv",
                "--literature-relationships",
                "data/literature/table_s4_retained_relationships.csv",
                "--out-dir",
                str(association_dir),
            ],
            workspace,
        )
    )
    graph_dir = workspace / "results" / "source_data" / "evidence_graph"
    logs.append(
        run_stage(
            "evidence_graph_data",
            [
                python,
                "scripts/literature/build_evidence_graph.py",
                "--nodes-out",
                str(graph_dir / "evidence_graph_nodes.csv"),
                "--edges-out",
                str(graph_dir / "evidence_graph_edges.csv"),
                "--summary-out",
                str(graph_dir / "evidence_graph_summary.json"),
            ],
            workspace,
        )
    )
    logs.append(
        run_stage(
            "literature_evidence_graph",
            [
                python,
                "scripts/literature/plot_evidence_graph.py",
                "--output-stem",
                str(
                    workspace
                    / "results"
                    / "figures"
                    / "supplementary"
                    / "FigureS6_evidence_graph"
                ),
            ],
            workspace,
        )
    )

    generated = workspace / "results"
    final_results = output / "results"
    shutil.move(str(generated), str(final_results))
    shutil.rmtree(workspace)
    summary = {
        "status": "passed",
        "source_counts": {
            "pl_library_rows": csv_rows(ROOT / "data" / "processed" / "pl_retention_metrics.csv"),
            "matched_pl_photokpfm_rows": csv_rows(
                ROOT / "data" / "processed" / "matched_pl_photokpfm_metrics.csv"
            ),
            "retained_literature_relationships": csv_rows(
                ROOT / "data" / "literature" / "table_s4_retained_relationships.csv"
            ),
        },
        "generated_files": sum(1 for path in final_results.rglob("*") if path.is_file()),
        "stages": logs,
    }
    (output / "reproduction_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"\nReproduction completed: {output}")


if __name__ == "__main__":
    main()
