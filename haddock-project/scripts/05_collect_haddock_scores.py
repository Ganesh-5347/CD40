#!/usr/bin/env python3
"""
Stage 5: HADDOCK 3 Score Collection

Collects HADDOCK 3 scores from capri_ss.tsv output files.
Parses numbered module directories for energy metrics.

Outputs:
    - refine_scores.csv: Best model (rank 1) per candidate/source
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect HADDOCK 3 scores from refinement outputs."
    )
    parser.add_argument(
        "--run-list",
        default="haddock-project/outputs/haddock_refine/run_dirs.txt",
        help="File containing refinement run directories.",
    )
    parser.add_argument(
        "--out-dir",
        default="haddock-project/outputs/reports",
        help="Output directory for score summaries.",
    )
    return parser.parse_args()


def find_capri_ss(run_dir: Path) -> Path | None:
    """Find the last capri_ss.tsv in a HADDOCK 3 run directory."""
    haddock_run = run_dir / "run"
    if not haddock_run.exists():
        return None

    # Find all caprieval dirs, take the last one (highest number = after refinement)
    capri_dirs = sorted(haddock_run.glob("*_caprieval"))
    if not capri_dirs:
        return None

    capri_ss = capri_dirs[-1] / "capri_ss.tsv"
    return capri_ss if capri_ss.exists() else None


def find_model_pdb(run_dir: Path, model_name: str) -> Path | None:
    """Resolve model PDB path from capri_ss model column."""
    # capri_ss model column contains relative paths like ../2_mdref/mdref_1.pdb
    # Try to resolve relative to the caprieval dir or run dir
    haddock_run = run_dir / "run"

    # Try direct path from run dir
    candidate = haddock_run / model_name.lstrip("../")
    if candidate.exists():
        return candidate

    # Try searching emref dirs for the filename
    basename = Path(model_name).name
    for emref_dir in haddock_run.glob("*_emref"):
        pdb = emref_dir / basename
        if pdb.exists():
            return pdb
        gz_pdb = emref_dir / (basename + ".gz")
        if gz_pdb.exists():
            return gz_pdb

    # Try searching mdref dirs for the filename
    for mdref_dir in haddock_run.glob("*_mdref"):
        pdb = mdref_dir / basename
        if pdb.exists():
            return pdb
        gz_pdb = mdref_dir / (basename + ".gz")
        if gz_pdb.exists():
            return gz_pdb

    # Try searching flexref dirs
    for flexref_dir in haddock_run.glob("*_flexref"):
        pdb = flexref_dir / basename
        if pdb.exists():
            return pdb
        gz_pdb = flexref_dir / (basename + ".gz")
        if gz_pdb.exists():
            return gz_pdb

    return None


def parse_capri_ss(capri_ss_path: Path, run_dir: Path) -> list[dict]:
    """Parse capri_ss.tsv and extract scoring metrics."""
    rows = []
    with capri_ss_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            model_name = row.get("model", "")
            model_pdb = find_model_pdb(run_dir, model_name)

            metrics = {
                "model_path": str(model_pdb) if model_pdb else model_name,
                "haddock_score": _float(row.get("score")),
                "evdw": _float(row.get("vdw")),
                "eelec": _float(row.get("elec")),
                "edesolv": _float(row.get("desolv")),
                "bsa": _float(row.get("bsa")),
                "eair": _float(row.get("air")),
            }
            rows.append(metrics)
    return rows


def _float(val) -> float | None:
    if val is None or val == "" or val == "nan":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def best_model(models: list[dict]):
    """Select best model by HADDOCK score (most negative)."""
    scored = [m for m in models if m.get("haddock_score") is not None]
    if not scored:
        return None
    return min(scored, key=lambda m: m["haddock_score"])


def main():
    args = parse_args()
    run_list = Path(args.run_list)
    if not run_list.exists():
        raise FileNotFoundError(f"Run list not found: {run_list}")

    rows = []

    for line in run_list.read_text().splitlines():
        run_dir = Path(line.strip())
        if not run_dir.exists():
            continue

        candidate = run_dir.parent.name
        source = run_dir.name

        capri_ss = find_capri_ss(run_dir)
        if not capri_ss:
            continue

        models = parse_capri_ss(capri_ss, run_dir)
        for m in models:
            m["candidate_id"] = candidate
            m["source"] = source

        best = best_model(models)
        if best:
            rows.append(best)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "candidate_id", "source", "model_path", "haddock_score",
        "evdw", "eelec", "edesolv", "bsa",
    ]

    out_path = out_dir / "refine_scores.csv"
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out_path} ({len(rows)} models)")


if __name__ == "__main__":
    main()
