#!/usr/bin/env python3
"""
Stage 5c: Merge HADDOCK and Rosetta scores — best model per candidate.

Reads:
    - refine_scores.csv      (HADDOCK energetics, 48 rows)
    - rosetta_interface.csv   (Rosetta characterisation, 48 rows)

For each candidate, selects the best source model by HADDOCK score
(most negative), then merges Rosetta characterisation columns.

Outputs:
    - combined_scores.csv     (one row per candidate, 21 rows)
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


HADDOCK_COLS = [
    "haddock_score",
]

ROSETTA_COLS = [
    "shape_complementarity",
    "packstat",
    "interface_hbonds",
    "interface_saltbridges",
    "nres_interface",
    "n_cdrs_engaged",
    "dominant_cdr",
    "h3_kink_type",
    "paratope_nres",
    "crd1_contacts",
    "crd2_contacts",
    "crd1_fraction",
    "crd2_fraction",
]

OUT_FIELDS = [
    "candidate_id",
    "source",
    *HADDOCK_COLS,
    *ROSETTA_COLS,
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge HADDOCK + Rosetta scores, best model per candidate.",
    )
    parser.add_argument(
        "--data-dir",
        default="haddock-project/outputs/reports",
        help="Directory containing refine_scores.csv and rosetta_interface.csv.",
    )
    parser.add_argument(
        "--out-dir",
        default="haddock-project/outputs/reports",
        help="Output directory for combined_scores.csv.",
    )
    return parser.parse_args()


def load_keyed(path: Path, key_cols: list[str]) -> dict[tuple, dict]:
    """Load a CSV into a dict keyed by (candidate_id, source)."""
    data = {}
    with path.open() as f:
        for row in csv.DictReader(f):
            key = tuple(row[k] for k in key_cols)
            data[key] = row
    return data


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    key_cols = ["candidate_id", "source"]

    # Load sources
    haddock_path = data_dir / "refine_scores.csv"
    rosetta_path = data_dir / "rosetta_interface.csv"

    if not haddock_path.exists():
        raise FileNotFoundError(f"Missing {haddock_path}")
    if not rosetta_path.exists():
        raise FileNotFoundError(f"Missing {rosetta_path}")

    haddock = load_keyed(haddock_path, key_cols)
    rosetta = load_keyed(rosetta_path, key_cols)

    # Merge all 48 rows first
    all_keys = sorted(set(haddock) | set(rosetta))
    merged = []
    for key in all_keys:
        row = dict(zip(key_cols, key))
        h = haddock.get(key, {})
        r = rosetta.get(key, {})
        for col in HADDOCK_COLS:
            row[col] = h.get(col, "NA")
        for col in ROSETTA_COLS:
            row[col] = r.get(col, "NA")
        merged.append(row)

    # Pick best source per candidate by Rosetta shape complementarity (highest)
    by_candidate: dict[str, list[dict]] = defaultdict(list)
    for row in merged:
        by_candidate[row["candidate_id"]].append(row)

    best_rows = []
    for cand in sorted(by_candidate):
        models = by_candidate[cand]
        scored = [m for m in models if m["shape_complementarity"] != "NA"]
        if scored:
            best = max(scored, key=lambda m: float(m["shape_complementarity"]))
        else:
            best = models[0]
        best_rows.append(best)

    out_path = out_dir / "combined_scores.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS)
        writer.writeheader()
        writer.writerows(best_rows)

    print(f"Wrote {out_path} ({len(best_rows)} candidates)")


if __name__ == "__main__":
    main()
