#!/usr/bin/env python3
"""
Stage 8: Final Ranking with CRD Domain Categories

Reads:
    - combined_scores.csv  (21 best models from 07, includes Rosetta CRD fractions)

Assigns CRD binding category (CRD1, CRD2, CRD1+CRD2, NEITHER) based on
post-refinement Rosetta contact fractions, then ranks by HADDOCK score.

Outputs:
    - reports/final_ranking.csv
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build final ranking with CRD domain categories.",
    )
    parser.add_argument(
        "--combined",
        default="haddock-project/outputs/reports/combined_scores.csv",
        help="Combined scores CSV from 07_select_best_models.",
    )
    parser.add_argument(
        "--min-fraction",
        type=float,
        default=0.10,
        help="Minimum CRD fraction (from Rosetta) to count as binding a domain.",
    )
    parser.add_argument(
        "--out-dir",
        default="haddock-project/outputs/reports",
        help="Output directory.",
    )
    return parser.parse_args()


OUT_FIELDS = [
    "rank",
    "candidate_id",
    "source",
    "crd_category",
    "crd1_fraction",
    "crd2_fraction",
    "haddock_score",
    "shape_complementarity",
    "packstat",
    "interface_hbonds",
    "interface_saltbridges",
    "nres_interface",
    "n_cdrs_engaged",
    "dominant_cdr",
    "h3_kink_type",
    "paratope_nres",
]


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    combined_path = Path(args.combined)
    if not combined_path.exists():
        raise FileNotFoundError(f"Missing {combined_path}")

    rows = []
    with combined_path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    # Assign CRD category from Rosetta-computed fractions
    thresh = args.min_fraction
    for row in rows:
        c1 = float(row.get("crd1_fraction", 0))
        c2 = float(row.get("crd2_fraction", 0))
        c1_pass = c1 >= thresh
        c2_pass = c2 >= thresh
        if c1_pass and c2_pass:
            row["crd_category"] = "CRD1+CRD2"
        elif c1_pass:
            row["crd_category"] = "CRD1"
        elif c2_pass:
            row["crd_category"] = "CRD2"
        else:
            row["crd_category"] = "NEITHER"

    # Sort by HADDOCK score (most negative first)
    rows.sort(key=lambda r: float(r.get("haddock_score", 0)))

    for i, row in enumerate(rows, 1):
        row["rank"] = i

    # Write CSV
    out_path = out_dir / "final_ranking.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Print summary
    cats = defaultdict(int)
    for row in rows:
        cats[row["crd_category"]] += 1

    print(f"Wrote {out_path} ({len(rows)} candidates)")
    print(f"\nCRD categories (threshold={thresh}):")
    for cat in ["CRD1+CRD2", "CRD1", "CRD2", "NEITHER"]:
        if cats.get(cat, 0) > 0:
            print(f"  {cat}: {cats[cat]}")


if __name__ == "__main__":
    main()
