#!/usr/bin/env python3
"""
Stage 2: Contact Filtering

Filters candidates by CRD1 contact fraction and categorizes by source agreement.

Outputs:
    - outputs/prescreen/prescreen_summary.csv
    - outputs/prescreen/prescreen_keep.txt
    - outputs/prescreen/source_analysis.csv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter candidates by CRD1 contact fraction."
    )
    parser.add_argument(
        "--contacts-csv",
        default="haddock-project/outputs/prescreen/prescreen_contacts.csv",
        help="Per-model contacts CSV from script 01.",
    )
    parser.add_argument(
        "--min-crd1-fraction",
        type=float,
        default=0.5,
        dest="min_crd1_fraction",
        help="Minimum CRD1 contact fraction to pass prescreen.",
    )
    parser.add_argument(
        "--min-sources",
        type=int,
        default=2,
        choices=[1, 2, 3],
        help="Minimum number of sources (AF2/AF3/Boltz) that must pass the fraction threshold (default: 2).",
    )
    parser.add_argument(
        "--output-dir",
        default="haddock-project/outputs/prescreen",
        help="Output directory for prescreen results.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    contacts_path = Path(args.contacts_csv)
    if not contacts_path.exists():
        raise FileNotFoundError(f"Missing contacts CSV: {contacts_path}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    min_frac = args.min_crd1_fraction

    # Read per-model contacts
    rows = []
    with contacts_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(row)

    # Group by candidate
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["candidate_id"]].append(row)

    # --- prescreen_summary.csv + prescreen_keep.txt ---
    summary_rows = []
    for candidate_id, items in sorted(grouped.items()):
        best = max(items, key=lambda r: float(r["crd1_fraction"]))
        models_total = len(items)
        models_passing = sum(
            1 for r in items if float(r["crd1_fraction"]) >= min_frac
        )

        def _source_passes(src):
            return any(
                r["source"] == src and float(r["crd1_fraction"]) >= min_frac
                for r in items
            )

        passing_sources = {s for s in ["af", "af2", "boltz"] if _source_passes(s)}

        keep = len(passing_sources) >= args.min_sources

        summary_rows.append(
            {
                "candidate_id": candidate_id,
                "models_total": models_total,
                "models_passing": models_passing,
                "sources_passing": len(passing_sources),
                "max_crd1_fraction": f"{float(best['crd1_fraction']):.4f}",
                "max_crd1_contacts": best["crd1_contacts"],
                "keep": keep,
                "best_model_path": best["model_path"],
            }
        )

    summary_path = out_dir / "prescreen_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys() if summary_rows else [])
        if summary_rows:
            writer.writeheader()
            writer.writerows(summary_rows)

    keep_ids = [row["candidate_id"] for row in summary_rows if row["keep"]]
    keep_path = out_dir / "prescreen_keep.txt"
    keep_path.write_text("\n".join(sorted(keep_ids)) + ("\n" if keep_ids else ""))

    # --- source_analysis.csv (3-source categorization) ---
    ALL_SOURCES = ["af", "af2", "boltz"]

    grouped_by_source = defaultdict(list)
    for row in rows:
        grouped_by_source[(row["candidate_id"], row["source"])].append(row)

    analysis_rows = []
    for candidate_id, items in sorted(grouped.items()):
        candidate_short = candidate_id.replace("_CD40", "")
        row_data = {"candidate": candidate_short}
        passing_sources = set()

        for src in ALL_SOURCES:
            src_items = grouped_by_source.get((candidate_id, src), [])
            src_total = len(src_items)
            src_passing = sum(
                1 for r in src_items if float(r["crd1_fraction"]) >= min_frac
            )
            src_max = max((float(r["crd1_fraction"]) for r in src_items), default=0.0)

            row_data[f"{src}_max_crd1"] = f"{src_max:.2f}"
            row_data[f"{src}_models_passing"] = src_passing
            row_data[f"{src}_total_models"] = src_total

            if src_passing > 0:
                passing_sources.add(src)

        if passing_sources == set(ALL_SOURCES):
            category = "ALL_THREE"
        elif len(passing_sources) == 0:
            category = "NEITHER"
        elif len(passing_sources) == 1:
            category = list(passing_sources)[0].upper() + "_ONLY"
        else:
            category = "+".join(sorted(s.upper() for s in passing_sources))

        row_data["category"] = category
        analysis_rows.append(row_data)

    analysis_fieldnames = ["candidate", "category"]
    for src in ALL_SOURCES:
        analysis_fieldnames.extend([
            f"{src}_max_crd1", f"{src}_models_passing", f"{src}_total_models",
        ])

    analysis_path = out_dir / "source_analysis.csv"
    with analysis_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=analysis_fieldnames)
        writer.writeheader()
        writer.writerows(analysis_rows)

    print(f"Wrote {summary_path}")
    print(f"Wrote {keep_path} ({len(keep_ids)} candidates)")
    print(f"Wrote {analysis_path}")


if __name__ == "__main__":
    main()
