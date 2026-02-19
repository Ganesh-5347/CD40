#!/usr/bin/env python3
"""
Stage 3: CIF to PDB Conversion

Converts AF3/AF2/Boltz mmCIF structures to PDB format for HADDOCK refinement.

Outputs:
    - complex_pdb/af/*.pdb: AlphaFold3 complex PDBs
    - complex_pdb/af2/*.pdb: AlphaFold2 complex PDBs
    - complex_pdb/boltz/*.pdb: Boltz complex PDBs
"""

import argparse
from pathlib import Path

from utils import load_candidates_from_dir, cif_to_pdb_lines


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build complex PDBs from AF3/AF2/Boltz CIFs for selected candidates."
    )
    parser.add_argument(
        "--keep-list",
        default="haddock-project/outputs/prescreen/prescreen_keep.txt",
        help="Text file with candidate IDs to include (one per line).",
    )
    parser.add_argument(
        "--candidate-dir",
        default="",
        help="Optional directory of YAML files (filenames used as candidates).",
    )
    parser.add_argument(
        "--af-root",
        default="haddock-project/inputs/af3",
        help="Root directory for AF3 outputs.",
    )
    parser.add_argument(
        "--af2-root",
        default="haddock-project/inputs/af2",
        help="Root directory for AF2 outputs.",
    )
    parser.add_argument(
        "--boltz-root",
        default="haddock-project/inputs/boltz",
        help="Root directory for Boltz outputs.",
    )
    parser.add_argument(
        "--contacts-csv",
        default="haddock-project/outputs/prescreen/prescreen_contacts.csv",
        help="Prescreen contacts CSV to select best CRD1 model per source.",
    )
    parser.add_argument(
        "--out-root",
        default="haddock-project/outputs/complex_pdb",
        help="Output directory for complex PDBs.",
    )
    parser.add_argument(
        "--antigen-chain",
        default="A",
        help="Antigen chain ID (AF3/Boltz convention).",
    )
    parser.add_argument(
        "--antibody-chains",
        default="H,L",
        help="Comma-separated antibody chain IDs (AF3/Boltz convention).",
    )
    return parser.parse_args()


def load_best_models(contacts_csv: Path, candidates: list[str], min_frac: float = 0.5) -> dict:
    """Select best CRD1 model per candidate/source from prescreen CSV.

    Only includes (candidate, source) pairs where the best model meets *min_frac*.
    """
    import csv
    best = {}  # (candidate, source) -> {"frac": float, "path": str}
    with contacts_csv.open(newline="") as f:
        for row in csv.DictReader(f):
            cand = row["candidate_id"]
            if cand not in candidates:
                continue
            src = row["source"]
            frac = float(row["crd1_fraction"])
            key = (cand, src)
            if key not in best or frac > best[key]["frac"]:
                best[key] = {"frac": frac, "path": row["model_path"]}
    return {k: v for k, v in best.items() if v["frac"] >= min_frac}


def main():
    args = parse_args()
    candidates = load_candidates_from_dir(
        Path(args.keep_list) if args.keep_list else None,
        Path(args.candidate_dir) if args.candidate_dir else None,
    )

    contacts_csv = Path(args.contacts_csv)
    best_models = load_best_models(contacts_csv, set(candidates))

    out_root = Path(args.out_root)

    # Standard chains for AF3/Boltz CIFs
    chains_keep = {args.antigen_chain}
    chains_keep.update({c.strip() for c in args.antibody_chains.split(",") if c.strip()})

    # AF2 CIFs use different chain IDs depending on the file type.
    # ranked/relaxed CIFs use A=heavy, B=light, C=antigen.
    # unrelaxed model CIFs use B=heavy, C=light, D=antigen.
    af2_relaxed_chains = {"A", "B", "C"}
    af2_relaxed_remap = {"A": "H", "B": "L", "C": "A"}
    af2_unrelaxed_chains = {"B", "C", "D"}
    af2_unrelaxed_remap = {"B": "H", "C": "L", "D": "A"}

    counts = {"af": 0, "af2": 0, "boltz": 0}

    for candidate in candidates:
        for source in ["af", "af2", "boltz"]:
            key = (candidate, source)
            if key not in best_models:
                continue

            cif_path = Path(best_models[key]["path"])
            if not cif_path.exists():
                print(f"  Warning: {cif_path} not found, skipping")
                continue

            source_out = out_root / source
            source_out.mkdir(parents=True, exist_ok=True)

            if source == "af2":
                # Determine chain mapping based on filename
                if "unrelaxed_" in cif_path.name:
                    lines = cif_to_pdb_lines(cif_path, af2_unrelaxed_chains, chain_remap=af2_unrelaxed_remap)
                else:
                    lines = cif_to_pdb_lines(cif_path, af2_relaxed_chains, chain_remap=af2_relaxed_remap)
            else:
                lines = cif_to_pdb_lines(cif_path, chains_keep)

            if not lines:
                print(f"  Warning: no atoms for {source} {candidate}, skipping")
                continue

            (source_out / f"{candidate}.pdb").write_text("".join(lines))
            counts[source] += 1

    for source in ["af", "af2", "boltz"]:
        print(f"Wrote {counts[source]} {source} PDBs to {out_root / source}")


if __name__ == "__main__":
    main()
