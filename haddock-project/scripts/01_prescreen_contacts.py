#!/usr/bin/env python3
"""
Stage 1: CRD Contact Prescreen

Screens AF/Boltz structure predictions for CRD-targeting antibody candidates.
Computes CRD contact fraction and coverage per model.

Outputs:
    - prescreen_contacts.csv: Per-model contact metrics
"""
from __future__ import annotations


import argparse
import csv
from pathlib import Path

import yaml

from utils import (
    parse_range,
    find_candidate_id,
    model_source,
    source_chains,
    load_atoms,
    count_contacts,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prescreen AF/Boltz CIFs for CRD-contacting models."
    )
    parser.add_argument(
        "--config",
        default="",
        help="YAML config containing antigen chain (default: config/crd1_region.yaml).",
    )
    parser.add_argument(
        "--inputs",
        nargs="*",
        default=[
            "haddock-project/inputs/af3",
            "haddock-project/inputs/af2",
            "haddock-project/inputs/boltz",
        ],
        help="Input directories to scan for *_model*.cif files.",
    )
    parser.add_argument(
        "--candidate-dir",
        default="haddock-project/config/yaml",
        help="Directory of YAML files to restrict candidates (filenames used).",
    )
    parser.add_argument(
        "--candidate-list",
        default="",
        help="Optional text file with candidate IDs to include (one per line).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Contact cutoff distance in Angstroms.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory for prescreen results (default: outputs/prescreen/).",
    )
    return parser.parse_args()


def load_candidates_filter(candidate_dir: Path, candidate_list: Path | None):
    """Load set of candidate IDs to include."""
    if candidate_list and candidate_list.is_file():
        return {
            line.strip()
            for line in candidate_list.read_text().splitlines()
            if line.strip()
        }
    if candidate_dir and candidate_dir.exists():
        return {p.stem for p in candidate_dir.glob("*.y*ml")}
    return None


def main():
    args = parse_args()

    # Load both CRD regions
    crd1_config = Path(args.config) if args.config else Path("haddock-project/config/crd1_region.yaml")
    crd2_config = Path("haddock-project/config/crd2_region.yaml")

    if not crd1_config.exists():
        raise FileNotFoundError(f"Missing config: {crd1_config}")

    config = yaml.safe_load(crd1_config.read_text())
    antigen_chain = config["antigen"]["chain"]
    crd1_set = parse_range(config["antigen"]["crd_range"])
    crd2_set = parse_range(yaml.safe_load(crd2_config.read_text())["antigen"]["crd_range"]) if crd2_config.exists() else set()

    candidate_list = Path(args.candidate_list) if args.candidate_list else None
    candidates_filter = load_candidates_filter(Path(args.candidate_dir), candidate_list)

    # Find all model CIF files
    model_paths = []
    for input_dir in args.inputs:
        base = Path(input_dir)
        if not base.exists():
            continue
        model_paths.extend(base.rglob("*_model*.cif"))

    rows = []
    for path in sorted(model_paths):
        candidate_id = find_candidate_id(path)
        if not candidate_id:
            continue
        if candidates_filter and candidate_id not in candidates_filter:
            continue
        source = model_source(path)

        # Use source-dependent chain IDs (AF2 uses B/C/D or A/B/C depending on file type)
        src_ag_chain, src_ab_chains = source_chains(source, str(path))
        if source not in ("af2",):
            src_ag_chain = antigen_chain  # use config value for af/boltz

        # Load atoms (all antigen residues, not just CRD)
        antigen_atoms, antibody_atoms = load_atoms(
            path, src_ag_chain, src_ab_chains, crd_set=None
        )
        total_contacts, crd1_ct, crd1_frac = count_contacts(antigen_atoms, antibody_atoms, crd1_set, args.cutoff)
        _, crd2_ct, crd2_frac = count_contacts(antigen_atoms, antibody_atoms, crd2_set, args.cutoff)

        if total_contacts == 0:
            binding_region = "no_contact"
        elif crd1_ct > 0 and crd2_ct > 0:
            binding_region = "crd1+crd2"
        elif crd1_ct > 0:
            binding_region = "crd1"
        elif crd2_ct > 0:
            binding_region = "crd2"
        else:
            binding_region = "other"

        rows.append(
            {
                "candidate_id": candidate_id,
                "source": source,
                "model_path": str(path),
                "total_contact_residues": total_contacts,
                "crd1_contacts": crd1_ct,
                "crd2_contacts": crd2_ct,
                "crd1_fraction": f"{crd1_frac:.4f}",
                "crd2_fraction": f"{crd2_frac:.4f}",
                "crd1_coverage": f"{crd1_ct / len(crd1_set):.4f}" if crd1_set else "0.0000",
                "crd2_coverage": f"{crd2_ct / len(crd2_set):.4f}" if crd2_set else "0.0000",
                "binding_region": binding_region,
            }
        )

    out_dir = Path(args.output_dir) if args.output_dir else Path("haddock-project/outputs/prescreen")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write per-model contacts
    contacts_path = out_dir / "prescreen_contacts.csv"
    with contacts_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys() if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)

    print(f"Wrote {contacts_path} ({len(rows)} models)")


if __name__ == "__main__":
    main()
