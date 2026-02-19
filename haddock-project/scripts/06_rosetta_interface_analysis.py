#!/usr/bin/env python3
"""
Stage 6: PyRosetta Interface Analysis

Runs InterfaceAnalyzerMover on HADDOCK3-refined PDBs to extract
hydrogen bond counts, salt bridges, and interface energetics.
Additionally performs antibody-specific CDR-level contact analysis
using AntibodyInfo to identify which CDR loops drive binding.

Requires the `rosetta` conda environment with PyRosetta installed.

Outputs:
    - rosetta_interface.csv: Best model (HADDOCK rank 1) per candidate/source
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tempfile
from pathlib import Path

import pyrosetta
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.scoring.hbonds import HBondSet
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.antibody import AntibodyInfo, CDRNameEnum


# ---------------------------------------------------------------------------
# PyRosetta setup
# ---------------------------------------------------------------------------

INTERFACE = "HL_A"  # antibody chains H+L vs antigen chain A

# Charged residue sets for salt bridge detection
NEGATIVE_RESIDUES = {"ASP", "GLU"}
POSITIVE_RESIDUES = {"ARG", "LYS", "HIS"}

# Atoms carrying formal charge used for salt bridge distance check
NEGATIVE_ATOMS = {"OD1", "OD2", "OE1", "OE2"}  # ASP/GLU carboxylate oxygens
POSITIVE_ATOMS = {
    "NH1", "NH2", "NE",   # ARG
    "NZ",                  # LYS
    "ND1", "NE2",         # HIS
}

SALT_BRIDGE_CUTOFF = 4.0  # Angstroms
CONTACT_CUTOFF = 5.0      # Angstroms, heavy-atom distance for CDR contacts

CDR_NAMES = ["H1", "H2", "H3", "L1", "L2", "L3"]

# CRD domain residue ranges (sequence numbering on antigen chain A)
CRD1_RANGE = set(range(1, 43))    # residues 1-42
CRD2_RANGE = set(range(43, 92))   # residues 43-91


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="PyRosetta interface analysis of HADDOCK3-refined models.",
    )
    parser.add_argument(
        "--run-list",
        default="haddock-project/outputs/haddock_refine/run_dirs.txt",
        help="File containing refinement run directories.",
    )
    parser.add_argument(
        "--out-dir",
        default="haddock-project/outputs/reports",
        help="Output directory for results.",
    )
    return parser.parse_args()


def find_capri_ss(run_dir: Path) -> Path | None:
    """Find the last capri_ss.tsv in a HADDOCK3 run directory."""
    haddock_run = run_dir / "run"
    if not haddock_run.exists():
        return None
    capri_dirs = sorted(haddock_run.glob("*_caprieval"))
    if not capri_dirs:
        return None
    capri_ss = capri_dirs[-1] / "capri_ss.tsv"
    return capri_ss if capri_ss.exists() else None


def resolve_model_path(run_dir: Path, model_ref: str) -> Path | None:
    """Resolve a capri_ss.tsv model reference to an actual file on disk.

    The TSV stores paths like ``../2_mdref/mdref_4.pdb`` but the files on
    disk are gzipped (``mdref_4.pdb.gz``).  We check for both variants.
    """
    haddock_run = run_dir / "run"
    basename = Path(model_ref).name  # e.g. mdref_4.pdb

    # Search emref directories for the file (plain or gzipped)
    for emref_dir in sorted(haddock_run.glob("*_emref")):
        plain = emref_dir / basename
        if plain.exists():
            return plain
        gzipped = emref_dir / (basename + ".gz")
        if gzipped.exists():
            return gzipped

    # Search mdref directories for the file (plain or gzipped)
    for mdref_dir in sorted(haddock_run.glob("*_mdref")):
        plain = mdref_dir / basename
        if plain.exists():
            return plain
        gzipped = mdref_dir / (basename + ".gz")
        if gzipped.exists():
            return gzipped

    return None


def parse_ranked_models(capri_ss_path: Path, run_dir: Path) -> list[dict]:
    """Parse capri_ss.tsv and return rows sorted by caprieval_rank."""
    rows = []
    with capri_ss_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            model_ref = row.get("model", "")
            model_path = resolve_model_path(run_dir, model_ref)
            if model_path is None:
                continue
            rows.append({
                "model_ref": model_ref,
                "model_path": model_path,
                "caprieval_rank": int(row.get("caprieval_rank", 9999)),
            })
    rows.sort(key=lambda r: r["caprieval_rank"])
    return rows


def load_pose(pdb_path: Path) -> "pyrosetta.Pose":
    """Load a PDB (or .pdb.gz) into a PyRosetta Pose."""
    if pdb_path.suffix == ".gz":
        # Decompress to a temp file and load
        with gzip.open(pdb_path, "rt") as gz:
            pdb_text = gz.read()
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tmp:
            tmp.write(pdb_text)
            tmp_path = tmp.name
        pose = pyrosetta.pose_from_pdb(tmp_path)
        Path(tmp_path).unlink()
        return pose
    return pyrosetta.pose_from_pdb(str(pdb_path))


# ---------------------------------------------------------------------------
# Interface analysis
# ---------------------------------------------------------------------------

def get_interface_residues(pose, interface: str) -> tuple[set[int], set[int]]:
    """Return (side1_resnums, side2_resnums) based on chain membership.

    Splits all pose residues into two groups by PDB chain letter,
    according to the interface string (e.g. "HL_A" -> chains H,L vs A).
    """
    left_chains, right_chains = interface.split("_")

    left_resnums: set[int] = set()
    right_resnums: set[int] = set()

    pi = pose.pdb_info()
    for i in range(1, pose.total_residue() + 1):
        chain = pi.chain(i)
        if chain in left_chains:
            left_resnums.add(i)
        elif chain in right_chains:
            right_resnums.add(i)

    return left_resnums, right_resnums


def count_interface_hbonds(pose, left_res: set[int], right_res: set[int]) -> int:
    """Count hydrogen bonds that cross the interface.

    Scores the pose, extracts the HBondSet, and counts bonds where
    donor and acceptor are on opposite sides of the interface.
    """
    sfxn = get_score_function()
    sfxn(pose)  # score to populate hbond data

    hbset = HBondSet(pose, bb_only=False)
    count = 0
    for i in range(1, hbset.nhbonds() + 1):
        hb = hbset.hbond(i)
        don_res = hb.don_res()
        acc_res = hb.acc_res()
        # Count only cross-interface hbonds
        if (don_res in left_res and acc_res in right_res) or \
           (don_res in right_res and acc_res in left_res):
            count += 1
    return count


def count_salt_bridges(pose, left_res: set[int], right_res: set[int]) -> int:
    """Count salt bridges across the interface geometrically.

    A salt bridge is defined as any O–N pair between a negatively charged
    residue (ASP/GLU) on one side and a positively charged residue
    (ARG/LYS/HIS) on the other side with distance < 4.0 A.
    """
    salt_bridges = 0

    # Collect charged atoms from each side
    def _charged_atoms(resnums, pose):
        neg_atoms = []  # (xyz,)
        pos_atoms = []
        for rnum in resnums:
            res = pose.residue(rnum)
            resname = res.name3().strip()
            if resname in NEGATIVE_RESIDUES:
                for ai in range(1, res.natoms() + 1):
                    if res.atom_name(ai).strip() in NEGATIVE_ATOMS:
                        neg_atoms.append(res.xyz(ai))
            elif resname in POSITIVE_RESIDUES:
                for ai in range(1, res.natoms() + 1):
                    if res.atom_name(ai).strip() in POSITIVE_ATOMS:
                        pos_atoms.append(res.xyz(ai))
        return neg_atoms, pos_atoms

    left_neg, left_pos = _charged_atoms(left_res, pose)
    right_neg, right_pos = _charged_atoms(right_res, pose)

    # Cross-interface pairs: left_neg–right_pos and left_pos–right_neg
    for o_xyz in left_neg:
        for n_xyz in right_pos:
            if o_xyz.distance(n_xyz) < SALT_BRIDGE_CUTOFF:
                salt_bridges += 1

    for n_xyz in left_pos:
        for o_xyz in right_neg:
            if n_xyz.distance(o_xyz) < SALT_BRIDGE_CUTOFF:
                salt_bridges += 1

    return salt_bridges


def analyse_cdr_contacts(
    pose, ab_info, antigen_res: set[int], left_res: set[int],
) -> dict:
    """Per-CDR contact counts and hbond breakdown against the antigen.

    Parameters
    ----------
    pose : pyrosetta.Pose
    ab_info : AntibodyInfo
    antigen_res : set[int]
        Pose residue numbers on the antigen side (chain A).
    left_res : set[int]
        Pose residue numbers on the antibody side (chains H+L).
    """
    # Build lookup: pose resnum -> CDR name (or None for framework)
    resnum_to_cdr: dict[int, str | None] = {}
    cdr_resnums: dict[str, set[int]] = {}

    for cdr_name in CDR_NAMES:
        cdr_enum = getattr(CDRNameEnum, cdr_name)
        if not ab_info.has_CDR(cdr_enum):
            continue
        start = ab_info.get_CDR_start(cdr_enum, pose)
        end = ab_info.get_CDR_end(cdr_enum, pose)
        resnums = set(range(start, end + 1))
        cdr_resnums[cdr_name] = resnums
        for r in resnums:
            resnum_to_cdr[r] = cdr_name

    # Pre-collect antigen heavy-atom coordinates per residue
    ag_atoms: list[object] = []  # list of xyz vectors
    for rnum in antigen_res:
        res = pose.residue(rnum)
        for ai in range(1, res.natoms() + 1):
            if not res.atom_is_hydrogen(ai):
                ag_atoms.append(res.xyz(ai))

    # Per-CDR contact counts (antibody residue has any heavy atom < CONTACT_CUTOFF of any antigen heavy atom)
    contacts: dict[str, int] = {cdr: 0 for cdr in CDR_NAMES}
    framework_contacts = 0
    paratope_residues: set[int] = set()

    for rnum in left_res:
        res = pose.residue(rnum)
        in_contact = False
        for ai in range(1, res.natoms() + 1):
            if res.atom_is_hydrogen(ai):
                continue
            ab_xyz = res.xyz(ai)
            for ag_xyz in ag_atoms:
                if ab_xyz.distance(ag_xyz) < CONTACT_CUTOFF:
                    in_contact = True
                    break
            if in_contact:
                break
        if not in_contact:
            continue
        paratope_residues.add(rnum)
        cdr = resnum_to_cdr.get(rnum)
        if cdr is not None:
            contacts[cdr] += 1
        else:
            framework_contacts += 1

    # Per-CDR hbond counts (cross-interface, antibody donor/acceptor in a CDR)
    sfxn = get_score_function()
    sfxn(pose)
    hbset = HBondSet(pose, bb_only=False)

    hbonds: dict[str, int] = {cdr: 0 for cdr in CDR_NAMES}
    for i in range(1, hbset.nhbonds() + 1):
        hb = hbset.hbond(i)
        don_res = hb.don_res()
        acc_res = hb.acc_res()
        # Must be cross-interface
        if don_res in left_res and acc_res in antigen_res:
            ab_resnum = don_res
        elif acc_res in left_res and don_res in antigen_res:
            ab_resnum = acc_res
        else:
            continue
        cdr = resnum_to_cdr.get(ab_resnum)
        if cdr is not None:
            hbonds[cdr] += 1

    # H3 kink type
    h3_kink_type = "NA"
    h3_enum = CDRNameEnum.H3
    if ab_info.has_CDR(h3_enum):
        try:
            h3_kink_type = ab_info.get_H3_kink_type_name()
        except Exception:
            h3_kink_type = "NA"

    # Summary metrics
    n_cdrs_engaged = sum(1 for cdr in CDR_NAMES if contacts[cdr] > 0)
    dominant_cdr = max(CDR_NAMES, key=lambda c: contacts[c]) if n_cdrs_engaged > 0 else "none"
    paratope_nres = len(paratope_residues)

    result = {}
    for cdr in CDR_NAMES:
        key = cdr.lower()
        result[f"{key}_contacts"] = contacts[cdr]
        result[f"{key}_hbonds"] = hbonds[cdr]
    result["framework_contacts"] = framework_contacts
    result["n_cdrs_engaged"] = n_cdrs_engaged
    result["dominant_cdr"] = dominant_cdr
    result["h3_kink_type"] = h3_kink_type
    result["paratope_nres"] = paratope_nres

    return result


def analyse_crd_domains(pose, left_res: set[int], right_res: set[int]) -> dict:
    """Count antigen-side contacts in CRD1 vs CRD2 domains.

    For each antigen residue in right_res, checks if any antibody heavy atom
    is within CONTACT_CUTOFF. Then classifies the contacted antigen residues
    by PDB residue number into CRD1 (1-42) or CRD2 (43-91).
    """
    pi = pose.pdb_info()

    # Pre-collect antibody heavy-atom coordinates
    ab_atoms = []
    for rnum in left_res:
        res = pose.residue(rnum)
        for ai in range(1, res.natoms() + 1):
            if not res.atom_is_hydrogen(ai):
                ab_atoms.append(res.xyz(ai))

    contacted_ag_resids: set[int] = set()
    for rnum in right_res:
        res = pose.residue(rnum)
        in_contact = False
        for ai in range(1, res.natoms() + 1):
            if res.atom_is_hydrogen(ai):
                continue
            ag_xyz = res.xyz(ai)
            for ab_xyz in ab_atoms:
                if ag_xyz.distance(ab_xyz) < CONTACT_CUTOFF:
                    in_contact = True
                    break
            if in_contact:
                break
        if in_contact:
            contacted_ag_resids.add(pi.number(rnum))

    crd1_contacts = len(contacted_ag_resids & CRD1_RANGE)
    crd2_contacts = len(contacted_ag_resids & CRD2_RANGE)
    total_contacts = len(contacted_ag_resids)
    crd1_fraction = crd1_contacts / total_contacts if total_contacts else 0.0
    crd2_fraction = crd2_contacts / total_contacts if total_contacts else 0.0

    return {
        "crd1_contacts": crd1_contacts,
        "crd2_contacts": crd2_contacts,
        "crd1_fraction": round(crd1_fraction, 4),
        "crd2_fraction": round(crd2_fraction, 4),
    }


def analyse_interface(pose, interface: str = INTERFACE) -> dict:
    """Run InterfaceAnalyzerMover and return extracted metrics."""
    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface)
    iam.set_pack_separated(True)
    iam.set_compute_interface_sc(True)
    iam.set_compute_packstat(True)
    iam.apply(pose)

    data = iam.get_all_data()
    left_res, right_res = get_interface_residues(pose, interface)

    metrics = {
        "rosetta_dG": iam.get_interface_dG(),
        "rosetta_dSASA": iam.get_interface_delta_sasa(),
        "interface_hbonds": count_interface_hbonds(pose, left_res, right_res),
        "interface_saltbridges": count_salt_bridges(pose, left_res, right_res),
        "shape_complementarity": data.sc_value,
        "packstat": iam.get_interface_packstat(),
        "nres_interface": iam.get_num_interface_residues(),
    }

    # CDR-level contact analysis
    try:
        ab_info = AntibodyInfo(pose)
        cdr_metrics = analyse_cdr_contacts(pose, ab_info, right_res, left_res)
        metrics.update(cdr_metrics)
    except Exception as exc:
        print(f"    CDR analysis failed (will use NA): {exc}")
        for cdr in CDR_NAMES:
            key = cdr.lower()
            metrics[f"{key}_contacts"] = "NA"
            metrics[f"{key}_hbonds"] = "NA"
        metrics["framework_contacts"] = "NA"
        metrics["n_cdrs_engaged"] = "NA"
        metrics["dominant_cdr"] = "NA"
        metrics["h3_kink_type"] = "NA"
        metrics["paratope_nres"] = "NA"

    # CRD domain contact analysis (antigen-side)
    crd_domain = analyse_crd_domains(pose, left_res, right_res)
    metrics.update(crd_domain)

    return metrics


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

CSV_FIELDS = [
    "candidate_id",
    "source",
    "model_path",
    "rosetta_dG",
    "rosetta_dSASA",
    "interface_hbonds",
    "interface_saltbridges",
    "shape_complementarity",
    "packstat",
    "nres_interface",
    "h1_contacts",
    "h2_contacts",
    "h3_contacts",
    "l1_contacts",
    "l2_contacts",
    "l3_contacts",
    "h1_hbonds",
    "h2_hbonds",
    "h3_hbonds",
    "l1_hbonds",
    "l2_hbonds",
    "l3_hbonds",
    "framework_contacts",
    "n_cdrs_engaged",
    "dominant_cdr",
    "h3_kink_type",
    "paratope_nres",
    "crd1_contacts",
    "crd2_contacts",
    "crd1_fraction",
    "crd2_fraction",
]


def main():
    args = parse_args()

    run_list = Path(args.run_list)
    if not run_list.exists():
        raise FileNotFoundError(f"Run list not found: {run_list}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Initialise PyRosetta (suppress most output)
    pyrosetta.init(
        "-mute all -ignore_unrecognized_res true -ignore_zero_occupancy false"
    )

    rows: list[dict] = []

    lines = [l.strip() for l in run_list.read_text().splitlines() if l.strip()]

    for idx, line in enumerate(lines, 1):
        run_dir = Path(line)
        if not run_dir.exists():
            print(f"  [{idx}/{len(lines)}] SKIP (dir missing): {run_dir}")
            continue

        candidate = run_dir.parent.name
        source = run_dir.name

        capri_ss = find_capri_ss(run_dir)
        if capri_ss is None:
            print(f"  [{idx}/{len(lines)}] SKIP (no capri_ss): {run_dir}")
            continue

        ranked = parse_ranked_models(capri_ss, run_dir)
        if not ranked:
            print(f"  [{idx}/{len(lines)}] SKIP (no models): {run_dir}")
            continue

        # Analyse the best (rank-1) model
        best_entry = ranked[0]
        model_path = best_entry["model_path"]
        print(f"  [{idx}/{len(lines)}] {candidate}/{source}  ->  {model_path.name}")

        try:
            pose = load_pose(model_path)
            metrics = analyse_interface(pose)
        except Exception as exc:
            print(f"    ERROR: {exc}")
            continue

        rows.append({
            "candidate_id": candidate,
            "source": source,
            "model_path": str(model_path),
            **metrics,
        })

    # Write Rosetta CSV (all columns)
    out_path = out_dir / "rosetta_interface.csv"
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {out_path} ({len(rows)} models)")


if __name__ == "__main__":
    main()
