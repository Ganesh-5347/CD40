#!/usr/bin/env python3
"""
Stage 4: HADDOCK 3 Refinement Setup

Prepares HADDOCK 3 refinement runs from complex PDBs.
Splits complexes into separate chain PDBs, generates antibody-aware AIR
restraint files, and generates TOML config files.

Outputs:
    - haddock_refine/CANDIDATE/SOURCE/refine.cfg
    - haddock_refine/CANDIDATE/SOURCE/antigen.pdb
    - haddock_refine/CANDIDATE/SOURCE/antibody_H.pdb
    - haddock_refine/CANDIDATE/SOURCE/antibody_L.pdb
    - haddock_refine/CANDIDATE/SOURCE/ambig.tbl
    - haddock_refine/CANDIDATE/SOURCE/unambig.tbl
    - haddock_refine/run_dirs.txt
"""

import argparse
import math
from pathlib import Path

# Antibody-aware refinement defaults.
AIR_CUTOFF = 5.0
CRD1_EPITOPE_RESIDUES = set(range(1, 43))  # CRD1 on antigen chain A
HL_RESTRAINT_PAIRS = 2
FLEXREF_TOLERANCE = 5
MDREF_TOLERANCE = 5
MDREF_SAMPLING_FACTOR = 20


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare HADDOCK 3 refinement runs from complex PDBs."
    )
    parser.add_argument(
        "--input-root",
        default="haddock-project/outputs/complex_pdb",
        help="Root directory containing af/ af2/ and boltz/ complex PDBs.",
    )
    parser.add_argument(
        "--out-root",
        default="haddock-project/outputs/haddock_refine",
        help="Output root for refinement runs.",
    )
    parser.add_argument(
        "--ncores",
        type=int,
        default=20,
        help="Number of CPU cores for HADDOCK 3.",
    )
    parser.add_argument(
        "--antigen-chain",
        default="A",
        help="Antigen chain ID.",
    )
    parser.add_argument(
        "--antibody-heavy-chain",
        default="H",
        help="Antibody heavy chain ID.",
    )
    parser.add_argument(
        "--antibody-light-chain",
        default="L",
        help="Antibody light chain ID.",
    )
    return parser.parse_args()


def parse_pdb_heavy_atoms(pdb_path: Path) -> tuple[list[tuple], dict[int, tuple[float, float, float]]]:
    """Parse heavy atoms and CA coordinates from a PDB file."""
    atoms = []
    ca_by_resid = {}
    with pdb_path.open() as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_name = line[12:16].strip()
            chain = line[21:22].strip() or line[72:76].strip()[:1]

            resid_text = line[22:26].strip()
            if not resid_text:
                continue
            try:
                resid = int(resid_text)
            except ValueError:
                continue

            element = line[76:78].strip() or atom_name[:1]
            if element.upper() == "H":
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            atoms.append((x, y, z, chain, resid))
            if atom_name == "CA":
                ca_by_resid[resid] = (x, y, z)

    return atoms, ca_by_resid


def build_grid(atoms: list[tuple], cell_size: float) -> dict[tuple[int, int, int], list[tuple]]:
    """Build a spatial hash grid for fast contact detection."""
    grid = {}
    inv = 1.0 / cell_size
    for x, y, z, chain, resid in atoms:
        cell = (math.floor(x * inv), math.floor(y * inv), math.floor(z * inv))
        grid.setdefault(cell, []).append((x, y, z, chain, resid))
    return grid


def find_contact_pairs(
    antigen_atoms: list[tuple],
    antibody_atoms: list[tuple],
    cutoff: float,
) -> set[tuple[int, str, int]]:
    """
    Find antigen-antibody heavy-atom contacts.

    Returns:
        set of (antigen_resid, antibody_chain, antibody_resid).
    """
    if not antigen_atoms or not antibody_atoms:
        return set()

    grid = build_grid(antibody_atoms, cutoff)
    cutoff2 = cutoff * cutoff
    pairs = set()

    for ax, ay, az, _, ag_resid in antigen_atoms:
        cell = (math.floor(ax / cutoff), math.floor(ay / cutoff), math.floor(az / cutoff))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neigh = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    for bx, by, bz, ab_chain, ab_resid in grid.get(neigh, []):
                        dist2 = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
                        if dist2 <= cutoff2:
                            pairs.add((ag_resid, ab_chain, ab_resid))
    return pairs


def pick_hl_anchor_pairs(
    heavy_ca: dict[int, tuple[float, float, float]],
    light_ca: dict[int, tuple[float, float, float]],
    n_pairs: int,
) -> list[tuple[int, int, float]]:
    """Pick shortest H-L CA pairs as unambiguous restraints."""
    if not heavy_ca or not light_ca:
        return []

    candidates = []
    for h_resid, h_xyz in heavy_ca.items():
        for l_resid, l_xyz in light_ca.items():
            dx = h_xyz[0] - l_xyz[0]
            dy = h_xyz[1] - l_xyz[1]
            dz = h_xyz[2] - l_xyz[2]
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)
            candidates.append((dist, h_resid, l_resid))
    candidates.sort(key=lambda t: t[0])

    used_h = set()
    used_l = set()
    anchors = []
    for dist, h_resid, l_resid in candidates:
        if h_resid in used_h or l_resid in used_l:
            continue
        anchors.append((h_resid, l_resid, dist))
        used_h.add(h_resid)
        used_l.add(l_resid)
        if len(anchors) >= n_pairs:
            break
    return anchors


def write_ambig_tbl(
    ambig_path: Path,
    paratope_residues: list[tuple[str, int]],
    epitope_residues: list[int],
):
    """Write ambiguous AIR restraints (paratope -> epitope)."""
    lines = [
        "! HADDOCK AIR restraints for antibody-antigen refinement\n",
    ]
    for ab_chain, ab_resid in paratope_residues:
        lines.append(f"assign (resid {ab_resid} and segid {ab_chain})\n")
        lines.append("       (\n")
        for i, ag_resid in enumerate(epitope_residues):
            if i > 0:
                lines.append("     or\n")
            lines.append(f"        (resid {ag_resid} and segid A)\n")
        lines.append("       ) 2.0 2.0 0.0\n")
        lines.append("!\n")
    ambig_path.write_text("".join(lines))


def write_unambig_tbl(unambig_path: Path, anchors: list[tuple[int, int, float]]):
    """Write unambiguous H-L chain restraints."""
    lines = ["! Restraints to keep antibody H and L chains together\n"]
    for h_resid, l_resid, dist in anchors:
        lines.append(
            f"assign (resid {h_resid} and name CA and segid H) "
            f"(resid {l_resid} and name CA and segid L) "
            f"{dist:.3f} 0.00 0.00\n"
        )
    unambig_path.write_text("".join(lines))


def generate_antibody_restraints(
    antigen_path: Path,
    heavy_path: Path,
    light_path: Path,
    ambig_path: Path,
    unambig_path: Path,
):
    """Generate AIR files from the starting complex geometry."""
    antigen_atoms, _ = parse_pdb_heavy_atoms(antigen_path)
    heavy_atoms, heavy_ca = parse_pdb_heavy_atoms(heavy_path)
    light_atoms, light_ca = parse_pdb_heavy_atoms(light_path)
    antibody_atoms = heavy_atoms + light_atoms

    contact_pairs = find_contact_pairs(antigen_atoms, antibody_atoms, AIR_CUTOFF)
    if not contact_pairs:
        raise ValueError(f"No antigen-antibody contacts found in {antigen_path.parent}")

    # Prefer CRD1-focused epitope contacts; fallback to all contacts if empty.
    crd1_pairs = {p for p in contact_pairs if p[0] in CRD1_EPITOPE_RESIDUES}
    if crd1_pairs:
        contact_pairs = crd1_pairs

    epitope_residues = sorted({ag_resid for ag_resid, _, _ in contact_pairs})
    paratope_residues = sorted(
        {(ab_chain, ab_resid) for _, ab_chain, ab_resid in contact_pairs},
        key=lambda item: (item[0], item[1]),
    )
    if not epitope_residues or not paratope_residues:
        raise ValueError(f"Failed to derive AIR residue sets for {antigen_path.parent}")

    anchors = pick_hl_anchor_pairs(heavy_ca, light_ca, HL_RESTRAINT_PAIRS)
    if not anchors:
        raise ValueError(f"Failed to derive H-L anchors for {antigen_path.parent}")

    write_ambig_tbl(ambig_path, paratope_residues, epitope_residues)
    write_unambig_tbl(unambig_path, anchors)


def split_complex_pdb(
    pdb_path: Path,
    antigen_path: Path,
    heavy_path: Path,
    light_path: Path,
    antigen_chain: str,
    heavy_chain: str,
    light_chain: str,
):
    """Split complex PDB into separate chain files."""
    antigen_lines = []
    heavy_lines = []
    light_lines = []

    with pdb_path.open() as f:
        for line in f:
            record = line[:6]
            if record.startswith("ATOM") or record.startswith("HETATM") or record.startswith("TER"):
                chain_id = line[21:22]
                if chain_id == antigen_chain:
                    antigen_lines.append(line)
                elif chain_id == heavy_chain:
                    heavy_lines.append(line)
                elif chain_id == light_chain:
                    light_lines.append(line)

    if not antigen_lines or not heavy_lines or not light_lines:
        raise ValueError(f"Missing chains in {pdb_path}")

    if not antigen_lines[-1].startswith("END"):
        antigen_lines.append("END\n")
    if not heavy_lines[-1].startswith("END"):
        heavy_lines.append("END\n")
    if not light_lines[-1].startswith("END"):
        light_lines.append("END\n")

    antigen_path.write_text("".join(antigen_lines))
    heavy_path.write_text("".join(heavy_lines))
    light_path.write_text("".join(light_lines))


def write_haddock3_cfg(cfg_path: Path, ncores: int):
    """Write HADDOCK 3 refinement TOML config."""
    cfg_path.write_text(
        f"""\
# HADDOCK 3 antibody-antigen refinement config (AIR-guided, refinement-only)
run_dir = "run"
mode = "local"
ncores = {ncores}

molecules = [
    "antigen.pdb",
    "antibody_H.pdb",
    "antibody_L.pdb"
]

[topoaa]
autohis = true

[flexref]
tolerance = {FLEXREF_TOLERANCE}
ambig_fname = "ambig.tbl"
unambig_fname = "unambig.tbl"
randremoval = false

[mdref]
tolerance = {MDREF_TOLERANCE}
sampling_factor = {MDREF_SAMPLING_FACTOR}
ambig_fname = "ambig.tbl"
unambig_fname = "unambig.tbl"
randremoval = false

[emref]
tolerance = {FLEXREF_TOLERANCE}
ambig_fname = "ambig.tbl"
unambig_fname = "unambig.tbl"
randremoval = false

[caprieval]

[contactmap]

[prodigyprotein]
chains = ["H,L", "A"]
"""
    )


def main():
    args = parse_args()
    input_root = Path(args.input_root)
    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    run_dirs = []

    for source_dir in sorted(input_root.glob("*")):
        if not source_dir.is_dir():
            continue
        source = source_dir.name
        for pdb_path in sorted(source_dir.glob("*.pdb")):
            candidate = pdb_path.stem
            run_dir = out_root / candidate / source
            run_dir.mkdir(parents=True, exist_ok=True)

            antigen_pdb = run_dir / "antigen.pdb"
            heavy_pdb = run_dir / "antibody_H.pdb"
            light_pdb = run_dir / "antibody_L.pdb"
            ambig_tbl = run_dir / "ambig.tbl"
            unambig_tbl = run_dir / "unambig.tbl"
            split_complex_pdb(
                pdb_path,
                antigen_pdb,
                heavy_pdb,
                light_pdb,
                args.antigen_chain,
                args.antibody_heavy_chain,
                args.antibody_light_chain,
            )
            generate_antibody_restraints(
                antigen_path=antigen_pdb,
                heavy_path=heavy_pdb,
                light_path=light_pdb,
                ambig_path=ambig_tbl,
                unambig_path=unambig_tbl,
            )

            cfg_path = run_dir / "refine.cfg"
            write_haddock3_cfg(cfg_path, args.ncores)

            run_dirs.append(run_dir)

    run_list = out_root / "run_dirs.txt"
    run_list.write_text("\n".join(str(p) for p in run_dirs) + ("\n" if run_dirs else ""))

    print(f"Prepared {len(run_dirs)} refinement runs in {out_root}")
    print(f"Run list: {run_list}")


if __name__ == "__main__":
    main()
