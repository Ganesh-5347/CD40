"""Contact detection utilities using spatial grid search."""
from __future__ import annotations


import math
from pathlib import Path

from .cif_parser import load_atoms


def build_grid(atoms: list, cell_size: float, include_chain: bool = False) -> dict:
    """
    Build spatial grid for efficient contact search.

    Args:
        atoms: List of atom tuples.
               If include_chain=False: [(x, y, z), ...]
               If include_chain=True: [(x, y, z, chain, resid), ...]
        cell_size: Grid cell size in Angstroms.
        include_chain: Whether atoms include chain/resid info.

    Returns:
        Grid dictionary mapping cell coordinates to atom lists.
    """
    grid = {}
    inv = 1.0 / cell_size
    for atom in atoms:
        x, y, z = atom[0], atom[1], atom[2]
        cell = (math.floor(x * inv), math.floor(y * inv), math.floor(z * inv))
        grid.setdefault(cell, []).append(atom)
    return grid


def contact_pairs(
    cif_path: Path,
    antigen_chain: str,
    antibody_chains: tuple[str, ...],
    crd_set: set[int],
    cutoff: float,
) -> set[tuple[int, str, int]]:
    """
    Find contact residue pairs between CRD and antibody.

    Args:
        cif_path: Path to mmCIF file.
        antigen_chain: Antigen chain ID.
        antibody_chains: Antibody chain IDs.
        crd_set: Set of CRD residue IDs.
        cutoff: Contact distance cutoff in Angstroms.

    Returns:
        Set of (antigen_resid, ab_chain, ab_resid) contact pairs.
    """
    antigen_atoms, antibody_atoms = load_atoms(
        cif_path, antigen_chain, antibody_chains, crd_set
    )
    if not antigen_atoms or not antibody_atoms:
        return set()

    grid = build_grid(antibody_atoms, cutoff, include_chain=True)
    cutoff2 = cutoff * cutoff
    pairs = set()

    for ax, ay, az, resid_a in antigen_atoms:
        cell = (math.floor(ax / cutoff), math.floor(ay / cutoff), math.floor(az / cutoff))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    for bx, by, bz, chain_b, resid_b in grid.get(neighbor, []):
                        dist2 = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
                        if dist2 <= cutoff2:
                            pairs.add((resid_a, chain_b, resid_b))
    return pairs


def count_contacts(
    antigen_atoms: list,
    antibody_atoms: list,
    crd_set: set[int],
    cutoff: float,
) -> tuple[int, int, float]:
    """
    Count contact residues between antigen and antibody.

    Args:
        antigen_atoms: List of (x, y, z, resid) tuples.
        antibody_atoms: List of (x, y, z) or (x, y, z, chain, resid) tuples.
        crd_set: Set of CRD residue IDs.
        cutoff: Contact distance cutoff in Angstroms.

    Returns:
        (total_contact_residues, crd_contact_residues, crd_fraction)
    """
    if not antigen_atoms or not antibody_atoms:
        return 0, 0, 0.0

    # Build grid from antibody atoms (just xyz)
    ab_xyz = [(a[0], a[1], a[2]) for a in antibody_atoms]
    grid = {}
    inv = 1.0 / cutoff
    for x, y, z in ab_xyz:
        cell = (math.floor(x * inv), math.floor(y * inv), math.floor(z * inv))
        grid.setdefault(cell, []).append((x, y, z))

    cutoff2 = cutoff * cutoff

    contact_residues = set()
    crd_contact_residues = set()

    for ax, ay, az, resid in antigen_atoms:
        cell = (math.floor(ax / cutoff), math.floor(ay / cutoff), math.floor(az / cutoff))
        found = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    for bx, by, bz in grid.get(neighbor, []):
                        dist2 = (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
                        if dist2 <= cutoff2:
                            contact_residues.add(resid)
                            if resid in crd_set:
                                crd_contact_residues.add(resid)
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break

    total_contacts = len(contact_residues)
    crd_contacts = len(crd_contact_residues)
    crd_fraction = (crd_contacts / total_contacts) if total_contacts else 0.0
    return total_contacts, crd_contacts, crd_fraction
