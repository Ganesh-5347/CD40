"""Confidence metric utilities (pLDDT, PAE) for AF2, AF3 and Boltz models."""
from __future__ import annotations


import json
from pathlib import Path

import numpy as np

from .cif_parser import iter_atom_site_records


def boltz_chain_lengths(records_path: Path) -> dict[str, int]:
    """
    Load chain lengths from Boltz records JSON.

    Args:
        records_path: Path to Boltz records JSON file.

    Returns:
        Dict mapping chain ID to number of residues.
    """
    with records_path.open() as handle:
        data = json.load(handle)
    return {ch["chain_name"]: ch["num_residues"] for ch in data["chains"]}


def af_residue_plddt(conf_path: Path, cif_path: Path) -> dict[tuple[str, int], float]:
    """
    Compute per-residue pLDDT from AlphaFold3 confidences.

    Args:
        conf_path: Path to *_confidences.json file.
        cif_path: Path to corresponding CIF file.

    Returns:
        Dict mapping (chain, resid) to average pLDDT.
    """
    with conf_path.open() as handle:
        data = json.load(handle)
    atom_plddt = data["atom_plddts"]

    # Collect atoms from CIF to match with pLDDT array
    all_atoms = []
    heavy_atoms = []
    idx_map = None
    for headers, parts in iter_atom_site_records(cif_path):
        if idx_map is None:
            idx_map = {h: i for i, h in enumerate(headers)}

        def get_value(key):
            idx = idx_map.get(key)
            return parts[idx] if idx is not None else None

        element = get_value("_atom_site.type_symbol")
        chain = get_value("_atom_site.auth_asym_id")
        if chain in (None, ".", "?"):
            chain = get_value("_atom_site.label_asym_id")
        seq = get_value("_atom_site.auth_seq_id")
        if seq in (None, ".", "?"):
            seq = get_value("_atom_site.label_seq_id")
        try:
            resid = int(seq)
        except (TypeError, ValueError):
            continue
        all_atoms.append((chain, resid, element))
        if element != "H":
            heavy_atoms.append((chain, resid, element))

    # Match atom records to pLDDT array length
    if len(heavy_atoms) == len(atom_plddt):
        atom_records = heavy_atoms
    elif len(all_atoms) == len(atom_plddt):
        atom_records = all_atoms
    else:
        raise ValueError(
            f"Atom count mismatch for {cif_path}: CIF atoms={len(all_atoms)}, "
            f"heavy atoms={len(heavy_atoms)}, pLDDT atoms={len(atom_plddt)}"
        )

    # Aggregate per residue
    residue_values = {}
    residue_counts = {}
    for (chain, resid, _), plddt in zip(atom_records, atom_plddt):
        if chain is None:
            continue
        key = (chain, resid)
        residue_values[key] = residue_values.get(key, 0.0) + float(plddt)
        residue_counts[key] = residue_counts.get(key, 0) + 1

    return {key: total / residue_counts[key] for key, total in residue_values.items()}


def boltz_residue_plddt(plddt_path: Path, records_path: Path) -> dict[tuple[str, int], float]:
    """
    Compute per-residue pLDDT from Boltz pLDDT array.

    Args:
        plddt_path: Path to plddt_*.npz file.
        records_path: Path to Boltz records JSON file.

    Returns:
        Dict mapping (chain, resid) to pLDDT value.
    """
    lengths = boltz_chain_lengths(records_path)
    h_len = lengths["H"]
    l_len = lengths["L"]
    a_len = lengths["A"]

    plddt = np.load(plddt_path)["plddt"].astype(float)
    # Boltz pLDDT may be 0-1 scale; convert to 0-100 if needed
    if np.nanmax(plddt) <= 1.5:
        plddt = plddt * 100.0

    residue_plddt = {}
    # Order is H, L, A in Boltz output
    for resid in range(1, h_len + 1):
        residue_plddt[("H", resid)] = float(plddt[resid - 1])
    for resid in range(1, l_len + 1):
        residue_plddt[("L", resid)] = float(plddt[h_len + resid - 1])
    for resid in range(1, a_len + 1):
        residue_plddt[("A", resid)] = float(plddt[h_len + l_len + resid - 1])
    return residue_plddt


def af_pae_index_map(conf_path: Path) -> tuple[np.ndarray, dict[tuple[str, int], int]]:
    """
    Load PAE matrix and build index map from AlphaFold3 confidences.

    Args:
        conf_path: Path to *_confidences.json file.

    Returns:
        (pae_matrix, index_map): PAE array and mapping of (chain, resid) to matrix index.
    """
    with conf_path.open() as handle:
        data = json.load(handle)
    pae = np.array(data["pae"], dtype=float)
    mapping = {}
    for i, (c, r) in enumerate(zip(data["token_chain_ids"], data["token_res_ids"])):
        mapping[(c, int(r))] = i
    return pae, mapping


def boltz_pae_index_map(records_path: Path) -> tuple[callable, dict[str, int]]:
    """
    Build PAE index mapper function for Boltz models.

    Args:
        records_path: Path to Boltz records JSON file.

    Returns:
        (mapper_func, lengths): Function mapping (chain, resid) to index,
                                and dict of chain lengths.
    """
    lengths = boltz_chain_lengths(records_path)
    offsets = {
        "H": 0,
        "L": lengths["H"],
        "A": lengths["H"] + lengths["L"],
    }

    def mapper(chain: str, resid: int) -> int | None:
        offset = offsets.get(chain)
        if offset is None:
            return None
        return offset + (resid - 1)

    return mapper, lengths


# ---------------------------------------------------------------------------
# AlphaFold2 helpers
# ---------------------------------------------------------------------------

def _af2_chain_residue_map(cif_path: Path) -> list[tuple[str, int]]:
    """Build ordered (chain, resid) list from an AF2 CIF file.

    AF2 confidence/PAE arrays use continuous indexing across chains in the
    order they appear in the CIF atom records.  This function returns that
    order so the arrays can be mapped back to (chain, resid) tuples.
    """
    seen = {}          # (chain, resid) -> first-seen index (for ordering)
    ordered = []
    for headers, parts in iter_atom_site_records(cif_path):
        idx_map = {h: i for i, h in enumerate(headers)}

        chain = parts[idx_map.get("_atom_site.auth_asym_id", -1)]
        if chain in (None, ".", "?"):
            chain = parts[idx_map.get("_atom_site.label_asym_id", -1)]
        seq = parts[idx_map.get("_atom_site.auth_seq_id", -1)]
        if seq in (None, ".", "?"):
            seq = parts[idx_map.get("_atom_site.label_seq_id", -1)]
        try:
            resid = int(seq)
        except (TypeError, ValueError):
            continue

        key = (chain, resid)
        if key not in seen:
            seen[key] = len(ordered)
            ordered.append(key)
    return ordered


def af2_residue_plddt(
    conf_path: Path,
    cif_path: Path,
) -> dict[tuple[str, int], float]:
    """Compute per-residue pLDDT from AlphaFold2 confidence JSON.

    AF2 ``confidence_model_*.json`` contains per-residue (not per-atom)
    scores in the ``confidenceScore`` array, indexed by continuous residue
    numbering across all chains.

    Args:
        conf_path: Path to confidence_model_*_multimer_v3_pred_*.json.
        cif_path:  Path to the corresponding CIF file (for chain mapping).

    Returns:
        Dict mapping (chain, resid) to pLDDT value (0-100 scale).
    """
    with conf_path.open() as handle:
        data = json.load(handle)
    scores = data["confidenceScore"]

    residue_order = _af2_chain_residue_map(cif_path)

    if len(residue_order) != len(scores):
        raise ValueError(
            f"AF2 residue count mismatch: CIF has {len(residue_order)} residues "
            f"but confidence JSON has {len(scores)} entries in {conf_path}"
        )

    return {key: float(scores[i]) for i, key in enumerate(residue_order)}


def af2_pae_index_map(
    pae_path: Path,
    cif_path: Path,
) -> tuple[np.ndarray, dict[tuple[str, int], int]]:
    """Load PAE matrix and build index map from AlphaFold2 PAE JSON.

    AF2 ``pae_model_*.json`` is a list containing one dict with keys
    ``predicted_aligned_error`` (NxN matrix) and ``max_predicted_aligned_error``.

    Args:
        pae_path: Path to pae_model_*_multimer_v3_pred_*.json.
        cif_path: Path to the corresponding CIF file (for chain mapping).

    Returns:
        (pae_matrix, index_map): PAE array and mapping of (chain, resid) to
        matrix index.
    """
    with pae_path.open() as handle:
        data = json.load(handle)
    # AF2 wraps the dict in a list
    if isinstance(data, list):
        data = data[0]
    pae = np.array(data["predicted_aligned_error"], dtype=float)

    residue_order = _af2_chain_residue_map(cif_path)
    mapping = {key: i for i, key in enumerate(residue_order)}
    return pae, mapping
