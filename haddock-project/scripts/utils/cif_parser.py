"""mmCIF file parsing utilities."""
from __future__ import annotations


from pathlib import Path
from typing import Iterator


def iter_atom_site_records(cif_path: Path) -> Iterator[tuple[list[str], list[str]]]:
    """
    Iterate over atom_site records in an mmCIF file.

    Yields:
        (headers, parts): Column headers and values for each atom record.
    """
    headers = []
    in_atom_loop = False
    with cif_path.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("loop_"):
                headers = []
                in_atom_loop = False
                continue
            if line.startswith("_atom_site."):
                headers.append(line.split()[0])
                in_atom_loop = True
                continue
            if in_atom_loop:
                if line.startswith("_") or line.startswith("#") or line.startswith("loop_"):
                    break
                parts = line.split()
                if len(parts) != len(headers):
                    continue
                yield headers, parts


def load_atoms(
    cif_path: Path,
    antigen_chain: str,
    antibody_chains: tuple[str, ...],
    crd_set: set[int] | None = None,
) -> tuple[list, list]:
    """
    Load atoms from CIF file for antigen and antibody chains.

    Args:
        cif_path: Path to mmCIF file.
        antigen_chain: Chain ID for antigen (e.g., "A").
        antibody_chains: Chain IDs for antibody (e.g., ("H", "L")).
        crd_set: Optional set of residue IDs to filter antigen atoms.
                 If None, all antigen residues are included.

    Returns:
        (antigen_atoms, antibody_atoms):
            antigen_atoms: List of (x, y, z, resid) tuples.
            antibody_atoms: List of (x, y, z, chain, resid) tuples.
    """
    antigen_atoms = []
    antibody_atoms = []
    idx_map = None

    for headers, parts in iter_atom_site_records(cif_path):
        if idx_map is None:
            idx_map = {h: i for i, h in enumerate(headers)}
            required = [
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.type_symbol",
            ]
            for req in required:
                if req not in idx_map:
                    raise ValueError(f"Missing {req} in {cif_path}")

        def get_value(key):
            idx = idx_map.get(key)
            return parts[idx] if idx is not None else None

        element = get_value("_atom_site.type_symbol")
        if element == "H":
            continue

        group = get_value("_atom_site.group_PDB")
        if group not in ("ATOM", "HETATM", None):
            continue

        chain = get_value("_atom_site.auth_asym_id")
        if chain in (None, ".", "?"):
            chain = get_value("_atom_site.label_asym_id")

        if chain not in (antigen_chain, *antibody_chains):
            continue

        seq = get_value("_atom_site.auth_seq_id")
        if seq in (None, ".", "?"):
            seq = get_value("_atom_site.label_seq_id")
        try:
            resid = int(seq)
        except (TypeError, ValueError):
            continue

        x = float(get_value("_atom_site.Cartn_x"))
        y = float(get_value("_atom_site.Cartn_y"))
        z = float(get_value("_atom_site.Cartn_z"))

        if chain == antigen_chain:
            if crd_set is None or resid in crd_set:
                antigen_atoms.append((x, y, z, resid))
        else:
            antibody_atoms.append((x, y, z, chain, resid))

    return antigen_atoms, antibody_atoms


def format_atom_name(atom_name: str, element: str) -> str:
    """Format atom name for PDB output."""
    if len(atom_name) >= 4:
        return atom_name[:4]
    if atom_name and atom_name[0].isdigit():
        return f"{atom_name:<4}"
    if len(element) == 1:
        return f" {atom_name:<3}"
    return f"{atom_name:>4}"


def cif_to_pdb_lines(
    cif_path: Path,
    chains_keep: set[str],
    chain_remap: dict[str, str] | None = None,
) -> list[str]:
    """
    Convert mmCIF file to PDB format lines.

    Args:
        cif_path: Path to mmCIF file.
        chains_keep: Set of chain IDs to include (matched against original CIF IDs).
        chain_remap: Optional mapping to rename chains in the output PDB,
                     e.g. ``{"A": "H", "B": "L", "C": "A"}``.

    Returns:
        List of PDB format lines.
    """
    idx_map = None
    serial = 1
    last_chain = None
    last_out_chain = None
    last_res_info = None
    lines = []

    for headers, parts in iter_atom_site_records(cif_path):
        if idx_map is None:
            idx_map = {h: i for i, h in enumerate(headers)}

        def get_value(key):
            idx = idx_map.get(key)
            return parts[idx] if idx is not None else None

        group = get_value("_atom_site.group_PDB")
        if group not in ("ATOM", "HETATM"):
            continue

        element = get_value("_atom_site.type_symbol") or ""
        if element == "H":
            continue

        chain = get_value("_atom_site.auth_asym_id")
        if chain in (None, ".", "?"):
            chain = get_value("_atom_site.label_asym_id")
        if chain not in chains_keep:
            continue
        out_chain = chain_remap.get(chain, chain) if chain_remap else chain

        atom_name = get_value("_atom_site.label_atom_id") or ""
        alt_id = get_value("_atom_site.label_alt_id")
        if alt_id in (None, ".", "?"):
            alt_id = " "

        res_name = get_value("_atom_site.label_comp_id") or "UNK"
        ins_code = get_value("_atom_site.pdbx_PDB_ins_code")
        if ins_code in (None, ".", "?"):
            ins_code = " "

        seq = get_value("_atom_site.auth_seq_id")
        if seq in (None, ".", "?"):
            seq = get_value("_atom_site.label_seq_id")
        try:
            res_seq = int(seq)
        except (TypeError, ValueError):
            continue

        x = float(get_value("_atom_site.Cartn_x"))
        y = float(get_value("_atom_site.Cartn_y"))
        z = float(get_value("_atom_site.Cartn_z"))
        occ = get_value("_atom_site.occupancy")
        try:
            occ_f = float(occ)
        except (TypeError, ValueError):
            occ_f = 1.00
        b_iso = get_value("_atom_site.B_iso_or_equiv")
        try:
            b_iso_f = float(b_iso)
        except (TypeError, ValueError):
            b_iso_f = 0.00

        if last_chain is not None and chain != last_chain and last_res_info:
            prev_res_name, prev_res_seq, prev_ins_code = last_res_info
            lines.append(
                f"TER   {serial:>5d}      {prev_res_name:>3s} {last_out_chain:1s}{prev_res_seq:>4d}{prev_ins_code:1s}\n"
            )
            serial += 1

        atom_field = format_atom_name(atom_name, element)
        line = (
            f"{group:<6s}{serial:>5d} {atom_field}{alt_id:1s}"
            f"{res_name:>3s} {out_chain:1s}{res_seq:>4d}{ins_code:1s}   "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{occ_f:>6.2f}{b_iso_f:>6.2f}          {element:>2s}\n"
        )
        lines.append(line)
        serial += 1
        last_chain = chain
        last_out_chain = out_chain
        last_res_info = (res_name, res_seq, ins_code)

    if lines:
        lines.append("END\n")
    return lines
