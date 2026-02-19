"""Shared utilities for CD40 antibody screening pipeline."""

from .cif_parser import iter_atom_site_records, load_atoms, cif_to_pdb_lines
from .contacts import build_grid, contact_pairs, count_contacts
from .candidates import (
    parse_range,
    load_both_candidates,
    load_candidates_from_dir,
    select_best_models,
    find_candidate_id,
    model_source,
    source_chains,
    CHAIN_MAP,
)
from .confidence import (
    af_residue_plddt,
    af2_residue_plddt,
    af2_pae_index_map,
    boltz_residue_plddt,
    boltz_chain_lengths,
    af_pae_index_map,
    boltz_pae_index_map,
)

__all__ = [
    # CIF parsing
    "iter_atom_site_records",
    "load_atoms",
    "cif_to_pdb_lines",
    # Contacts
    "build_grid",
    "contact_pairs",
    "count_contacts",
    # Candidates
    "parse_range",
    "load_both_candidates",
    "load_candidates_from_dir",
    "select_best_models",
    "find_candidate_id",
    "model_source",
    "source_chains",
    "CHAIN_MAP",
    # Confidence
    "af_residue_plddt",
    "af2_residue_plddt",
    "af2_pae_index_map",
    "boltz_residue_plddt",
    "boltz_chain_lengths",
    "af_pae_index_map",
    "boltz_pae_index_map",
]
