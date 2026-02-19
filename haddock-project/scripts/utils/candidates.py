"""Candidate loading and filtering utilities."""
from __future__ import annotations


import csv
import re
from pathlib import Path


# Chain ID mapping per source: maps canonical names to source-specific CIF chain IDs.
# AF2 unrelaxed CIFs use B(heavy), C(light), D(antigen).
# AF2 ranked/relaxed CIFs use A(heavy), B(light), C(antigen).
CHAIN_MAP = {
    "af":          {"antigen": "A", "heavy": "H", "light": "L"},
    "af2":         {"antigen": "D", "heavy": "B", "light": "C"},
    "af2_relaxed": {"antigen": "C", "heavy": "A", "light": "B"},
    "boltz":       {"antigen": "A", "heavy": "H", "light": "L"},
}


def source_chains(source: str, model_path: str = "") -> tuple[str, tuple[str, str]]:
    """Return (antigen_chain, (heavy_chain, light_chain)) for a given source.

    For AF2, ranked/relaxed CIFs have different chain IDs than unrelaxed.
    """
    key = source
    if source == "af2" and model_path and "unrelaxed" not in model_path:
        key = "af2_relaxed"
    m = CHAIN_MAP.get(key, CHAIN_MAP["af"])
    return m["antigen"], (m["heavy"], m["light"])


def parse_range(range_text: str) -> set[int]:
    """
    Parse a residue range string (e.g., "1-42") into a set of integers.

    Args:
        range_text: Range string like "1-42" or single number.

    Returns:
        Set of residue numbers.
    """
    text = range_text.strip()
    if not text:
        raise ValueError("Empty range string")
    if "-" in text:
        start_s, end_s = text.split("-", 1)
        start = int(start_s.strip())
        end = int(end_s.strip())
    else:
        start = end = int(text)
    if end < start:
        raise ValueError(f"Invalid range: {text}")
    return set(range(start, end + 1))


def find_candidate_id(path: Path) -> str:
    """Extract candidate ID (e.g., Ab-1.1_CD40) from file path."""
    match = re.search(r"(Ab-[^/]+_CD40)", str(path))
    return match.group(1) if match else ""


def model_source(path: Path) -> str:
    """Determine model source (af/af2/boltz) from file path."""
    s = str(path)
    if "/af2/" in s:
        return "af2"
    if "/af3/" in s:
        return "af"
    if "/boltz/" in s:
        return "boltz"
    return "unknown"


def load_both_candidates(
    analysis_path: Path,
    append_cd40: bool = True,
    category: str = "BOTH",
) -> set[str]:
    """
    Load candidates matching a category from analysis CSV.

    Args:
        analysis_path: Path to source_analysis.csv (or legacy af_vs_boltz_analysis.csv).
        append_cd40: Whether to append "_CD40" suffix.
        category: Category to filter for (e.g. "BOTH", "ALL_THREE").

    Returns:
        Set of candidate IDs.
    """
    result = set()
    with analysis_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("category") == category:
                cand = row["candidate"]
                if append_cd40 and not cand.endswith("_CD40"):
                    cand = cand + "_CD40"
                result.add(cand)
    return result


def load_candidates_from_dir(
    keep_list: Path | None = None,
    candidate_dir: Path | None = None,
) -> list[str]:
    """
    Load candidate IDs from keep list or YAML directory.

    Args:
        keep_list: Path to text file with candidate IDs (one per line).
        candidate_dir: Directory containing YAML files (stems used as IDs).

    Returns:
        List of candidate IDs.
    """
    if keep_list and keep_list.exists():
        return [
            line.strip()
            for line in keep_list.read_text().splitlines()
            if line.strip()
        ]
    if candidate_dir and candidate_dir.exists():
        return [p.stem for p in candidate_dir.glob("*.y*ml")]
    raise FileNotFoundError("No keep list or candidate dir found")


def select_best_models(
    prescreen_path: Path,
    candidates: set[str],
    crd_prefix: str = "crd1",
) -> dict[tuple[str, str], dict]:
    """
    Select best model per candidate/source based on CRD metrics.

    Args:
        prescreen_path: Path to prescreen_contacts.csv.
        candidates: Set of candidate IDs to include.
        crd_prefix: Column prefix for CRD metrics (e.g. "crd1" or "crd2").

    Returns:
        Dict mapping (candidate, source) to {"score": tuple, "model_path": str, "row": dict}.
    """
    best = {}
    with prescreen_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cand = row["candidate_id"]
            if cand not in candidates:
                continue
            source = row["source"]
            if source not in ("af", "af2", "boltz"):
                continue
            score = (
                float(row[f"{crd_prefix}_fraction"]),
                float(row[f"{crd_prefix}_coverage"]),
                int(row[f"{crd_prefix}_contact_residues"]),
                int(row["total_contact_residues"]),
            )
            key = (cand, source)
            if key not in best or score > best[key]["score"]:
                best[key] = {"score": score, "model_path": row["model_path"], "row": row}
    return best
