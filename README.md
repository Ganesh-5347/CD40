# CD40 CRD1 Antibody Shortlisting Pipeline

Computational pipeline to identify and rank antibody candidates that bind the CRD1 domain (residues 1-42) of CD40, using structure predictions from AlphaFold2, AlphaFold3, and Boltz followed by physics-based refinement and independent validation.

## Pipeline Overview

```
01  Prescreen CRD1 contacts          (AF2/AF3/Boltz CIFs -> contact metrics)
     |
02  Filter candidates                (CRD1 fraction >= 0.5, require >= 2 sources)
     |
03  Convert best CIF -> PDB          (standardise chains: A/H/L)
     |
04  Prepare HADDOCK3 runs            (CDR-epitope AIRs + H-L anchoring + emref)
     |
    run_haddock_refine_batch.sh      (execute HADDOCK3 refinement)
     |
05  Collect HADDOCK3 scores          (parse capri_ss.tsv energy terms)
     |
06  Rosetta interface analysis       (shape complementarity, H-bonds, CDR contacts)
     |
07  Select best models               (merge HADDOCK + Rosetta, best per candidate)
     |
08  Final ranking                    (CRD category assignment, rank by HADDOCK score)
```

**Result:** 38 candidates -> 21 pass screening -> top 5 shortlisted for experimental follow-up.


## Prerequisites

- **Python 3.10+**
- **BioPython** -- CIF/PDB parsing (`pip install biopython`)
- **HADDOCK3** -- physics-based refinement ([haddock3 docs](https://www.bonvinlab.org/haddock3/))
- **PyRosetta** -- interface analysis (`conda install pyrosetta` in a dedicated env)
- Standard scientific Python: `numpy`, `csv`, `pathlib`

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/Ganesh-5347/CD40.git
cd CD40
```

### 2. Download structure predictions and HADDOCK3

```bash
bash haddock-project/scripts/download_predictions.sh
```

This downloads and sets up:
- `haddock-project/inputs/` -- 1,330 CIF files (82 MB compressed, 333 MB extracted)
  - `af3/` -- AlphaFold3: 38 candidates x 5 seeds = 190 CIFs
  - `boltz/` -- Boltz: 38 candidates x 5 samples = 190 CIFs
  - `af2/` -- AlphaFold2-multimer: 38 candidates x 5 architectures x 5 seeds = 950 CIFs
- `haddock3/` -- HADDOCK3 source (cloned from GitHub)

### 3. Install dependencies

```bash
# HADDOCK3
pip install -e haddock3/

# Python packages
pip install biopython numpy

# PyRosetta (for step 06, requires separate conda env)
conda create -n rosetta python=3.10
conda activate rosetta
conda install pyrosetta -c https://conda.rosetta.com
```

## Running the Pipeline

All scripts are run from the repository root. Each script has `--help` for full argument documentation.

```bash
# Stage 1: Prescreen all CIF models for CRD1 contacts
python haddock-project/scripts/01_prescreen_contacts.py

# Stage 2: Filter to candidates with CRD1 signal in >= 2 sources
python haddock-project/scripts/02_filter_candidates.py

# Stage 3: Convert best CIF per candidate/source to PDB
python haddock-project/scripts/03_build_complex_pdb.py

# Stage 4: Prepare HADDOCK3 refinement (split chains, generate AIRs + configs)
python haddock-project/scripts/04_prepare_haddock_runs.py

# Execute HADDOCK3 refinement (requires haddock3 on PATH)
bash haddock-project/scripts/run_haddock_refine_batch.sh

# Stage 5: Collect HADDOCK3 energy scores
python haddock-project/scripts/05_collect_haddock_scores.py

# Stage 6: Rosetta interface analysis (requires PyRosetta)
python haddock-project/scripts/06_rosetta_interface_analysis.py

# Stage 7: Merge HADDOCK + Rosetta scores, select best model per candidate
python haddock-project/scripts/07_select_best_models.py

# Stage 8: Assign CRD categories and produce final ranking
python haddock-project/scripts/08_final_ranking.py
```

## Key Results

| Metric | Value |
|--------|-------|
| Candidates screened | 38 |
| Models scanned (CIF) | 1330 (950 AF2 + 190 AF3 + 190 Boltz) |
| Passed CRD1 filter | 21 |
| HADDOCK3 refinements | 48 (one per candidate/source) |
| Rosetta validations | 48 |
| Shortlisted (top 5) | Ab-36.1, Ab-9.1, Ab-26.1, Ab-1.1, Ab-23.1 |

**Top candidate:** Ab-36.1 -- HADDOCK score -378.2, shape complementarity 0.80, 25 H-bonds, 6/6 CDR loops engaged, 97% CRD1-specific.

## Output Files

All outputs are written under `haddock-project/outputs/` (excluded from git):

| Path | Description |
|------|-------------|
| `outputs/prescreen/prescreen_contacts.csv` | Per-model contact metrics (1330 rows) |
| `outputs/prescreen/prescreen_summary.csv` | Per-candidate summary with keep/reject |
| `outputs/prescreen/prescreen_keep.txt` | Passing candidate IDs (21) |
| `outputs/complex_pdb/{af,af2,boltz}/*.pdb` | Best CIF converted to PDB (48 files) |
| `outputs/haddock_refine/*/run/` | HADDOCK3 refinement outputs |
| `outputs/reports/refine_scores.csv` | HADDOCK3 energy scores |
| `outputs/reports/rosetta_interface.csv` | Rosetta interface metrics |
| `outputs/reports/combined_scores.csv` | Merged best model per candidate (21 rows) |
| `outputs/reports/final_ranking.csv` | Final ranked list with CRD categories |
