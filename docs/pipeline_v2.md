# Pipeline Documentation

CD40-CRD1 antibody shortlisting pipeline. Takes structure predictions from AF2, AF3, and Boltz, screens for CRD1-binding candidates, refines in HADDOCK3 with antibody-aware restraints, validates with Rosetta, and ranks by binding energy.

## Overview

```
01 Prescreen CRD1 contacts (AF2 + AF3 + Boltz CIFs -> contact metrics)
    |
02 Filter candidates (CRD1 fraction >= 0.5, require >= 2/3 sources)
    |
03 Convert best CIF -> PDB (standardise chains A/H/L for HADDOCK3)
    |
04 Prepare HADDOCK3 runs (CDR-epitope AIRs + H-L anchoring restraints)
    |
run_haddock_refine_batch.sh (HADDOCK3: topoaa -> flexref -> emref -> caprieval)
    |
05 Collect HADDOCK3 scores (parse capri_ss.tsv energy terms)
    |
06 Rosetta interface analysis (shape complementarity, H-bonds, CDR contacts)
    |
07 Select best models (merge HADDOCK + Rosetta, one per candidate)
    |
08 Final ranking (CRD category assignment, rank by HADDOCK score)
```

**Key parameters:**
- CRD1 region: Residues 1-42 of chain A (CD40 antigen)
- Contact cutoff: 5.0 Angstroms
- Refinement: HADDOCK3 with CDR-epitope AIRs + H-L anchoring restraints
- Validation: Rosetta InterfaceAnalyzer + AntibodyInfo (independent of HADDOCK)
- Sources: AlphaFold2 (ColabFold), AlphaFold3, Boltz

**Results:** 38 candidates screened, 1330 models scanned -> 21 passed CRD1 filter -> 48 HADDOCK3 refinements -> top 5 shortlisted (Ab-36.1, Ab-9.1, Ab-26.1, Ab-1.1, Ab-23.1).

---

## Step-by-Step Details

### Step 01: Contact Prescreen

**Script:** `01_prescreen_contacts.py`

Scans all CIF model files across three input directories and computes antibody-antigen contacts at 5A heavy-atom cutoff. For each model, calculates CRD1 and CRD2 contact counts, fractions, and coverage. Labels each model with a binding region: `crd1`, `crd2`, `crd1+crd2`, `other`, or `no_contact`.

**Input:** CIF files from `inputs/af2/`, `inputs/af3/`, `inputs/boltz/`
**Output:** `outputs/prescreen/prescreen_contacts.csv`

**Result:** Scanned 1330 models (950 AF2 + 190 AF3 + 190 Boltz) across 38 candidates. Most models bind CRD2; a minority show CRD1 or CRD1+CRD2 contacts.

---

### Step 02: Candidate Filtering

**Script:** `02_filter_candidates.py`

Reads the prescreen CSV and applies a CRD1 fraction >= 0.5 threshold per model. Counts how many sources (AF2, AF3, Boltz) have at least one passing model for each candidate. Keeps candidates that pass in `--min-sources` or more sources (default 2). Categorizes each candidate by source agreement.

**Input:** `outputs/prescreen/prescreen_contacts.csv`
**Output:**
- `outputs/prescreen/prescreen_summary.csv` -- per-candidate summary (models total, models passing, sources passing, keep/reject)
- `outputs/prescreen/prescreen_keep.txt` -- passing candidate IDs
- `outputs/prescreen/source_analysis.csv` -- per-source breakdown with categories (ALL_THREE, AF+BOLTZ, etc.)

**Result:** 21 candidates passed (6 ALL_THREE + 15 with 2/3 sources). 17 eliminated (15 single-source + 2 no CRD1 signal).

---

### Step 03: CIF to PDB Conversion

**Script:** `03_build_complex_pdb.py`

Converts the best CRD1-scoring CIF model per candidate per source into PDB format for HADDOCK3. Selects the model with the highest CRD1 fraction from each passing source. Standardizes chain IDs to A (antigen), H (heavy), L (light).

**Input:** CIF files + `prescreen_contacts.csv`
**Output:** `outputs/complex_pdb/{af,af2,boltz}/*.pdb`

**Result:** 48 PDB files (20 AF + 7 AF2 + 21 Boltz). Only sources that passed CRD1 fraction >= 0.5 get a PDB. AF2 chain remapping: ranked/relaxed CIFs (A/B/C -> H/L/A), unrelaxed CIFs (B/C/D -> H/L/A).

---

### Step 04: HADDOCK3 Run Setup

**Script:** `04_prepare_haddock_runs.py`

Splits each complex PDB into separate chain files (antigen.pdb, antibody_H.pdb, antibody_L.pdb) and generates antibody-aware restraint files:
- **ambig.tbl** -- CDR-epitope Ambiguous Interaction Restraints (AIRs) between CDR loop residues and CRD1 epitope residues within 5A
- **unambig.tbl** -- H-L chain anchoring restraints to prevent non-physical chain separation during simulation

Generates a HADDOCK3 TOML config (`refine.cfg`) for each run with the pipeline: topoaa -> flexref -> emref -> caprieval -> contactmap -> prodigyprotein.

**Input:** `outputs/complex_pdb/{af,af2,boltz}/*.pdb`
**Output:**
- `outputs/haddock_refine/CANDIDATE/SOURCE/refine.cfg`
- `outputs/haddock_refine/CANDIDATE/SOURCE/{antigen,antibody_H,antibody_L}.pdb`
- `outputs/haddock_refine/CANDIDATE/SOURCE/{ambig,unambig}.tbl`
- `outputs/haddock_refine/run_dirs.txt`

**Result:** 48 refinement runs prepared (one per candidate/source PDB), each with CDR-epitope AIRs and H-L anchoring restraints.

---

### HADDOCK3 Refinement

**Script:** `run_haddock_refine_batch.sh`

Executes HADDOCK3 refinement for all prepared runs. Reads `run_dirs.txt`, runs `haddock3 refine.cfg` in each directory. Skips already-completed runs (checks for `*_caprieval` directory).

**Input:** `outputs/haddock_refine/run_dirs.txt`
**Output:** HADDOCK3 run directories with refined models, energy scores, contact maps, and PRODIGY binding affinity.

**Result:** All 48 refinement runs completed. Each produced numbered module directories (0_topoaa/, 1_flexref/, 2_mdref/ or 2_emref/, 3_caprieval/, 4_contactmap/, 5_prodigyprotein/).

---

### Step 05: Score Collection

**Script:** `05_collect_haddock_scores.py`

Parses HADDOCK3 output (`capri_ss.tsv`) from the last caprieval directory in each run. Extracts HADDOCK score, Evdw, Eelec, Edesolv, and BSA. Selects the best-scoring (rank 1) model per candidate/source.

**Input:** HADDOCK3 run directories via `run_dirs.txt`
**Output:** `outputs/reports/refine_scores.csv` -- best model per candidate/source (48 rows)

**Result:** Energy scores collected for all 48 refined models across 21 candidates.

---

### Step 06: Rosetta Interface Analysis

**Script:** `06_rosetta_interface_analysis.py`

Runs PyRosetta InterfaceAnalyzerMover and AntibodyInfo on each HADDOCK3-refined PDB (rank 1 model). Extracts:
- Shape complementarity (Sc)
- Packing quality (packstat)
- Hydrogen bond count
- Salt bridge count (charged atom pairs within 4A)
- Per-CDR loop contact decomposition (H1-H3, L1-L3)
- H3 kink type classification
- CRD1 vs CRD2 contact mapping

**Input:** HADDOCK3 refined PDBs (emref outputs)
**Output:** `outputs/reports/rosetta_interface.csv` -- one row per candidate/source (48 rows)

**Result:** Interface characterisation for all 48 models. Shape complementarity range: 0.50-0.82. H-bond range: 3-25.

---

### Step 07: Select Best Models

**Script:** `07_select_best_models.py`

Merges HADDOCK energy scores with Rosetta interface metrics. For each candidate, selects the best source model by most negative HADDOCK score, then attaches all Rosetta columns.

**Input:**
- `outputs/reports/refine_scores.csv` (HADDOCK energetics)
- `outputs/reports/rosetta_interface.csv` (Rosetta characterisation)

**Output:** `outputs/reports/combined_scores.csv` -- one row per candidate (21 rows)

**Result:** 21 best models selected (one per candidate). Sources: 3 AF2, 4 AF3, 14 Boltz.

---

### Step 08: Final Ranking

**Script:** `08_final_ranking.py`

Assigns CRD binding category (CRD1, CRD1+CRD2, CRD2, NEITHER) based on Rosetta contact fractions (threshold: 10%). Ranks candidates by HADDOCK score (most negative first).

**Input:** `outputs/reports/combined_scores.csv`
**Output:** `outputs/reports/final_ranking.csv` -- ranked list with CRD categories

**Result:** 21 candidates ranked. 16 classified as CRD1-specific (>90% CRD1 contacts), 5 as CRD1+CRD2 (12-17% CRD2). Top 5: Ab-36.1, Ab-9.1, Ab-26.1, Ab-1.1, Ab-23.1.

---

## Notes

### Model counts

1330 total across 38 candidates. AF2 gives 25 per candidate (5 architectures x 5 seeds, ranked models). AF3 gives 5 diffusion samples per candidate. Boltz gives 5 diffusion samples per candidate.

### CRD fraction vs CRD coverage

Fraction = what proportion of the antibody's total antigen contacts land on CRD1. Tells you specificity. A fraction of 1.0 means the antibody only touches CRD1. A fraction of 0.5 means half the contacts are CRD1, half elsewhere.

Coverage = how many of the 42 CRD1 residues are contacted. Tells you breadth. Coverage of 0.52 means 22 out of 42 CRD1 residues are touched.

### Are any candidates pure CRD1 binders?

Yes. 9 of the 21 candidates have 100% CRD1 contacts (zero CRD2) in their best refined model: Ab-26.1, Ab-6.1, Ab-34.1, Ab-37.1, Ab-15.1, Ab-12.1, Ab-16.1, Ab-30.1, Ab-28.1. An additional 7 are CRD1-dominant (93-97% CRD1), bringing the total CRD1-specific count to 16. The remaining 5 candidates contact CRD2 at 12-17%.

### Cross-source agreement

6 candidates have CRD1 fraction >= 0.5 in all 3 sources (AF2, AF3, Boltz): Ab-12.1, Ab-25.1, Ab-29.1, Ab-30.1, Ab-36.1, Ab-38.1.

5 candidates have at least one pure CRD1 model (zero CRD2 contacts) in 2+ sources: Ab-12.1, Ab-30.1, Ab-36.1, Ab-37.1, Ab-38.1.

4 overlap between both lists: Ab-12.1, Ab-30.1, Ab-36.1, Ab-38.1.

### Why cross-source agreement over seed count?

Seed count within a single method is noise. Cross-source agreement is the real signal -- if AF3 and Boltz both independently predict CRD1 for the same candidate, two different algorithms arrived at the same answer.

### Why no pLDDT/PAE confidence filtering?

HADDOCK3 refinement handles quality assessment through energy scoring. Structures with poor geometry get poor HADDOCK scores, so a separate confidence filter is redundant.

### Why 21 candidates?

CRD1 fraction >= 0.5 in 2 or more sources. This eliminates 17 candidates -- 15 that only pass in 1 source and 2 with no CRD1 signal (Ab-5.1, Ab-7.1). The 2+ sources threshold gives 21 candidates, which is manageable for HADDOCK refinement.

### Why best model per source for HADDOCK?

Running all seeds through HADDOCK would be hundreds of runs on mostly CRD2-binding poses. Best-per-source gives up to 3 runs per candidate (best AF2, best AF3, best Boltz), 48 total. This lets us compare across methods while keeping the run count practical. Final ranking selects the best-scoring model per candidate.

---

## CRD Residue Sequences

CRD1 (CD40 20-62) -> positions 1-42
EPPTACREKQYLINSQCCSLCQPGQKLVSDCTEFTETECLPC

CRD2 (CD40 63-111) -> positions 43-91
GESEFLDTWNRETHCHQHKYCDPNLGLRVQQKGTSETDTICTCEEGWHC

CRD3 (CD40 112-161) -> positions 92-141
TSEACESCVLHRSCSPGFGVKQIATGVSDTICEPCPVGFFSNVSSAFEKC

CRD4 (CD40 162-193) -> positions 142-173
HPWTSCETKDLVVQQAGTNKTDVVCGPQDRLR
