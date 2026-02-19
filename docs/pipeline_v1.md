# Pipeline v1

Initial version of the CD40-CRD1 antibody screening pipeline (January 2026).

Integrated AlphaFold3 + Boltz predictions, CRD1-focused contact analysis, HADDOCK refinement and scoring, and final ranking by averaged HADDOCK scores.

---

## Step 01: Input Preparation

Consolidated AF3 and Boltz predictions into a single project workspace. Standardised chain IDs to A (antigen), H (heavy), L (light). Converted CIF to PDB where needed for refinement.

**Input:** AF3 predictions, Boltz predictions 
**Output:** Standardised PDB files

---

## Step 02: CRD1 Contact Prescreen (Specificity)

Scanned all models for antibody-antigen contacts at 5A heavy-atom cutoff. Computed CRD1 contact fraction (specificity) per model. Applied CRD1 fraction >= 0.5 threshold.



**Candidate-level results (passes if any model from that source passes):**
- AF pass: 25
- Boltz pass: 27
- Both AF + Boltz pass: 18
- AF-only: 7
- Boltz-only: 9
- Neither: 4

---

## Step 03: CRD1 Coverage (Breadth)

Computed CRD1 coverage = crd1_contact_residues / 42 for each model. Applied coverage >= 0.30 threshold (~13-15 residues on CRD1). This ensured the interface footprint spanned a meaningful portion of CRD1 rather than a small edge contact.

---

## Step 04: Local Confidence Filtering (Interface Quality)

Computed local pLDDT averaged over CRD1-contacting residues and their contacting antibody residues. Computed local PAE for CRD1-CDR contact pairs.

**Thresholds:**
- AF local pLDDT >= 67
- Boltz local pLDDT >= 70


---

## Step 05: HADDOCK Refinement and Scoring

Refined predicted complexes using HADDOCK refine-complex protocol with N_COMP=3 (chains A/H/L preserved). This was refinement of the predicted pose, not blind docking.

**HADDOCK stages used:**
- Stage it1 (semi-flexible refinement): interface-side flexibility and optimisation of the existing pose
- Stage water (explicit solvent refinement): final relaxation in water to remove clashes and refine energetics
- Stage it0 (rigid-body docking) was skipped since AF/Boltz already provided plausible binding poses

**HADDOCK score formula:**
```
HADDOCK score = 1.0*Evdw + 0.2*Eelec + 1.0*Edesolv + 0.1*Eair
```
In refinement-only runs without AIR restraints, Eair was typically 0. Score dominated by Evdw, Eelec, and Edesolv.

**Metrics extracted:** HADDOCK score, Evdw, Eelec, Edesolv, BSA (FreeSASA), H-bonds, salt bridges.

---

## Step 06: Binding-Mode RMSD vs 6FAX

Computed binding-mode RMSD by aligning the antigen to the 6FAX crystal structure (a CRD1 reference) and measuring antibody orientation displacement. Used to compare binding modes across candidates.

---

## Step 07: Final Ranking

Applied combined gate: CRD1 fraction + coverage + local pLDDT must pass in **both** AF and Boltz. Ranked by HADDOCK average score (AF + Boltz), with Boltz local PAE as tie-breaker / confidence context.

---

## Results

**10 candidates** passed all filters (CRD1 specificity + breadth + local pLDDT in both AF and Boltz).

| Rank | Candidate | CRD1 cov (AF) | CRD1 cov (Boltz) | pLDDT (AF) | pLDDT (Boltz) | PAE (Boltz) | HADDOCK avg | BSA (Boltz) |
|------|-----------|---------------|-------------------|------------|---------------|-------------|-------------|-------------|
| 1 | Ab-38.1 | 0.33 | 0.50 | 68.0 | 92.6 | 8.8 | -297.5 | 3640 |
| 2 | Ab-36.1 | 0.38 | 0.50 | 69.2 | 84.8 | 11.9 | -279.6 | 3508 |
| 3 | Ab-24.1 | 0.48 | 0.55 | 69.2 | 79.4 | 9.8 | -279.6 | 3888 |
| 4 | Ab-3.1 | 0.33 | 0.33 | 68.4 | 77.5 | 12.3 | -277.7 | 3340 |
| 5 | Ab-28.1 | 0.52 | 0.45 | 67.8 | 70.5 | 12.5 | -268.1 | 3726 |
| 6 | Ab-25.1 | 0.38 | 0.40 | 71.0 | 75.9 | 20.5 | -267.8 | 3398 |
| 7 | Ab-29.1 | 0.48 | 0.43 | 69.0 | 81.5 | 8.6 | -260.3 | 3508 |
| 8 | Ab-35.1 | 0.31 | 0.43 | 72.4 | 73.2 | 7.1 | -253.8 | 3426 |
| 9 | Ab-34.1 | 0.38 | 0.38 | 72.3 | 82.0 | 20.0 | -250.0 | 3424 |
| 10 | Ab-12.1 | 0.31 | 0.33 | 73.0 | 89.7 | 4.7 | -245.6 | 3084 |

---

## Output Files

- `outputs/prescreen/prescreen_contacts.csv` -- per-model contact metrics
- `outputs/prescreen/prescreen_summary.csv` -- per-candidate summary
- `outputs/ranks/final_ranking.csv` -- ranked candidates with full metrics
- `outputs/ranks/final_report_simple.csv` -- simplified ranking table
