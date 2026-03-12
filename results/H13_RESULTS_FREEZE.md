# H13 Results Freeze: From Trehalose Vulnerability to AKT-Sensitive Metabolic Subtype

**Date:** 2026-03-12 (updated with H13g convergence analysis)
**Status:** COMPLETE — Ready for paper drafting
**Pipeline:** H13a → H13b → H13c → H13d → H13e → H13f → H13g → H13h

---

## Executive Summary

H13 began as a trehalose/autophagy vulnerability hypothesis and, through systematic
falsification and pivoting, converged on a cleaner and stronger finding:

> **SLC2A8 (GLUT8) expression marks a metabolic subtype of KRAS-mutant cancers
> that is selectively sensitive to AKT inhibitors (capivasertib, ipatasertib,
> afuresertib), independent of tissue lineage.**

The evidence chain spans four independent data layers: unbiased drug screening
(PRISM), lineage-controlled pharmacological validation, transcriptional
co-expression, and CRISPR genetic dependencies. Convergence testing (H13g)
confirmed that the biomarker finding is robust but that the lipogenesis
dependencies represent independent features of the subtype, not a unified
targetable pathway — clarifying the claim and strengthening the paper's
scientific integrity.

---

## 1. Discovery Arc

### H13a — Initial Hypothesis Testing

**Question:** Do SLC2A8-high / RAS-mutant cells show enriched autophagy/lysosomal
dependencies?

**Results:**
- F1 PASS: 69% of PDAC lines express SLC2A8 above pan-cancer 25th percentile
- F2 PASS: Unsupervised clustering yields coherent k=2 subgroups (silhouette=0.355)
- F3 FAIL: No gene set (autophagy/lysosomal/ER stress) reached set-level significance
- F4 PASS: Within-PDAC stratification significant (p=0.0215)

**Verdict:** PASS (3/4) but F3 failure motivated deeper investigation.

### H13b — Autophagy Initiation Deep Dive

**Question:** Does splitting autophagy into initiation vs execution reveal a masked signal?

**Key Results:**
- F3b PASS: Initiation enrichment p=3.5e-05 vs execution p=0.70 (20,000x separation)
- Driver genes: PIK3C3/VPS34 (p=1.4e-05), ATG13 (p=0.039), BECN1 (p=0.039)
- F5 PASS: PRISM drug sensitivity — top compounds near p<0.05
- F6 FAIL: TCGA-PAAD survival — bootstrap HR=1.020 (CI includes 1.0)
- F7 FAIL: Signal is pan-GI, not PDAC-specific

**Verdict:** PARTIAL (2/4). PIK3C3 identified as dominant driver gene.

### H13c — Unbiased PRISM Screen + MOA Enrichment

**Question:** What compounds actually kill SLC2A8-high/RAS-mut cells?

**Key Results:**
- 6,790 compounds tested. Only 2 survive FDR < 0.10:
  - Alpelisib (PI3K inhibitor): p=3.77e-06
  - Afuresertib/GSK2110183 (AKT inhibitor): p=5.20e-06
- MOA enrichment (Fisher exact test):
  - AKT inhibitor: 10/25 sig, OR=19.93, p=3.28e-09
  - PI3K inhibitor: 13/66 sig, OR=7.39, p=2.18e-07
  - MEK inhibitor: 5/26 sig, OR=6.95, p=1.54e-03
  - mTOR inhibitor: 1/38 sig, OR=0.77, p=0.73 (NOT enriched)
- VPS34-specific inhibitors (SAR405, PIK-93): NOT significant
- Autophagy inducer control: null (STF-62247 p=0.95)

**Verdict:** Drug sensitivity is PI3K/AKT, not autophagy/VPS34.

### H13d — PI3K/AKT Pathway Validation

**Question:** Do target cells show elevated PI3K/AKT pathway expression?

**Key Results:**
- AKT1: d=+0.522, p=2.58e-11 (HIGHER in target)
- PDPK1: d=+0.420, p=1.17e-07
- PIK3CB: d=+0.360, p=1.31e-05
- MTOR: d=+0.200, p=0.012
- RPS6KB1: d=+0.155, p=0.029
- PTEN: NOT different (d=-0.096, p=0.20)
- SLC2A8 vs AKT1: rho=+0.331, p=1.51e-29
- SLC2A8 vs MTOR: rho=+0.284, p=7.05e-22

**But:** PI3K pathway expression predicts RESISTANCE to PI3K drugs (rho=+0.18),
not sensitivity. Expression ≠ mechanism.

### H13e — Sanity Checks (Critical)

**Three checks that decomposed the signal:**

1. **PRISM orientation confirmed:** Bortezomib 100% negative. Lower = more killing.

2. **Within-group correlation:** PI3K pathway paradox partially resolves. But SLC2A8
   expression directly predicts AKT drug sensitivity:
   - GSK2110183: rho=-0.165, p=4.0e-04
   - GDC-0068: rho=-0.141, p=2.4e-03
   - AZD5363: rho=-0.137, p=3.1e-03
   - PI3K inhibitors: null (alpelisib rho=-0.044, p=0.34)

3. **Lineage control (decisive):** After regressing out lineage:
   - AKT inhibitors SURVIVE: GSK2110183 partial r=-0.144 (p=0.002),
     GDC-0068 r=-0.147 (p=0.002), AZD5363 r=-0.131 (p=0.005)
   - PI3K inhibitors GO NULL: alpelisib r=-0.011 (p=0.82)
   - PDAC lineage alone predicts alpelisib sensitivity (d=-0.55, p=0.002)

**Verdict:** Two independent signals were conflated in H13c. The real biomarker
finding is SLC2A8 → AKT inhibitor sensitivity (lineage-independent).

### H13f — Metabolic Mechanism

**Question:** Why does SLC2A8 predict AKT inhibitor sensitivity?

**Co-expression fingerprint (Analysis 1):**
- Endosomal trafficking: VPS35 rho=+0.36, RAB7A +0.19, LAMP1 +0.29
- Lipid biosynthesis: FASN rho=+0.29 (d=+0.39), SREBF1 +0.18 (d=+0.39),
  ACACA +0.25 (d=+0.27)
- Mid-glycolysis: PFKL +0.33 (d=+0.46), PKM +0.21 (d=+0.36)
- NOT Warburg: HK2, LDHA, GAPDH all null. PDK1 LOWER (d=-0.21)

**Drug predictors (Analysis 2):**
- SLC2A4/GLUT4 matches SLC2A8 as AKT drug predictor (avg rho ≈ -0.15)
- SCD1 also predicts (avg rho ≈ -0.11)
- AKT1 and MTOR expression have ZERO predictive power (rho ≈ 0)

**Head-to-head (Analysis 3):**
- SLC2A8 alone: best predictor (rho=-0.165 for GSK2110183)
- Lysosomal composite: diluted by MCOLN1/LAMTOR1 (predict resistance)
- PI3K pathway: null

**CRISPR dependencies (Analysis 4):**
- SCD: d=-0.316, p=5.5e-06 (as strong as PIK3C3)
- FASN: d=-0.173, p=0.019
- AKT1: d=-0.277, p=4.5e-05
- PIK3C3: d=-0.302, p=7.5e-05
- Lysosomal machinery: NO differential dependency

### H13g — Convergence Testing (Falsification of Mechanistic Loop)

**Question:** Do AKT drug sensitivity and lipid gene CRISPR dependency co-occur
at the individual cell-line level? Does the SLC2A8 → AKT → lipogenesis pipeline
function as a single targetable chain?

**Key Results:**

1. **Cell-line-level convergence: NULL.** SCD CRISPR dependency shows zero
   correlation with AKT drug sensitivity, both pan-cancer and within the
   SLC2A8-high/RAS-mut target group:
   - SCD vs afuresertib (target only): rho=-0.030, p=0.72
   - SCD vs ipatasertib (target only): rho=+0.013, p=0.88
   - SCD vs capivasertib (target only): rho=+0.064, p=0.44
   - FASN: equally null across all drugs (all p>0.4)

2. **Interaction model: SLC2A8 is the sole predictor.** In the regression model
   `drug_sens ~ SLC2A8 + SCD_dep + SLC2A8*SCD_dep`:
   - SLC2A8: significant for all 3 AKT drugs (p=0.0003 to 0.003)
   - SCD dependency: null (p=0.27 to 0.72)
   - Interaction term: null (p=0.70 to 0.997)
   - Same pattern with FASN dependency

3. **SCD inhibitor PRISM data: weak and contradictory.**
   - A-939572 (SCD inh): diff=-0.122, p=0.071 (trending, not significant)
   - Aramchol (SCD inh): diff=+0.113, p=0.014 (WRONG direction — target resistant)
   - FASN/ACC inhibitors: all null

4. **Quadrant analysis: group effect, not synergy.** Target cells are enriched
   in the dual-vulnerability quadrant (AKT-sens + SCD-dep):
   - Median split: OR=2.32, p=0.0002
   - 25th percentile: OR=4.52, p=0.0002
   - But target cells are DEPLETED from SCD-dep-only (OR=0.52, p=0.008)
   - Interpretation: independent marginal enrichments, not cell-level synergy

5. **Lineage-controlled SLC2A8 vs lipid dependency: ALL NULL.**
   - SCD: partial r=+0.009, p=0.77
   - FASN: partial r=+0.018, p=0.55
   - ACLY, ACACA, SREBF1: all p>0.37
   - Transcriptional co-expression is real; functional dependency is not linked

**Verdict:** The metabolic dependency loop does NOT close at the cell-line level.
SLC2A8-high/RAS-mut cells show both AKT sensitivity and lipogenesis dependency
as independent features of the metabolic subtype, not as a unified causal chain.
This falsifies the "synthetic lethal metabolic therapy" hypothesis but strengthens
the core biomarker claim by clarifying exactly what SLC2A8 does and does not predict.

### H13h — Reviewer Analyses (RAS Interaction, Classification, TCGA)

**Question:** Is the RAS mutation context biologically meaningful? How does SLC2A8
perform as a classifier? Do co-expression patterns replicate in patient tumors?

**Key Results:**

1. **RAS mutation interaction: SIGNIFICANT.** SLC2A8:KRAS interaction term is
   significant for all 3 AKT drugs:
   - Afuresertib: β=-0.216, p=0.004
   - Ipatasertib: β=-0.153, p=0.040
   - Capivasertib: β=-0.211, p=0.005
   - SLC2A8 *alone* is NULL (p=0.69, 0.56, 0.997)
   - F-test for KRAS terms: p=0.0003, 0.026, 0.004
   - **SLC2A8 predicts AKT sensitivity only within the RAS-mutant context**

2. **Classification performance: modest but consistent.**
   - ROC AUC (median split): 0.550-0.593
   - Top vs bottom quartile: Cohen's d=-0.37 to -0.47, all p<0.006
   - Top-quartile SLC2A8 shows 0.29-0.37 greater drug sensitivity

3. **TCGA-PAAD partial replication (n=178).**
   - SLC2A8 vs AKT1: rho=+0.212, p=0.005 (replicates)
   - SLC2A8 vs FASN: rho=+0.180, p=0.016 (replicates)
   - SLC2A8 vs VPS35: rho=-0.171, p=0.023 (OPPOSITE direction from DepMap)
   - SLC2A8 vs MTOR: rho=-0.022, p=0.77 (null)
   - SLC2A8 vs SCD: rho=+0.031, p=0.68 (null)
   - Distribution differs from DepMap PDAC (KS p<0.0001) — expected due to TME

**Verdict:** RAS context is biologically meaningful (interaction, not filter).
Classification is modest but consistent with single-gene pharmacogenomic biomarkers.
TCGA partially replicates key co-expression patterns (AKT1, FASN) despite bulk
tumor RNA-seq limitations.

---

## 2. The Final Model

```
SLC2A8-high expression (GLUT8, lysosomal glucose transporter)
    ↕ marks a metabolic tumor state characterized by
Co-expression of:
  ├─ Endosomal trafficking (VPS35/retromer, RAB7A, LAMP1)
  ├─ Lipid biosynthesis (FASN, SCD, ACACA, SREBF1)
  └─ Oxidative glycolysis (PFKL, PKM, GPI; NOT Warburg)
    ↕ this state is
Selectively sensitive to AKT inhibition
  (capivasertib, ipatasertib, afuresertib)
    ↕ confirmed by
Lineage-independent pharmacological correlation (partial r=-0.131 to -0.147)
```

**What the model IS:** SLC2A8 is a biomarker that identifies a metabolic
subtype selectively sensitive to AKT inhibitors.

**What the model is NOT:** A unified causal chain (SLC2A8 → AKT → lipogenesis
→ metabolic collapse). H13g showed that AKT drug sensitivity and lipogenesis
dependency are independent features of the subtype, not a single targetable
pathway. The co-expression pattern describes the metabolic state; the drug
sensitivity is the actionable finding.

The PIK3C3/VPS34 CRISPR dependency (from H13b) and SCD/FASN dependencies
(from H13f) are real genetic features of the SLC2A8-high subtype but do not
predict AKT drug sensitivity at the individual cell-line level (H13g).
They characterize the biology of the subtype without forming a druggable
causal chain.

---

## 3. Key Statistics (Citable)

### Biomarker

| Metric | Value |
|--------|-------|
| SLC2A8 vs AKT drug sensitivity (GSK2110183) | rho=-0.165, p=4.0e-4, n=457 |
| SLC2A8 vs AKT drug sensitivity (GDC-0068) | rho=-0.141, p=2.4e-3, n=461 |
| SLC2A8 vs AKT drug sensitivity (AZD5363) | rho=-0.137, p=3.1e-3, n=461 |
| Lineage-controlled partial r (GSK2110183) | r=-0.144, p=0.002 |
| Lineage-controlled partial r (GDC-0068) | r=-0.147, p=0.002 |
| Lineage-controlled partial r (AZD5363) | r=-0.131, p=0.005 |

### Co-expression (SLC2A8 Spearman, n=1103)

| Gene | rho | p |
|------|-----|---|
| VPS35 (retromer) | +0.364 | 7.9e-36 |
| AKT1 | +0.331 | 1.5e-29 |
| PFKL (glycolysis) | +0.335 | 3.0e-30 |
| SLC2A6 (lysosomal GLUT) | +0.323 | 3.8e-28 |
| FASN (lipogenesis) | +0.288 | 1.7e-22 |
| MTOR | +0.284 | 7.1e-22 |
| LAMP1 | +0.287 | 2.2e-22 |
| ACACA (lipogenesis) | +0.249 | 4.7e-17 |
| SREBF1 (lipogenesis TF) | +0.178 | 2.5e-09 |
| RPS6KB1 | +0.163 | 5.1e-08 |

### CRISPR Dependencies (target vs background, Cohen's d)

| Gene | d | p |
|------|---|---|
| SCD (lipogenesis) | -0.316 | 5.5e-06 |
| PIK3C3 (VPS34) | -0.302 | 7.5e-05 |
| AKT1 | -0.277 | 4.5e-05 |
| FASN (lipogenesis) | -0.173 | 0.019 |

### Unbiased PRISM MOA Enrichment (6,790 compounds)

| MOA | Sig/Total | OR | p |
|-----|-----------|----|----|
| AKT inhibitor | 10/25 | 19.93 | 3.3e-09 |
| PI3K inhibitor* | 13/66 | 7.39 | 2.2e-07 |
| MEK inhibitor | 5/26 | 6.95 | 1.5e-03 |

*PI3K signal is lineage-confounded; AKT signal survives lineage correction.

---

## 4. What This Is NOT

- NOT a trehalose/autophagy finding (original hypothesis falsified at H13a/H13c)
- NOT PI3K pathway addiction (pathway expression predicts resistance, not sensitivity)
- NOT PDAC-specific for drug sensitivity (AKT signal is lineage-independent)
- NOT a causal mechanism (expression ≠ activation; correlation ≠ causation)
- NOT a synthetic-lethal metabolic therapy target (H13g: AKT sensitivity and
  lipogenesis dependency do not converge at the cell-line level)
- NOT a combination therapy rationale for AKT + SCD inhibitors (no interaction
  effect, SCD inhibitor PRISM data weak/contradictory)

## 5. What This IS

- A **biomarker hypothesis**: SLC2A8 marks a metabolic state sensitive to AKT inhibitors
- A **metabolic subtype description**: co-expression of endosomal trafficking,
  lipogenesis, and oxidative glycolysis programs
- A **lineage-independent signal**: survives correction for tissue type
- A **translatable finding**: AKT inhibitors are FDA-approved / in clinical trials
- A **scientifically honest claim**: tested the mechanistic causal chain (H13g)
  and reported the negative result, narrowing the claim to what the data supports

## 6. Experimental Validation Roadmap

| Priority | Experiment | Tests | Note |
|----------|-----------|-------|------|
| 1 | SLC2A8-high vs low cell lines + AKT inhibitors | Biomarker validation | Core claim |
| 2 | Patient tumor SLC2A8 IHC + AKT inhibitor response | Clinical translation | Core claim |
| 3 | AKT inhibition → lipogenesis measurement | Subtype characterization | Context, not prediction |
| 4 | SCD + AKT inhibitor combination | Combination testing | Exploratory only — H13g shows no cell-level convergence; test empirically |
| 5 | VPS35 knockdown → AKT signaling | Endosomal recycling model | Mechanistic, not translational |

---

## 7. Data Provenance

| Resource | Version | Files |
|----------|---------|-------|
| DepMap | 24Q4 | Expression, CRISPR, Mutations, Model metadata |
| PRISM | 24Q2 Extended | Sensitivity matrix, Compound list |
| TCGA-PAAD | GDC current | Expression + clinical (H13b F6 only) |

All analyses use the same target/background split: SLC2A8-high (above median) +
RAS-mutant (KRAS/NRAS/HRAS/BRAF/PIK3CA) vs all others, across 1,103 cell lines
with matched expression + CRISPR + metadata.

---

## 8. Scripts

| Script | Phase | What it does |
|--------|-------|-------------|
| `scripts/h13c_deep_dive.py` | H13c | Unbiased PRISM screen, MOA enrichment, per-gene drivers |
| `scripts/h13d_pi3k_validation.py` | H13d | PI3K/AKT pathway expression validation |
| `scripts/h13e_sanity_checks.py` | H13e | PRISM orientation, within-group correlation, lineage control |
| `scripts/h13f_metabolic_clue.py` | H13f | Metabolic co-expression, drug predictors, CRISPR deps |
| `scripts/h13g_convergence.py` | H13g | Cell-level convergence testing, interaction models, quadrant analysis |
| `scripts/h13h_reviewer_analyses.py` | H13h | RAS interaction model, ROC AUC, TCGA descriptive |

Earlier phases (H13a, H13b) are in `src/pdac/h13_trehalose_vulnerability/`.

---

*Results frozen 2026-03-12 (H13a-f). Updated 2026-03-12 with H13g convergence
analysis. Scientific arc complete: discovery → validation → falsification of
overclaim → refined biomarker hypothesis. Ready for paper drafting.*
