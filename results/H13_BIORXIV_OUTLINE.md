# bioRxiv Preprint Outline

**Working Title:**
SLC2A8 expression identifies a metabolic tumor subtype selectively sensitive to AKT inhibitors across KRAS-mutant cancers

**Authors:** [TBD]
**Target:** bioRxiv (Cancer Biology section)
**Format:** Research Article

---

## Abstract (~250 words)

**Background:** AKT inhibitors (capivasertib, ipatasertib) are entering clinical use but lack robust predictive biomarkers beyond PIK3CA mutation and PTEN loss. Identifying patients who benefit requires new stratification approaches.

**Methods:** We performed an integrative computational analysis across 1,103 cancer cell lines using DepMap (24Q4) expression, CRISPR dependency, and mutation data combined with PRISM (24Q2) drug sensitivity screening of 6,790 compounds. We stratified cell lines by SLC2A8 (GLUT8) expression and RAS-pathway mutation status, then tested for differential drug sensitivity, transcriptional co-expression, and genetic dependencies.

**Results:** Unbiased PRISM screening identified AKT inhibitors as the most enriched drug class among SLC2A8-high/RAS-mutant-sensitive compounds (OR=19.93, p=3.3e-9). SLC2A8 expression directly predicted AKT inhibitor sensitivity (Spearman rho=-0.165, p=4.0e-4 for afuresertib), and this signal survived lineage correction (partial r=-0.144, p=0.002), unlike PI3K inhibitor signals which were confounded by tissue type. The SLC2A8-high subtype co-expressed endosomal trafficking (VPS35, RAB7A, LAMP1) and lipid biosynthesis programs (FASN, SCD, SREBF1), and showed differential CRISPR dependencies on SCD (d=-0.316, p=5.5e-6) and AKT1 (d=-0.277, p=4.5e-5). However, convergence testing revealed that these genetic dependencies and AKT drug sensitivity are independently distributed at the cell-line level, indicating that SLC2A8 marks a metabolic state rather than a single targetable pathway.

**Conclusions:** SLC2A8 expression identifies a pan-cancer metabolic subtype selectively sensitive to AKT inhibitors, independent of tissue lineage. This biomarker hypothesis has immediate translational relevance given that capivasertib is already FDA-approved.

---

## 1. Introduction (~800 words)

### Paragraph 1: Clinical problem
- AKT inhibitors are entering oncology practice (capivasertib approved 2023 for breast cancer)
- Current biomarkers (PIK3CA, AKT1, PTEN) enrich for response but miss many responders
- KRAS-mutant cancers (PDAC, colorectal, lung) lack effective targeted therapies
- Need for biomarkers that cut across lineage and genotype

### Paragraph 2: SLC2A8/GLUT8 biology
- SLC2A8 (GLUT8) is a type III facilitative glucose transporter localized to lysosomes and endosomes
- Unlike plasma membrane GLUTs (GLUT1, GLUT4), GLUT8 transports glucose across intracellular membranes
- Functions in trehalose metabolism in model organisms; role in cancer metabolism poorly characterized
- Recent evidence links lysosomal glucose metabolism to mTOR signaling and autophagy
- GLUT8 is highly expressed in RAS-mutant gastrointestinal cancers (cite H13a: 69% of PDAC lines above pan-cancer 25th percentile)

### Paragraph 3: Rationale and approach
- We initially hypothesized that SLC2A8-high/RAS-mutant cells would show autophagy/lysosomal vulnerability (hypothesis-driven)
- Through systematic testing and falsification, we discovered a different and stronger signal: selective sensitivity to AKT inhibitors
- This paper presents the full evidence chain: unbiased drug screening, lineage-controlled pharmacological validation, transcriptional co-expression, and CRISPR genetic dependencies
- We also tested whether the observed vulnerabilities form a single targetable metabolic pathway; they do not, clarifying the biomarker claim

### Paragraph 4: Summary of findings (brief)
- SLC2A8 marks a metabolic subtype selectively sensitive to AKT inhibitors
- The signal is lineage-independent and validated across four independent data layers
- The subtype is characterized by distinct biology (endosomal trafficking, lipogenesis) but the actionable finding is the AKT drug sensitivity

**Key references to cite:**
- Capivasertib approval (Turner et al., NEJM 2023)
- GLUT8 biology (Mueckler & Thorens, Mol Aspects Med 2013; DeBosch et al., various)
- DepMap/PRISM methodology (Tsherniak et al., Cell 2017; Corsello et al., Nat Cancer 2020)
- AKT/SREBP1/lipogenesis axis (Porstmann et al., PNAS 2008; Li et al., PNAS 2010)
- Retromer/VPS35 and receptor recycling (Seaman, J Cell Sci 2012)

---

## 2. Results (~2800 words)

### 2.1 SLC2A8 expression defines a coherent subgroup within RAS-mutant cancers
- **Figure 1A**: SLC2A8 expression by lineage (violin plot, PDAC highlighted)
- 69% of PDAC lines express SLC2A8 above pan-cancer 25th percentile
- Unsupervised clustering yields k=2 subgroups (silhouette=0.355)
- SLC2A8-high/RAS-mutant stratification is significant within PDAC (p=0.0215)

### 2.2 Unbiased drug screening identifies AKT inhibitors as the dominant sensitivity
- **Figure 2**: MOA enrichment barplot
- 6,790 PRISM compounds tested for differential sensitivity
- Only 2 compounds survive FDR<0.10: alpelisib (PI3K) and afuresertib (AKT)
- MOA enrichment: AKT inhibitor class OR=19.93 (p=3.3e-9), the strongest enrichment
- VPS34-specific inhibitors (SAR405, PIK-93): not significant
- Autophagy inducers: null (STF-62247 p=0.95)
- mTOR inhibitors: not enriched (OR=0.77, p=0.73)

### 2.3 PI3K inhibitor sensitivity is a lineage confound; AKT signal is genuine
- **Figure 3**: Lineage decomposition (raw vs corrected correlations)
- SLC2A8 expression predicts AKT drug sensitivity continuously:
  - Afuresertib (GSK2110183): rho=-0.165, p=4.0e-4
  - Ipatasertib (GDC-0068): rho=-0.141, p=2.4e-3
  - Capivasertib (AZD5363): rho=-0.137, p=3.1e-3
- PI3K inhibitors show no continuous correlation (alpelisib: rho=-0.044, p=0.34)
- Lineage-controlled partial correlation:
  - AKT inhibitors SURVIVE: partial r=-0.131 to -0.147 (all p<0.005)
  - PI3K inhibitors GO NULL: alpelisib partial r=-0.011 (p=0.82)
- PDAC lineage alone predicts alpelisib sensitivity (d=-0.55, p=0.002) -- a confound
- PI3K pathway expression paradoxically predicts RESISTANCE to PI3K drugs (rho=+0.18)
- Conclusion: two independent signals were conflated; the true biomarker is SLC2A8 -> AKT sensitivity

### 2.4 SLC2A8-high cells co-express an endosomal trafficking and lipogenic program
- **Figure 4**: Metabolic co-expression fingerprint
- Endosomal/lysosomal trafficking: VPS35 rho=+0.36, LAMP1 +0.29, RAB7A +0.19
- Lipid biosynthesis: FASN rho=+0.29, ACACA +0.25, SREBF1 +0.18
- Mid-glycolysis (oxidative, NOT Warburg): PFKL +0.33, PKM +0.21, GPI +0.21
- Warburg genes null: HK2, LDHA, GAPDH show no correlation
- PDK1 is LOWER in target (d=-0.21), suggesting oxidative rather than fermentative metabolism
- SLC2A4/GLUT4 co-expression supports AKT->AS160->GLUT4 translocation axis
- Framing: this characterizes the metabolic state that SLC2A8 marks

### 2.5 CRISPR dependencies confirm distinct biology of the SLC2A8-high subtype
- **Figure 5**: CRISPR dependency box plots
- SCD: d=-0.316, p=5.5e-6 (the strongest metabolic dependency)
- PIK3C3/VPS34: d=-0.302, p=7.5e-5
- AKT1: d=-0.277, p=4.5e-5
- FASN: d=-0.173, p=0.019
- Lysosomal machinery (MCOLN1, LAMTOR1, CTSD): NO differential dependency
- These dependencies characterize the biology of the subtype

### 2.6 SLC2A8 outperforms pathway-level and composite biomarkers
- SLC2A8 alone: best AKT drug predictor (rho=-0.165)
- Lysosomal composite score: diluted by MCOLN1/LAMTOR1 (predict resistance)
- PI3K pathway composite: null predictive power
- AKT1 and MTOR expression: zero predictive power despite being co-expressed
- Expression =/= activation: pathway expression marks the subtype but doesn't drive drug sensitivity
- SLC2A8 marks a downstream metabolic state rather than the signaling pathway itself

### 2.7 Convergence testing: drug sensitivity and genetic dependencies are independently distributed
- **Supplementary Figure S3**: Convergence analysis
- Tested whether AKT drug sensitivity and SCD/FASN CRISPR dependency co-occur
  at the individual cell-line level
- Cell-line-level correlation within target group: null
  - SCD vs afuresertib: rho=-0.030, p=0.72
  - SCD vs ipatasertib: rho=+0.013, p=0.88
  - SCD vs capivasertib: rho=+0.064, p=0.44
- Interaction model (drug_sens ~ SLC2A8 + SCD_dep + SLC2A8*SCD_dep):
  - SLC2A8: significant (p<0.001 for all 3 AKT drugs)
  - SCD dependency: null (p=0.27-0.72)
  - Interaction: null (p=0.70-0.997)
- SLC2A8-high/RAS-mut cells are enriched in dual-vulnerability quadrant
  (AKT-sens + SCD-dep, OR=2.32, p=0.0002) but depleted from SCD-dep-only
  quadrant (OR=0.52, p=0.008), indicating independent group-level features
- Interpretation: SLC2A8 marks a metabolic state with multiple independent
  biological features. The AKT drug sensitivity is the actionable finding;
  the lipogenesis dependencies describe the subtype biology without forming
  a single targetable chain

---

## 3. Discussion (~1200 words)

### Paragraph 1: Summary of key finding
- SLC2A8 identifies a lineage-independent metabolic subtype selectively sensitive to AKT inhibitors
- Four independent evidence layers converge on the biomarker (drug screen, lineage-corrected correlation, co-expression, CRISPR)
- This is a biomarker hypothesis, not a causal mechanism (clearly state)

### Paragraph 2: Metabolic subtype characterization
- **Figure 6**: Summary model diagram
- SLC2A8-high cells co-express endosomal trafficking, lipid biosynthesis, and oxidative glycolysis programs
- This is NOT Warburg metabolism (HK2, LDHA null; PDK1 lower)
- The co-expression and CRISPR data characterize the biology of the subtype but
  convergence testing (Section 2.7) showed these are independently distributed features,
  not a unified causal pathway
- The metabolic state may create a cellular context in which AKT signaling is
  particularly important, without requiring a direct mechanistic chain from
  lipogenesis to drug sensitivity

### Paragraph 3: Clinical context and translational relevance
- Capivasertib (AZD5363) is FDA-approved for PIK3CA/AKT1/PTEN-altered HR+/HER2- breast cancer
- Ipatasertib and afuresertib are in clinical trials for multiple tumor types
- Current biomarkers (PIK3CA mutation, PTEN loss) are necessary but not sufficient
- SLC2A8 could complement existing biomarkers as a transcriptomic predictor
- Particularly relevant for KRAS-mutant cancers (PDAC, colorectal) where AKT inhibitors are being tested in combination with MEK/KRAS inhibitors
- MEK inhibitor enrichment (OR=6.95, p=1.5e-3) suggests potential AKT+MEK combination rationale

### Paragraph 4: The lineage confound as a cautionary finding
- Our analysis initially conflated PI3K and AKT signals because PDAC cells are both SLC2A8-high and alpelisib-sensitive
- Lineage correction decomposed this into two independent signals
- This has implications for biomarker studies that don't control for tissue type
- The false-positive rate for lineage-confounded biomarkers in pan-cancer analyses is likely substantial
- Recommendation: all pan-cancer drug sensitivity biomarkers should be tested with lineage correction

### Paragraph 5: Limitations
- Computational study only; no wet-lab validation
- Expression =/= protein levels =/= activity (transcriptomic biomarker, not mechanistic proof)
- Cell line-based; tumor microenvironment effects not captured
- Effect sizes are modest (rho ~ -0.14 to -0.17); clinically meaningful threshold unknown
- PRISM sensitivity measures growth inhibition in 2D culture; in vivo efficacy may differ
- PTEN was not differentially expressed (d=-0.096, p=0.20), but PTEN protein loss (post-translational) is not captured by mRNA data
- Correlational design; causal inference requires experimental validation
- Convergence testing showed that AKT drug sensitivity and lipogenesis dependencies
  are independently distributed; SLC2A8 predicts AKT sensitivity, not a metabolic
  pathway that can be targeted at multiple points simultaneously

### Paragraph 6: Experimental validation roadmap
- Priority 1: Dose-response curves in SLC2A8-high vs low cell lines + AKT inhibitors (core biomarker validation)
- Priority 2: Patient tumor SLC2A8 IHC correlated with AKT inhibitor response data (clinical translation)
- Priority 3: Lipogenesis assays (14C-acetate incorporation) after AKT inhibition (subtype characterization)
- Priority 4: SCD inhibitor + AKT inhibitor combination (exploratory -- H13g shows no cell-level convergence; empirical testing needed)

---

## 4. Methods (~1500 words)

### 4.1 Data sources
- DepMap 24Q4: OmicsExpressionProteinCodingGenesTPMLogp1.csv (expression), CRISPRGeneEffect.csv (CRISPR), OmicsSomaticMutations.csv (mutations), Model.csv (metadata)
- PRISM 24Q2 Extended: secondary-screen-replicate-treatment-info.csv (compound info), PRISM_24Q2_extended.csv (sensitivity)
- TCGA-PAAD: expression + clinical from GDC (used for survival analysis in early phases only; negative result, HR=1.020)
- All data downloaded from DepMap portal (https://depmap.org/portal/)

### 4.2 Cell line stratification
- Expression matrix: 1,103 cell lines x 18,450 protein-coding genes (log2 TPM+1)
- Target gene: SLC2A8 (GLUT8)
- Target group: SLC2A8 expression above pan-cancer median AND RAS-pathway mutation (KRAS, NRAS, HRAS, BRAF, or PIK3CA hotspot mutations)
- Background: all remaining cell lines
- Split yields 223 target and 880 background lines

### 4.3 Unbiased drug sensitivity screen
- PRISM dataset: 6,790 compounds x ~480 cell lines with matched expression data
- For each compound, Welch's t-test comparing sensitivity (log2 FC viability) between target and background groups
- Multiple testing correction: Benjamini-Hochberg FDR
- Sensitivity convention: more negative = more killing (confirmed with bortezomib control)
- MOA enrichment: Fisher's exact test comparing proportion of significant hits (p<0.05) within each MOA class vs overall rate

### 4.4 Lineage-controlled regression
- Lineage encoded as one-hot dummy variables for top 15 OncotreeLineage categories
- Ordinary least squares regression: drug sensitivity ~ lineage dummies (residualized drug response)
- Ordinary least squares regression: SLC2A8 expression ~ lineage dummies (residualized expression)
- Partial correlation = Spearman correlation of residuals
- Significance by parametric p-value from partial correlation test

### 4.5 Co-expression analysis
- Spearman rank correlation between SLC2A8 and each candidate gene across all 1,103 cell lines
- No multiple testing correction applied to individual gene tests (hypothesis-driven candidate approach)
- Cohen's d for group comparisons (target vs background) with Welch's t-test p-values

### 4.6 CRISPR dependency analysis
- CRISPRGeneEffect.csv: Chronos-normalized gene effect scores
- More negative = more essential (growth-inhibiting when knocked out)
- Cohen's d and Welch's t-test comparing target vs background for candidate genes
- No genome-wide screen; candidate gene approach focused on metabolic and signaling genes identified by prior analyses

### 4.7 Convergence testing
- Cell-line-level Spearman correlation between CRISPR dependency scores and drug
  sensitivity values, tested pan-cancer, within target group, and within background
- Interaction regression: drug_sensitivity ~ SLC2A8 + gene_dep + SLC2A8*gene_dep
  (standardized predictors)
- Quadrant analysis: Fisher's exact test for target enrichment in dual-vulnerability
  (AKT-sensitive AND SCD-dependent) vs single-vulnerability quadrants
- Lineage-controlled partial correlations between SLC2A8 expression and lipid gene
  CRISPR dependencies

### 4.8 Statistical considerations
- All p-values are two-sided unless otherwise stated
- Effect sizes reported as Spearman rho (correlations), Cohen's d (group comparisons), or odds ratio (enrichment)
- Sample sizes range from 457-1,103 depending on data availability per analysis
- All analyses performed in Python 3.12 using scipy.stats, pandas, and numpy

### 4.9 Code availability
- All analysis scripts available at [GitHub repository URL]
- Key scripts: h13c_deep_dive.py, h13d_pi3k_validation.py, h13e_sanity_checks.py, h13f_metabolic_clue.py, h13g_convergence.py, h13_paper_figures.py

---

## 5. Figures

| Figure | Title | Type | Content |
|--------|-------|------|---------|
| 1 | SLC2A8 expression landscape and AKT inhibitor correlation | Two-panel | A: Violin plot of SLC2A8 by lineage (PDAC highlighted). B: Scatter of SLC2A8 vs afuresertib sensitivity |
| 2 | Unbiased PRISM drug screen identifies AKT inhibitors | Bar chart | MOA enrichment odds ratios with p-values (AKT, PI3K*, MEK, mTOR) |
| 3 | Lineage decomposition separates AKT signal from PI3K confound | Two-panel | A: Raw vs lineage-corrected SLC2A8-drug correlations. B: PDAC lineage effect per drug |
| 4 | SLC2A8 metabolic co-expression fingerprint | Horizontal bars | Spearman rho by gene, grouped by metabolic category (endosomal, lipid, glycolysis, Warburg) |
| 5 | CRISPR genetic dependencies in SLC2A8-high/RAS-mut cells | Box plots | 2x3 grid: SCD, PIK3C3, AKT1, FASN + controls (SLC2A1, MCOLN1) |
| 6 | SLC2A8 marks a metabolic state sensitive to AKT inhibition | Diagram | SLC2A8 marks subtype with co-expressed programs (endosomal, lipogenic, glycolytic); AKT sensitivity is the actionable arm; lipogenesis deps are independent subtype features |

| Supp Figure | Title | Content |
|-------------|-------|---------|
| S1 | PIK3C3 drives autophagy initiation signal | Per-gene decomposition showing PIK3C3 dominance |
| S2 | PRISM sanity check | Bortezomib/doxorubicin distributions confirming sign convention |
| S3 | Convergence testing: AKT sensitivity and lipogenesis dependency are independently distributed | A: SCD dep vs afuresertib sensitivity scatter (within target; rho=-0.030, p=0.72). B: Quadrant analysis (target enrichment in dual-vuln OR=2.32 but depletion from SCD-only OR=0.52). C: Interaction model coefficients (SLC2A8 sig, SCD null, interaction null) |

---

## 6. Supplementary Tables

| Table | Content |
|-------|---------|
| S1 | Full PRISM screen results: all 6,790 compounds with t-statistic, p-value, FDR |
| S2 | Complete MOA enrichment: all MOA classes with Fisher test results |
| S3 | SLC2A8 co-expression: genome-wide Spearman correlations (top 200 genes) |
| S4 | CRISPR dependencies: all tested genes with Cohen's d and p-values |
| S5 | Lineage correction: per-drug raw and partial correlations with SLC2A8 |
| S6 | Cell line classification: target vs background assignment for all 1,103 lines |
| S7 | Convergence analysis: cell-line-level correlations, interaction models, quadrant enrichment |

---

## 7. Narrative Arc

The paper follows a discovery narrative (hypothesis -> falsification -> pivot -> validation -> self-challenge):

1. **Hook** (Introduction): AKT inhibitors are clinically approved but lack biomarkers. KRAS-mutant cancers need new approaches.

2. **Entry point** (Results 2.1): SLC2A8 expression defines a coherent subgroup of RAS-mutant cancers. Initially studied for autophagy/trehalose vulnerability (cite briefly, don't dwell).

3. **Unbiased discovery** (Results 2.2): Drug screening reveals AKT inhibitors, not autophagy modulators. This is the "twist" -- the original hypothesis was wrong, but the subgroup stratification revealed something better.

4. **Critical decomposition** (Results 2.3): The PI3K signal was a lineage confound. Only AKT survives correction. This section demonstrates rigor and builds trust.

5. **Subtype characterization** (Results 2.4-2.5): Co-expression and CRISPR data describe the metabolic biology of the SLC2A8-high subtype. Framed as characterization, not as a targetable chain.

6. **Biomarker performance** (Results 2.6): SLC2A8 outperforms pathway scores and composite markers. Simple, measurable, translatable.

7. **Self-challenge** (Results 2.7): We tested whether the metabolic dependencies and drug sensitivity form a unified pathway. They don't. This narrows the claim to what the data supports and demonstrates scientific integrity.

8. **Clinical bridge** (Discussion): Capivasertib is already approved. SLC2A8 IHC is feasible. The claim is focused: SLC2A8 predicts AKT drug sensitivity. Experimental roadmap is clear.

---

## 8. Writing Notes

### Tone
- Measured, evidence-driven. Not hype.
- Acknowledge limitations prominently (computational only, effect sizes modest, correlational)
- Frame as hypothesis-generating, not confirmatory
- Emphasize the lineage correction as methodological contribution
- Present the H13g negative result as a strength: active self-falsification

### Key phrases to use
- "biomarker hypothesis" not "biomarker"
- "computational evidence" not "proof"
- "suggests" and "is consistent with" not "demonstrates" or "proves"
- "lineage-independent" -- the key differentiator from prior biomarker studies
- "metabolic subtype" -- describes what SLC2A8 marks
- "independently distributed features" -- for the lipogenesis/AKT non-convergence
- "characterizes the subtype" -- for co-expression and CRISPR data (not "validates the pathway")

### Key phrases to avoid
- "novel" (overused in cancer papers)
- "druggable" (commercial connotation)
- "precision medicine" (too broad)
- "paradigm" (no)
- "metabolic collapse" (overclaims mechanism; H13g refutes this)
- "synthetic lethality" (no cell-level convergence to support this)
- "targetable pathway" for the lipogenesis arm (only AKT drug sensitivity is the actionable claim)

### Word count targets
- Abstract: ~250 words
- Introduction: ~800 words
- Results: ~2800 words (added Section 2.7)
- Discussion: ~1200 words
- Methods: ~1500 words
- **Total main text: ~6550 words**

---

*Outline updated 2026-03-12 with H13g convergence analysis. Narrative refined
from "AKT-dependent lipogenic subtype" to cleaner "metabolic subtype selectively
sensitive to AKT inhibition." Scientific arc now complete: discovery -> validation
-> self-challenge -> refined claim. Ready for manuscript drafting.*
