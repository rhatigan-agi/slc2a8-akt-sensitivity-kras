# SLC2A8 Identifies a Metabolic Tumor Subtype Sensitive to AKT Inhibitors across RAS-Pathway-Mutant Cancers

Code and data for the preprint:

> **SLC2A8 expression identifies a metabolic tumor subtype selectively sensitive to AKT inhibitors across RAS-pathway-mutant cancers**
>
> Jeffrey Rhatigan (2026). *bioRxiv* [preprint]. DOI: [pending]

*Dedicated to Jennifer J. Rhatigan (1964--2025).
She met her illness with unwavering optimism and the conviction that everything---even her death---was meant to be.
Her passing redirected my focus to this research.
If these findings help even one family going through this battle, she was right: it was meant to be.*

---

## TL;DR — Why This Matters

**The problem:** AKT inhibitors (capivasertib, ipatasertib, afuresertib) work, but we don't know which patients to give them to beyond PIK3CA/PTEN mutations — and those mutations are rare in KRAS-mutant cancers like PDAC.

**What we found:** SLC2A8 (GLUT8) expression, in RAS-mutant contexts, stratifies tumors with differential sensitivity to AKT inhibitors across cancer types. Notably, in stratifier-free analysis (no mutation grouping), AKT inhibitors rank 1--5 out of 6,790 compounds by correlation with SLC2A8 expression alone. In PDAC:
- **~69% of cell lines** are SLC2A8-high (the target population)
- **~95% of PDAC** is KRAS-mutant (the RAS-pathway context where SLC2A8 stratifies)
- **<10% of PDAC** has PIK3CA/PTEN mutations (what current tests catch)

**The clinical gap:** As of March 2026, SLC2A8 is **not on any major clinical genomic panel** — not FoundationOne CDx (324 genes), not Tempus xT, not Caris MI Profile.

**Clinical context.** AKT inhibitors such as capivasertib are clinically available (FDA-approved 2023 for breast cancer), but robust biomarkers for patient selection in RAS-pathway-mutant cancers remain limited.

### Validation Depth (Not "Just Another AI Hypothesis")

| Layer | What We Tested | Result |
|-------|---------------|--------|
| Unbiased drug screen | 6,790 PRISM compounds | AKT inhibitors: OR=19.93, p=3.3e-9 |
| Lineage correction | Regressed out tissue type | AKT signal survives (r=-0.131 to -0.147, all p<0.005); PI3K signal collapses |
| PIK3CA/PTEN confound | Added mutation covariates (Davies 2012 challenge) | ~15% attenuation, all 3 drugs remain p<0.01 |
| Co-expression | Transcriptional fingerprint | Endosomal trafficking + lipogenesis programs |
| CRISPR dependency | Genetic knockouts | SCD d=-0.316, AKT1 d=-0.277 |
| Self-falsification | Tested causal chain hypothesis | Rejected overclaim (H13g); lipogenesis deps are independent |
| RAS interaction | Interaction regression | SLC2A8 alone null; SLC2A8 x KRAS significant (p=0.004) |

### Population Impact

| Cancer Type | KRAS Mutation Rate | SLC2A8-High Prevalence | Current PIK3CA Testing Captures |
|---|---|---|---|
| PDAC | ~95% | ~69% of lines | ~5-7% of patients |
| Colorectal | ~45% | Present in subset | ~15-20% of patients |
| Non-small cell lung | ~30% | Variable by subtype | ~5-10% of patients |

*SLC2A8 identifies a large subset of tumors with increased AKT inhibitor sensitivity that would not be captured by PIK3CA/PTEN mutation testing.*

---

## Key Finding

SLC2A8 (GLUT8) expression marks a metabolic subtype within RAS-mutant cancers — characterized by endosomal trafficking, lipogenesis, and oxidative glycolysis programs — that is selectively sensitive to AKT inhibitors (capivasertib, ipatasertib, afuresertib). The signal is:

- **Lineage-independent** (survives correction across 15 tissue types)
- **RAS-context-dependent** (interaction p=0.004--0.040; null without RAS mutation)
- **PIK3CA/PTEN-independent** (~15% attenuation after mutation adjustment; all drugs remain p<0.01)
- **Validated across five data layers** (PRISM drug screening, lineage-corrected pharmacology, co-expression, CRISPR dependency, PIK3CA/PTEN confound control)
- **Self-challenged** (convergence testing shows AKT sensitivity and lipogenesis dependency are independently distributed)
- **Not currently captured by clinical panels** (SLC2A8 absent from FoundationOne CDx, Tempus xT, and other panels as of March 2026)

## Discovery Narrative

This project began with a different hypothesis: that SLC2A8-high, RAS-mutant cells would exhibit autophagy-related vulnerability based on GLUT8's role in trehalose metabolism. Systematic testing falsified that hypothesis and revealed a stronger signal — selective sensitivity to AKT inhibitors. The analysis scripts preserve this full discovery arc, including negative results (H13g convergence testing), because the falsification process is part of the scientific contribution.

## Repository Structure

```
paper/
  main.tex                  # LaTeX manuscript
  main.pdf                  # Compiled preprint
  references.bib            # Bibliography
  generate_figures.py        # Regenerate all figures from source data
  figures/                   # Publication figures (6 main + 7 supplementary)

scripts/
  h13c_deep_dive.py          # Unbiased PRISM screen + MOA enrichment
  h13d_pi3k_validation.py    # PI3K/AKT pathway validation
  h13e_sanity_checks.py      # PRISM orientation + lineage correction
  h13f_metabolic_clue.py     # Metabolic co-expression + CRISPR dependencies
  h13g_convergence.py        # Convergence testing (self-falsification)
  h13h_reviewer_analyses.py  # RAS interaction model + ROC + TCGA descriptive
  h13i_pik3ca_confound.py    # PIK3CA/PTEN mutation confound check (Davies 2012)
  h13j_kras_only_sensitivity.py  # KRAS-only stratifier sensitivity analysis
  h13k_slc2a8_independence.py    # Stratifier-free validation + subgroup analysis

src/pdac/
  h13_trehalose_vulnerability/  # Core analysis modules
    data.py                     # DepMap/PRISM data loading
    enrichment.py               # Dependency enrichment tests
    gene_sets.py                # Pre-specified gene set definitions
    scoring.py                  # Vulnerability scoring + stratification
    pipeline.py                 # Full H13 pipeline
    h13b_prism.py               # PRISM sensitivity analysis
    h13b_lineage.py             # Lineage correction utilities
    ...
  shared/
    logging.py                  # Structured logging
    cox.py                      # Cox regression utilities

tests/                          # Unit tests
results/                        # Analysis outputs + documentation
```

## Reproducing the Analysis

### Prerequisites

- Python >= 3.12
- Download source datasets (see [DATA_SOURCES.md](DATA_SOURCES.md))

### Setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### Quick Test

To verify the environment is working:

```bash
python scripts/h13c_deep_dive.py
```

This runs the unbiased PRISM screen (~2 min with data in place) and produces logged output.

### Running Analyses

Scripts are designed to run sequentially. Each produces logged output and writes results to `results/`.

```bash
# Core discovery pipeline
python scripts/h13c_deep_dive.py        # Unbiased PRISM screen
python scripts/h13d_pi3k_validation.py  # PI3K vs AKT decomposition
python scripts/h13e_sanity_checks.py    # Lineage correction
python scripts/h13f_metabolic_clue.py   # Metabolic fingerprint
python scripts/h13g_convergence.py      # Self-falsification
python scripts/h13h_reviewer_analyses.py # RAS interaction + TCGA
python scripts/h13i_pik3ca_confound.py  # PIK3CA/PTEN confound check

# Regenerate paper figures
python paper/generate_figures.py
```

### Compiling the Paper

```bash
cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## Data

No data is included in this repository. All analyses use publicly available datasets from:

- **DepMap** (24Q4): Gene expression, CRISPR dependencies, mutations, cell line metadata
- **PRISM** (24Q2): Drug sensitivity screening (6,790 compounds across ~480 cell lines)
- **TCGA-PAAD**: Patient tumor RNA-seq (n=178, optional validation)

See [DATA_SOURCES.md](DATA_SOURCES.md) for download instructions.

## Methods Summary

| Analysis | Script | Key Result |
|----------|--------|------------|
| Unbiased PRISM screen | `h13c` | AKT inhibitors most enriched MOA (OR=19.93, p=3.3e-9) |
| Lineage correction | `h13e` | AKT signal survives (partial r=-0.144); PI3K signal collapses |
| Co-expression | `h13f` | Endosomal trafficking + lipogenesis programs |
| CRISPR dependencies | `h13f` | SCD (d=-0.316), PIK3C3 (d=-0.302), AKT1 (d=-0.277) |
| Convergence testing | `h13g` | Drug sensitivity and CRISPR deps independently distributed |
| RAS interaction | `h13h` | Signal exists only in RAS-mutant context (interaction p=0.004) |
| PIK3CA/PTEN confound | `h13i` | Signal survives mutation adjustment (~15% attenuation; all p<0.01) |
| Stratifier sensitivity | `h13j` | Excluding PIK3CA collapses group MOA enrichment; MEK enrichment surges |
| Stratifier-free validation | `h13k` | AKT inhibitors rank 1–5 among 6,790 compounds by SLC2A8 correlation alone |

## Acknowledgments

This work used publicly available data from the DepMap and PRISM projects at the Broad Institute.
Parts of the computational workflow were assisted by a proprietary research system
(Rhatigan Labs) integrating large language models with custom analysis pipelines. All analyses,
code execution, and results validation were performed and verified by the author.

## License

MIT

## Citation

If you use this code or find the biomarker hypothesis useful, please cite:

```bibtex
@article{rhatigan2026slc2a8,
  title={SLC2A8 expression identifies a metabolic tumor subtype selectively
         sensitive to AKT inhibitors across RAS-pathway-mutant cancers},
  author={Rhatigan, Jeffrey},
  journal={bioRxiv},
  year={2026},
  note={bioRxiv preprint},
  doi={pending}
}
```

## About Rhatigan Labs

Rhatigan Labs develops AI-assisted research infrastructure combining large language models with structured data analysis pipelines.

This project represents an early demonstration of AI-assisted hypothesis discovery validated through reproducible computational analysis.
