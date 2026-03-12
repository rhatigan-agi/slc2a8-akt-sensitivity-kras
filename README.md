# SLC2A8 Identifies a Metabolic Tumor Subtype Sensitive to AKT Inhibitors

Code and data for the preprint:

> **SLC2A8 expression identifies a metabolic tumor subtype selectively sensitive to AKT inhibitors across KRAS-mutant cancers**
>
> Jeffrey Rhatigan (2026). *bioRxiv* [preprint]. DOI: [pending]

## Key Finding

SLC2A8 (GLUT8) expression marks a metabolic subtype within RAS-mutant cancers — characterized by endosomal trafficking, lipogenesis, and oxidative glycolysis programs — that is selectively sensitive to AKT inhibitors (capivasertib, ipatasertib, afuresertib). The signal is:

- **Lineage-independent** (survives correction across 15 tissue types)
- **RAS-context-dependent** (interaction p=0.004--0.040; null without RAS mutation)
- **Validated across four data layers** (PRISM drug screening, lineage-corrected pharmacology, co-expression, CRISPR dependency)
- **Self-challenged** (convergence testing shows AKT sensitivity and lipogenesis dependency are independently distributed)

## Discovery Narrative

This project began with a different hypothesis: that SLC2A8-high, RAS-mutant cells would exhibit autophagy-related vulnerability based on GLUT8's role in trehalose metabolism. Systematic testing falsified that hypothesis and revealed a stronger signal — selective sensitivity to AKT inhibitors. The analysis scripts preserve this full discovery arc, including negative results (H13g convergence testing), because the falsification process is part of the scientific contribution.

## Repository Structure

```
paper/
  main.tex                  # LaTeX manuscript
  main.pdf                  # Compiled preprint
  references.bib            # Bibliography
  generate_figures.py        # Regenerate all figures from source data
  figures/                   # Publication figures (6 main + 3 supplementary)

scripts/
  h13c_deep_dive.py          # Unbiased PRISM screen + MOA enrichment
  h13d_pi3k_validation.py    # PI3K/AKT pathway validation
  h13e_sanity_checks.py      # PRISM orientation + lineage correction
  h13f_metabolic_clue.py     # Metabolic co-expression + CRISPR dependencies
  h13g_convergence.py        # Convergence testing (self-falsification)
  h13h_reviewer_analyses.py  # RAS interaction model + ROC + TCGA descriptive

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
         sensitive to AKT inhibitors across KRAS-mutant cancers},
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
