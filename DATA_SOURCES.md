# Data Sources

All analyses use publicly available datasets. No data files are included in this repository.
Download the following datasets and place them in `data/raw/depmap/` before running scripts.

## Required Datasets

### DepMap (Release 24Q4)

Download from [depmap.org/portal](https://depmap.org/portal/):

| File | Description |
|------|-------------|
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | Gene expression (log2 TPM+1, 18,450 genes) |
| `CRISPRGeneEffect.csv` | CRISPR gene effect scores (Chronos-normalized) |
| `OmicsSomaticMutations.csv` | Somatic mutation calls |
| `Model.csv` | Cell line metadata (lineage, disease, etc.) |

### PRISM Repurposing (Release 24Q2)

Download from [depmap.org/repurposing](https://depmap.org/repurposing/):

| File | Description |
|------|-------------|
| `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv` | Drug sensitivity (6,790 compounds) |
| `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv` | Compound annotations (MOA, target) |

### TCGA-PAAD (Optional — for TCGA validation in h13h)

Download from [portal.gdc.cancer.gov](https://portal.gdc.cancer.gov/):

- TCGA-PAAD STAR-Counts gene expression (n=178 tumors)
- Place processed expression matrix in `data/processed/`

## Directory Structure After Download

```
data/
  raw/
    depmap/
      OmicsExpressionProteinCodingGenesTPMLogp1.csv
      CRISPRGeneEffect.csv
      OmicsSomaticMutations.csv
      Model.csv
      Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv
      Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv
  processed/          # Created by scripts
```

## Data Access

All datasets are open access. DepMap and PRISM require no application.
TCGA-PAAD is open-access through the GDC Data Portal (free account may be required for bulk downloads).

## References

- Tsherniak, A. et al. Defining a Cancer Dependency Map. *Cell* 170, 564-576 (2017).
- Corsello, S.M. et al. Discovering the anti-cancer potential of non-oncology drugs. *Nat Cancer* 1, 235-248 (2020).
- Dempster, J.M. et al. Chronos: a CRISPR cell population dynamics model. *bioRxiv* (2021).
