"""Prepare GEO replication cohort (GSE79668) for H01 analysis.

GSE79668: 66 PDAC samples with bulk RNA-seq and overall survival data.
Used as independent replication for the TCGA-PAAD primary analysis.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from pdac.shared.logging import get_logger

logger = get_logger(__name__)


def load_geo_series_matrix(matrix_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load a GEO series matrix file, extracting expression and phenotype data.

    GEO series matrix files have a metadata header (lines starting with !)
    followed by the expression matrix. Phenotype/clinical data is embedded
    in the header lines starting with !Sample_.

    Args:
        matrix_path: Path to the series matrix .txt file.

    Returns:
        Tuple of (expression_df, phenotype_df).
        expression_df: genes x samples (will need transposing for our pipeline).
        phenotype_df: samples x phenotype fields.
    """
    pheno_lines: list[str] = []
    header_line: str | None = None
    data_start: int = 0

    with matrix_path.open() as f:
        for i, line in enumerate(f):
            if line.startswith("!Sample_"):
                pheno_lines.append(line.strip())
            elif line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                header_line = line.strip()
                data_start = i
                break

    if header_line is None:
        msg = f"Could not find expression data header in {matrix_path}"
        raise ValueError(msg)

    # Parse expression matrix
    expr_df = pd.read_csv(
        matrix_path, sep="\t", skiprows=data_start, index_col=0,
        comment="!", low_memory=False,
    )
    # Drop the trailing !series_matrix_table_end row if present
    expr_df = expr_df[~expr_df.index.astype(str).str.startswith("!")]

    # Parse phenotype data from header
    pheno_df = _parse_pheno_lines(pheno_lines)

    logger.info(
        "GEO series matrix loaded",
        extra={
            "expression_shape": list(expr_df.shape),
            "phenotype_fields": len(pheno_df.columns),
            "samples": len(pheno_df),
        },
    )
    return expr_df, pheno_df


def _parse_pheno_lines(lines: list[str]) -> pd.DataFrame:
    """Parse !Sample_ lines from series matrix into a phenotype DataFrame."""
    records: dict[str, list[str]] = {}
    for line in lines:
        parts = line.split("\t")
        key = parts[0].replace("!Sample_", "").strip('"')
        values = [v.strip('"') for v in parts[1:]]
        if key in records:
            # Multiple lines with same key -> concatenate
            records[key] = [f"{a}; {b}" for a, b in zip(records[key], values)]
        else:
            records[key] = values

    return pd.DataFrame(records)


def prepare_gse79668_clinical(
    phenotype_df: pd.DataFrame,
    survival_col_patterns: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Extract standardized clinical data from GSE79668 phenotype table.

    GEO phenotype data varies by dataset. This function searches for
    survival-related columns using pattern matching on column names
    and characteristic fields.

    Args:
        phenotype_df: Phenotype DataFrame from load_geo_series_matrix.
        survival_col_patterns: Optional dict mapping our schema fields
            to regex patterns for matching GEO column names.

    Returns:
        Standardized clinical DataFrame with patient_id, os_days, os_event.
    """
    if survival_col_patterns is None:
        survival_col_patterns = {
            "os_days": r"(?i)(overall.?survival|os).*(days|time|month)",
            "os_event": r"(?i)(vital|status|event|dead|alive|os.*(status|event))",
            "age": r"(?i)(age|years)",
        }

    clinical = pd.DataFrame()
    clinical["patient_id"] = (
        phenotype_df.get("geo_accession", phenotype_df.get("title", pd.Series(dtype=str)))
    )

    # Search characteristics columns for survival data
    char_cols = [c for c in phenotype_df.columns if "characteristics" in c.lower()]
    all_searchable = list(phenotype_df.columns) + char_cols

    for target, pattern in survival_col_patterns.items():
        found = _find_matching_column(phenotype_df, pattern, all_searchable)
        if found is not None:
            clinical[target] = _extract_numeric_from_field(phenotype_df[found])

    # If os_event is text-based, convert
    if "os_event" in clinical.columns:
        if clinical["os_event"].dtype == object:
            clinical["os_event"] = (
                clinical["os_event"].astype(str).str.lower()
                .map({"dead": 1, "alive": 0, "deceased": 1, "living": 0, "1": 1, "0": 0})
            )

    # If survival reported in months, convert to days
    if "os_days" in clinical.columns:
        vals = pd.to_numeric(clinical["os_days"], errors="coerce")
        if vals.median() < 100:  # likely months
            logger.info("Converting OS from months to days")
            clinical["os_days"] = (vals * 30.44).round(0)
        else:
            clinical["os_days"] = vals

    pre = len(clinical)
    clinical = clinical.dropna(subset=["os_days", "os_event"])
    clinical = clinical[clinical["os_days"] > 0]
    logger.info(
        "GSE79668 clinical prepared",
        extra={"kept": len(clinical), "dropped": pre - len(clinical)},
    )
    return clinical.reset_index(drop=True)


def _find_matching_column(
    df: pd.DataFrame, pattern: str, col_list: list[str]
) -> str | None:
    """Find first column matching a regex pattern."""
    import re

    for col in col_list:
        if col in df.columns and re.search(pattern, col):
            return col
        # Also search within characteristic values
        if col in df.columns:
            sample_vals = df[col].astype(str).head(3).str.cat(sep=" ")
            if re.search(pattern, sample_vals):
                return col
    return None


def _extract_numeric_from_field(series: pd.Series) -> pd.Series:
    """Extract numeric values from GEO characteristic fields.

    GEO stores values like "overall survival (months): 12.5" — we need
    to extract just the numeric part.
    """
    import re

    def _extract(val: str) -> float | None:
        if not isinstance(val, str):
            return None
        # Try to find a number after a colon
        match = re.search(r":\s*([\d.]+)", val)
        if match:
            return float(match.group(1))
        # Try plain numeric
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    return series.apply(_extract)
