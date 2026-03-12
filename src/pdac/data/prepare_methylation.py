"""Load and process GDC Illumina 450K methylation beta-value files.

Expects the GDC download directory structure where each subdirectory
is named by file UUID and contains a single methylation beta-value TSV.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from pdac.shared.logging import get_logger

logger = get_logger(__name__)

# Chen et al. 2013 cross-reactive probes (450K) — abbreviated list of prefixes
# In practice, a full list would be loaded from a resource file.
_CROSS_REACTIVE_PREFIXES = ("ch.", "rs")

# Chromosomes to exclude (sex chromosomes introduce sex-based bias)
_EXCLUDE_CHROMOSOMES = {"chrX", "chrY"}


def load_gdc_methylation(
    download_dir: Path,
    uuid_to_barcode: dict[str, str],
) -> pd.DataFrame:
    """Merge per-sample GDC 450K methylation beta-value files into a matrix.

    Each subdirectory in *download_dir* is named by its GDC file UUID and
    contains a single ``*.methylation_array.sesame.level3betas.tsv`` file
    (or similar naming from the Sesame pipeline).

    Args:
        download_dir: Path to GDC methylation download directory.
        uuid_to_barcode: Mapping from directory UUID to TCGA patient barcode.

    Returns:
        DataFrame with rows = TCGA barcodes, columns = CpG probe IDs,
        values = beta values (0–1).

    Raises:
        FileNotFoundError: If no methylation files found.
    """
    beta_files = sorted(download_dir.rglob("*.level3betas.txt"))
    if not beta_files:
        beta_files = sorted(download_dir.rglob("*.level3betas.tsv"))
    if not beta_files:
        beta_files = sorted(download_dir.rglob("*methylation*.txt"))
    if not beta_files:
        msg = f"No methylation beta-value files found in {download_dir}"
        raise FileNotFoundError(msg)

    logger.info("Loading methylation files", extra={"n_files": len(beta_files)})

    sample_frames: dict[str, pd.Series] = {}
    skipped = 0

    for fpath in beta_files:
        file_uuid = fpath.parent.name
        barcode = uuid_to_barcode.get(file_uuid, file_uuid)

        # Skip duplicate patients (keep first aliquot)
        if barcode in sample_frames:
            skipped += 1
            continue

        df = pd.read_csv(fpath, sep="\t", index_col=0, header=None, names=["probe_id", "beta"])
        sample_frames[barcode] = df["beta"]

    if not sample_frames:
        msg = "No valid methylation samples loaded"
        raise FileNotFoundError(msg)

    beta_matrix = pd.DataFrame(sample_frames).T
    beta_matrix.index.name = "patient_id"

    logger.info(
        "Methylation matrix assembled",
        extra={
            "samples": beta_matrix.shape[0],
            "probes": beta_matrix.shape[1],
            "skipped_duplicates": skipped,
        },
    )
    return beta_matrix


def filter_methylation_probes(
    beta_matrix: pd.DataFrame,
    max_na_fraction: float = 0.2,
    min_variance: float = 0.001,
) -> pd.DataFrame:
    """Filter methylation probes for quality and informativeness.

    Removes:
    1. Cross-reactive probes (Chen et al. 2013 prefix heuristic)
    2. Probes with high missingness (> max_na_fraction)
    3. Low-variance probes (< min_variance)

    Remaining NaN values are imputed with probe median.
    All heavy computation uses numpy arrays to avoid pandas column-wise overhead.

    Args:
        beta_matrix: Raw beta-value matrix (patients x probes).
        max_na_fraction: Maximum fraction of missing values per probe.
        min_variance: Minimum variance threshold for probe retention.

    Returns:
        Filtered beta-value matrix with imputed NaN values (probe median).
    """
    n_start = beta_matrix.shape[1]
    probes = beta_matrix.columns.to_numpy()

    # 1. Remove cross-reactive probes (fast string filter on probe names)
    keep_xr = np.array([
        not any(str(p).startswith(pref) for pref in _CROSS_REACTIVE_PREFIXES)
        for p in probes
    ])
    probes = probes[keep_xr]
    n_after_xreact = len(probes)
    logger.info("Cross-reactive probes removed", extra={"kept": n_after_xreact})

    # Work in numpy from here — orders of magnitude faster than pandas on wide matrices
    values = beta_matrix.values[:, keep_xr].astype(np.float64)
    sample_ids = beta_matrix.index

    # 2. Remove high-missingness probes
    na_frac = np.isnan(values).mean(axis=0)
    keep_na = na_frac <= max_na_fraction
    values = values[:, keep_na]
    probes = probes[keep_na]
    n_after_na = len(probes)
    logger.info("High-NA probes removed", extra={"kept": n_after_na})

    # 3. Impute remaining NaN with probe median (numpy column-wise)
    nan_mask = np.isnan(values)
    if nan_mask.any():
        col_medians = np.nanmedian(values, axis=0)
        row_idx, col_idx = np.where(nan_mask)
        values[row_idx, col_idx] = col_medians[col_idx]
    logger.info("NaN imputation complete")

    # 4. Remove low-variance probes
    variances = np.var(values, axis=0)
    keep_var = variances >= min_variance
    values = values[:, keep_var]
    probes = probes[keep_var]
    n_final = len(probes)

    logger.info(
        "Methylation probes filtered",
        extra={
            "n_start": n_start,
            "after_cross_reactive": n_after_xreact,
            "after_na_filter": n_after_na,
            "n_final": n_final,
            "removed_total": n_start - n_final,
        },
    )

    return pd.DataFrame(values, index=sample_ids, columns=probes)
