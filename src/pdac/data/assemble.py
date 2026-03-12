"""Assemble raw GDC downloads into analysis-ready parquet files.

This script:
1. Resolves GDC file UUIDs to TCGA patient barcodes via the GDC API
2. Merges per-sample STAR count files into a single expression matrix
3. Loads and cleans clinical data
4. Computes TMB from per-sample MAF.gz files
5. Prepares GSE79668 replication cohort
6. Writes clean parquet files to data/processed/

Usage:
    python -m pdac.data.assemble --data-root data/
"""

import argparse
import json
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from pdac.data.prepare_methylation import filter_methylation_probes, load_gdc_methylation
from pdac.data.prepare_tcga import (
    compute_tmb_from_maf_dir,
    load_gdc_clinical,
)
from pdac.shared.logging import get_logger

logger = get_logger(__name__)

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"


def resolve_uuids_to_barcodes(file_uuids: list[str]) -> dict[str, str]:
    """Query GDC API to map file UUIDs to TCGA patient barcodes.

    Args:
        file_uuids: List of GDC file UUIDs (directory names from cart download).

    Returns:
        Dict mapping file_uuid -> TCGA patient barcode (12 chars, e.g. TCGA-IB-7897).
    """
    uuid_to_barcode: dict[str, str] = {}

    # GDC API accepts batches of up to 100 filters via POST
    batch_size = 50
    for i in range(0, len(file_uuids), batch_size):
        batch = file_uuids[i : i + batch_size]
        payload = {
            "filters": {
                "op": "in",
                "content": {
                    "field": "files.file_id",
                    "value": batch,
                },
            },
            "fields": "file_id,cases.submitter_id",
            "size": len(batch),
            "format": "JSON",
        }

        req = Request(
            GDC_FILES_ENDPOINT,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )

        with urlopen(req) as resp:  # noqa: S310
            result = json.loads(resp.read().decode("utf-8"))

        for hit in result.get("data", {}).get("hits", []):
            file_id = hit.get("file_id", "")
            cases = hit.get("cases", [])
            if cases:
                barcode = cases[0].get("submitter_id", "")
                if barcode:
                    uuid_to_barcode[file_id] = barcode

        logger.info(
            "GDC API batch resolved",
            extra={"batch": i // batch_size + 1, "resolved": len(uuid_to_barcode)},
        )

    return uuid_to_barcode


def merge_star_counts_with_barcodes(
    download_dir: Path,
    uuid_to_barcode: dict[str, str],
    value_column: str = "tpm_unstranded",
    gene_type_filter: str = "protein_coding",
) -> pd.DataFrame:
    """Merge per-sample STAR count files using resolved TCGA barcodes as sample IDs.

    Args:
        download_dir: Path to GDC gene expression download directory.
        uuid_to_barcode: Mapping from directory UUID to TCGA patient barcode.
        value_column: Which expression value to extract (tpm_unstranded recommended).
        gene_type_filter: Gene biotype filter (protein_coding recommended).

    Returns:
        DataFrame with rows = TCGA barcodes, columns = gene symbols.
    """
    tsv_files = sorted(download_dir.rglob("*.rna_seq.augmented_star_gene_counts.tsv"))
    if not tsv_files:
        msg = f"No STAR count files found in {download_dir}"
        raise FileNotFoundError(msg)

    logger.info("Merging STAR count files", extra={"n_files": len(tsv_files)})

    # Read first file to build gene reference
    ref = pd.read_csv(tsv_files[0], sep="\t", comment="#")
    ref = ref[ref["gene_id"].str.startswith("ENSG")]
    if gene_type_filter:
        ref = ref[ref["gene_type"] == gene_type_filter]

    gene_names = ref["gene_name"].values

    # Handle duplicate gene names
    seen: dict[str, int] = {}
    unique_names: list[str] = []
    for name in gene_names:
        if name in seen:
            seen[name] += 1
            unique_names.append(f"{name}_dup{seen[name]}")
        else:
            seen[name] = 0
            unique_names.append(name)

    n_genes = len(unique_names)

    # Process each file
    sample_data: dict[str, np.ndarray] = {}
    skipped = 0
    for tsv_path in tsv_files:
        # Parent directory name is the file UUID
        file_uuid = tsv_path.parent.name
        barcode = uuid_to_barcode.get(file_uuid, file_uuid)

        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        df = df[df["gene_id"].str.startswith("ENSG")]
        if gene_type_filter:
            df = df[df["gene_type"] == gene_type_filter]

        values = pd.to_numeric(df[value_column], errors="coerce").values
        if len(values) == n_genes:
            # If multiple aliquots per patient, keep first
            if barcode not in sample_data:
                sample_data[barcode] = values
        else:
            skipped += 1

    expression = pd.DataFrame(sample_data, index=unique_names).T
    expression.index.name = "patient_id"

    logger.info(
        "Expression matrix assembled",
        extra={
            "samples": expression.shape[0],
            "genes": expression.shape[1],
            "skipped": skipped,
        },
    )
    return expression


def prepare_gse79668(
    counts_path: Path,
    series_matrix_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare the GSE79668 replication cohort expression and clinical data.

    Args:
        counts_path: Path to GSE79668_51_tumors_sharedgenecounts.txt.
        series_matrix_path: Path to GSE79668_series_matrix.txt.

    Returns:
        Tuple of (expression_df, clinical_df).
        expression_df: rows = samples, columns = gene symbols.
        clinical_df: standardized clinical table.
    """
    # --- Expression ---
    expr = pd.read_csv(counts_path, sep="\t", index_col=0)
    # Gene IDs are like "CTU2_ENSG00000174177.7" — extract gene name before _ENSG
    gene_names = []
    for idx in expr.index:
        name = str(idx)
        # Split on _ENSG to get gene name
        parts = name.rsplit("_ENSG", maxsplit=1)
        gene_names.append(parts[0])

    expr.index = gene_names

    # Handle duplicate gene names
    seen: dict[str, int] = {}
    unique_names: list[str] = []
    for name in gene_names:
        if name in seen:
            seen[name] += 1
            unique_names.append(f"{name}_dup{seen[name]}")
        else:
            seen[name] = 0
            unique_names.append(name)
    expr.index = unique_names

    # Transpose so rows = samples
    expr = expr.T
    expr.index.name = "sample_id"
    logger.info(
        "GSE79668 expression loaded",
        extra={"samples": expr.shape[0], "genes": expr.shape[1]},
    )

    # --- Clinical from series matrix ---
    clinical = _parse_gse79668_clinical(series_matrix_path, list(expr.index))

    return expr, clinical


def _parse_gse79668_clinical(
    series_matrix_path: Path,
    sample_ids: list[str],
) -> pd.DataFrame:
    """Parse clinical data from GSE79668 series matrix header.

    The series matrix has !Sample_characteristics_ch1 rows with structured
    values like "patient survival status: Dead", "survival time (days): 393", etc.
    """
    characteristics: dict[str, list[str]] = {}
    sample_titles: list[str] = []

    with series_matrix_path.open() as f:
        for line in f:
            line = line.strip()
            if line.startswith("!Sample_title"):
                parts = line.split("\t")[1:]
                sample_titles = [p.strip('"') for p in parts]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")[1:]
                values = [p.strip('"') for p in parts]
                if values:
                    # Determine the key from the first value (e.g., "patient survival status: Dead")
                    key_part = values[0].split(":")[0].strip() if ":" in values[0] else values[0]
                    extracted = []
                    for v in values:
                        if ":" in v:
                            extracted.append(v.split(":", 1)[1].strip())
                        else:
                            extracted.append(v)
                    characteristics[key_part] = extracted

    # Map sample titles to the column names in the expression matrix
    # sample_titles are like "T_07_07_A082a_SL11228"
    # expression columns are like "T_07_07_A082a"
    # Strip the _SLxxxxx suffix to match
    title_to_expr_id: dict[str, str] = {}
    for title in sample_titles:
        # Remove _SLxxxxx suffix
        base = "_".join(title.split("_")[:-1]) if "_SL" in title else title
        # Find matching expression sample ID
        for sid in sample_ids:
            if sid == base or sid.startswith(base):
                title_to_expr_id[title] = sid
                break
        if title not in title_to_expr_id:
            title_to_expr_id[title] = base

    n = len(sample_titles)
    clinical = pd.DataFrame({"patient_id": [title_to_expr_id.get(t, t) for t in sample_titles]})

    # Extract known characteristics
    if "patient survival status" in characteristics:
        vals = characteristics["patient survival status"]
        clinical["os_event"] = [1.0 if v.lower() == "dead" else 0.0 for v in vals]

    if "survival time (days)" in characteristics:
        clinical["os_days"] = pd.to_numeric(characteristics["survival time (days)"], errors="coerce")

    if "age" in characteristics:
        clinical["age"] = pd.to_numeric(characteristics["age"], errors="coerce")

    # TNM staging -> approximate AJCC stage
    if "t (tnm score)" in characteristics:
        t_scores = characteristics["t (tnm score)"]
        n_scores = characteristics.get("n (tnm score)", [""] * n)
        m_scores = characteristics.get("m (tnm score)", [""] * n)
        stages = []
        stage_nums = []
        for t, ns, m in zip(t_scores, n_scores, m_scores):
            stage = _tnm_to_stage(t, ns, m)
            stages.append(stage)
            stage_nums.append(_stage_str_to_numeric(stage))
        clinical["ajcc_stage"] = stages
        clinical["stage_numeric"] = stage_nums

    pre = len(clinical)
    clinical = clinical.dropna(subset=["os_days", "os_event"])
    clinical = clinical[clinical["os_days"] > 0]

    logger.info(
        "GSE79668 clinical prepared",
        extra={"kept": len(clinical), "dropped": pre - len(clinical)},
    )
    return clinical.reset_index(drop=True)


def _tnm_to_stage(t: str, n: str, m: str) -> str:
    """Approximate AJCC stage from TNM components for PDAC."""
    if "M1" in m.upper():
        return "Stage IV"
    t_val = t.upper().replace("T", "")
    n_val = n.upper().replace("N", "")
    if t_val in ("1", "2") and n_val == "0":
        return "Stage I" if t_val == "1" else "Stage IB"
    if t_val == "3" and n_val == "0":
        return "Stage IIA"
    if n_val == "1":
        return "Stage IIB"
    if t_val == "4":
        return "Stage III"
    return "Stage II"


def _stage_str_to_numeric(stage: str) -> float | None:
    """Convert stage string to numeric."""
    if "IV" in stage:
        return 4.0
    if "III" in stage:
        return 3.0
    if "II" in stage:
        return 2.0
    if "I" in stage:
        return 1.0
    return None


def assemble_all(data_root: Path) -> None:
    """Run the full data assembly pipeline.

    Args:
        data_root: Root data directory (contains raw/ and processed/).
    """
    raw = data_root / "raw"
    processed = data_root / "processed"
    processed.mkdir(parents=True, exist_ok=True)

    # --- Locate raw directories ---
    tcga_dir = raw / "tcga_paad"
    geo_dir = raw / "gse79668"

    expr_dir = _find_subdir(tcga_dir, "gene_gdc_download")
    clinical_path = _find_file(tcga_dir, "clinical.tsv", recursive=True)
    maf_dir = _find_subdir(tcga_dir, "mutations_gdc_download")

    # --- 1. Resolve UUIDs to barcodes ---
    logger.info("Step 1: Resolving GDC file UUIDs to TCGA barcodes")
    barcode_cache = processed / "uuid_to_barcode.json"

    if barcode_cache.exists():
        logger.info("Loading cached UUID-to-barcode mapping")
        with barcode_cache.open() as f:
            uuid_to_barcode = json.load(f)
    else:
        # Get all directory UUIDs from the expression download
        file_uuids = [
            d.name for d in expr_dir.iterdir()
            if d.is_dir() and d.name != "__MACOSX"
        ]
        uuid_to_barcode = resolve_uuids_to_barcodes(file_uuids)

        # Cache the mapping
        with barcode_cache.open("w") as f:
            json.dump(uuid_to_barcode, f, indent=2)
        logger.info("UUID-to-barcode mapping cached", extra={"n_mapped": len(uuid_to_barcode)})

    # --- 2. Merge expression ---
    logger.info("Step 2: Merging expression matrix")
    expression = merge_star_counts_with_barcodes(
        expr_dir, uuid_to_barcode, value_column="tpm_unstranded"
    )
    expr_out = processed / "tcga_paad_expression_tpm.parquet"
    expression.to_parquet(expr_out)
    logger.info("Expression saved", extra={"path": str(expr_out), "shape": list(expression.shape)})

    # --- 3. Clinical ---
    logger.info("Step 3: Loading clinical data")
    clinical = load_gdc_clinical(clinical_path)
    clinical_out = processed / "clinical" / "tcga_paad_clinical.parquet"
    clinical_out.parent.mkdir(exist_ok=True)
    clinical.to_parquet(clinical_out)
    logger.info("Clinical saved", extra={"path": str(clinical_out), "patients": len(clinical)})

    # --- 4. TMB (optional) ---
    if maf_dir and maf_dir.exists():
        logger.info("Step 4: Computing TMB from MAF files")
        tmb = compute_tmb_from_maf_dir(maf_dir)
        tmb_out = processed / "clinical" / "tcga_paad_tmb.parquet"
        tmb.to_frame().to_parquet(tmb_out)
        logger.info("TMB saved", extra={"path": str(tmb_out), "patients": len(tmb)})

        # Merge TMB into clinical
        clinical = clinical.merge(
            tmb.reset_index().rename(columns={"index": "patient_id"}),
            on="patient_id",
            how="left",
        )
        clinical.to_parquet(clinical_out)
        logger.info("Clinical updated with TMB", extra={"patients_with_tmb": clinical["tmb_log"].notna().sum()})

    # --- 5. Methylation (H03) ---
    meth_dir = _find_subdir_optional(tcga_dir, "methylation_gdc_download")
    if meth_dir and meth_dir.exists():
        logger.info("Step 5: Assembling 450K methylation data")

        # Methylation files have different UUIDs than expression — resolve them
        meth_uuids = [
            d.name for d in meth_dir.iterdir()
            if d.is_dir() and d.name != "__MACOSX" and not d.name.startswith(".")
        ]
        new_uuids = [u for u in meth_uuids if u not in uuid_to_barcode]
        if new_uuids:
            logger.info(
                "Resolving methylation UUIDs to barcodes",
                extra={"n_new": len(new_uuids)},
            )
            meth_barcode_map = resolve_uuids_to_barcodes(new_uuids)
            uuid_to_barcode.update(meth_barcode_map)
            # Update cache
            with barcode_cache.open("w") as f:
                json.dump(uuid_to_barcode, f, indent=2)
            logger.info(
                "Barcode cache updated with methylation UUIDs",
                extra={"n_total": len(uuid_to_barcode)},
            )

        beta_matrix = load_gdc_methylation(meth_dir, uuid_to_barcode)
        beta_filtered = filter_methylation_probes(beta_matrix)
        meth_out = processed / "tcga_paad_methylation_450k.parquet"
        beta_filtered.to_parquet(meth_out)
        logger.info(
            "Methylation saved",
            extra={"path": str(meth_out), "shape": list(beta_filtered.shape)},
        )
    else:
        logger.info("Methylation download dir not found; skipping (see data/README.md)")

    # --- 6. GSE79668 replication ---
    counts_path = geo_dir / "GSE79668_51_tumors_sharedgenecounts.txt"
    matrix_path = geo_dir / "GSE79668_series_matrix.txt"

    if counts_path.exists() and matrix_path.exists():
        logger.info("Step 5: Preparing GSE79668 replication cohort")
        geo_expr, geo_clinical = prepare_gse79668(counts_path, matrix_path)

        geo_expr_out = processed / "gse79668_expression_counts.parquet"
        geo_clinical_out = processed / "clinical" / "gse79668_clinical.parquet"

        geo_expr.to_parquet(geo_expr_out)
        geo_clinical.to_parquet(geo_clinical_out)

        logger.info(
            "GSE79668 saved",
            extra={
                "expr_shape": list(geo_expr.shape),
                "clinical_patients": len(geo_clinical),
            },
        )
    else:
        logger.info("GSE79668 files not found; skipping replication cohort")

    # --- Summary ---
    logger.info("=" * 60)
    logger.info("DATA ASSEMBLY COMPLETE")
    logger.info("=" * 60)
    logger.info(
        "TCGA-PAAD",
        extra={
            "expression_samples": expression.shape[0],
            "expression_genes": expression.shape[1],
            "clinical_patients": len(clinical),
        },
    )
    if counts_path.exists():
        logger.info(
            "GSE79668",
            extra={
                "expression_samples": geo_expr.shape[0],
                "clinical_patients": len(geo_clinical),
            },
        )
    logger.info("=" * 60)


def _find_subdir_optional(parent: Path, prefix: str) -> Path | None:
    """Find a subdirectory starting with prefix, or None if not found."""
    if not parent.exists():
        return None
    for d in sorted(parent.iterdir()):
        if d.is_dir() and d.name.startswith(prefix):
            return d
    return None


def _find_subdir(parent: Path, prefix: str) -> Path:
    """Find a subdirectory starting with prefix."""
    for d in sorted(parent.iterdir()):
        if d.is_dir() and d.name.startswith(prefix):
            return d
    msg = f"No directory starting with '{prefix}' found in {parent}"
    raise FileNotFoundError(msg)


def _find_file(parent: Path, name: str, recursive: bool = False) -> Path:
    """Find a file by name, optionally recursively."""
    if recursive:
        matches = list(parent.rglob(name))
    else:
        matches = list(parent.glob(name))
    if not matches:
        msg = f"File '{name}' not found in {parent}"
        raise FileNotFoundError(msg)
    return matches[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble PDAC research data")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data"),
        help="Root data directory (contains raw/ and processed/)",
    )
    args = parser.parse_args()

    assemble_all(args.data_root)
