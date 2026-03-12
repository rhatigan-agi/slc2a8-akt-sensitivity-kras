"""Prepare TCGA-PAAD clinical and expression data for H01 analysis.

Reads raw GDC per-sample downloads and produces clean parquet files with
standardized column names ready for deconvolution and Cox modeling.

GDC downloads have this structure:
  gene_gdc_download_.../
    <uuid>/<hash>.rna_seq.augmented_star_gene_counts.tsv  (one per sample)
    MANIFEST.txt

  clinical_<date>/
    clinical.tsv  (wide TSV with dot-separated column names)

  mutations_gdc_download_.../
    <uuid>/<hash>.wxs.aliquot_ensemble_masked.maf.gz  (one per sample)
"""

import gzip
import json
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from pdac.shared.logging import get_logger

logger = get_logger(__name__)

# Standard column names used throughout the pipeline
CLINICAL_COLS = [
    "patient_id",
    "os_days",
    "os_event",
    "age",
    "stage_numeric",
    "ajcc_stage",
]

# GDC sentinel for missing values
GDC_MISSING = {"'--", "--", "not reported", "Not Reported", "nan", ""}


def merge_gdc_star_counts(
    download_dir: Path,
    value_column: str = "tpm_unstranded",
    gene_type_filter: str = "protein_coding",
) -> pd.DataFrame:
    """Merge per-sample STAR gene count files into a single expression matrix.

    Each GDC file has columns: gene_id, gene_name, gene_type, unstranded,
    stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded,
    fpkm_uq_unstranded.

    We extract one value column per sample, using gene_name as the gene index.

    Args:
        download_dir: Path to the GDC gene expression download directory
            containing UUID subdirectories with .rna_seq.augmented_star_gene_counts.tsv files.
        value_column: Which count column to extract. Options: 'unstranded' (raw counts),
            'tpm_unstranded' (TPM), 'fpkm_unstranded' (FPKM), 'fpkm_uq_unstranded' (FPKM-UQ).
            TPM recommended for deconvolution.
        gene_type_filter: Only keep genes of this type (e.g. 'protein_coding').
            Set to None to keep all.

    Returns:
        DataFrame with rows = samples (TCGA barcodes), columns = gene symbols.
    """
    tsv_files = sorted(download_dir.rglob("*.rna_seq.augmented_star_gene_counts.tsv"))
    if not tsv_files:
        msg = f"No STAR count files found in {download_dir}"
        raise FileNotFoundError(msg)

    logger.info("Merging STAR count files", extra={"n_files": len(tsv_files)})

    # Read first file to get gene list
    ref = pd.read_csv(tsv_files[0], sep="\t", comment="#")
    # Skip the N_unmapped, N_multimapping, N_noFeature, N_ambiguous rows
    ref = ref[ref["gene_id"].str.startswith("ENSG")]
    if gene_type_filter:
        ref = ref[ref["gene_type"] == gene_type_filter]
    gene_names = ref["gene_name"].values
    gene_ids = ref["gene_id"].values

    # Handle duplicate gene names by appending Ensembl ID suffix
    seen: dict[str, int] = {}
    unique_names: list[str] = []
    for name in gene_names:
        if name in seen:
            seen[name] += 1
            unique_names.append(f"{name}_dup{seen[name]}")
        else:
            seen[name] = 0
            unique_names.append(name)

    # Extract sample barcode from MANIFEST.txt
    manifest_path = download_dir / "MANIFEST.txt"
    file_to_barcode = _build_file_to_barcode_map(manifest_path, tsv_files)

    # Build matrix: one column per sample
    sample_data: dict[str, np.ndarray] = {}
    for tsv_path in tsv_files:
        sample_id = file_to_barcode.get(tsv_path.name)
        if not sample_id:
            # Fall back to parent directory UUID
            sample_id = tsv_path.parent.name

        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        df = df[df["gene_id"].str.startswith("ENSG")]
        if gene_type_filter:
            df = df[df["gene_type"] == gene_type_filter]

        if value_column not in df.columns:
            msg = f"Column {value_column} not found in {tsv_path.name}"
            raise ValueError(msg)

        values = pd.to_numeric(df[value_column], errors="coerce").values
        if len(values) == len(unique_names):
            sample_data[sample_id] = values
        else:
            logger.info(
                "Gene count mismatch, skipping file",
                extra={"file": tsv_path.name, "expected": len(unique_names), "got": len(values)},
            )

    expression = pd.DataFrame(sample_data, index=unique_names).T
    logger.info(
        "Expression matrix assembled",
        extra={"samples": expression.shape[0], "genes": expression.shape[1]},
    )
    return expression


def _build_file_to_barcode_map(
    manifest_path: Path, tsv_files: list[Path]
) -> dict[str, str]:
    """Map TSV filenames to TCGA patient barcodes (12-char) via the GDC API.

    The manifest only contains file UUIDs, not TCGA barcodes.  We query the
    GDC ``/files`` endpoint in batches to resolve each file UUID to its
    case ``submitter_id`` (e.g. ``TCGA-2J-AAB1``).

    A local cache (``_uuid_barcode_cache.json``) next to the manifest avoids
    repeated API calls on subsequent runs.

    Args:
        manifest_path: Path to MANIFEST.txt inside the expression directory.
        tsv_files: List of STAR count TSV paths (used for UUID extraction).

    Returns:
        Dict mapping TSV **filename** (not full path) to 12-char TCGA barcode.
    """
    cache_path = manifest_path.parent / "_uuid_barcode_cache.json"
    if cache_path.exists():
        with open(cache_path) as fh:
            cached: dict[str, str] = json.load(fh)
        logger.info("Loaded UUID→barcode cache", extra={"n_entries": len(cached)})
        return cached

    # Collect file UUIDs — the parent directory name of each TSV
    uuid_to_filename: dict[str, str] = {}
    for tsv in tsv_files:
        file_uuid = tsv.parent.name
        uuid_to_filename[file_uuid] = tsv.name

    uuids = list(uuid_to_filename.keys())
    uuid_to_barcode: dict[str, str] = {}

    # Query GDC API in batches of 100
    batch_size = 100
    for start in range(0, len(uuids), batch_size):
        batch = uuids[start : start + batch_size]
        payload = json.dumps({
            "filters": {
                "op": "in",
                "content": {"field": "file_id", "value": batch},
            },
            "fields": "file_id,cases.submitter_id",
            "size": batch_size,
        }).encode()

        req = Request(
            "https://api.gdc.cancer.gov/files",
            data=payload,
            headers={"Content-Type": "application/json"},
        )
        try:
            with urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode())
        except (URLError, TimeoutError) as exc:
            logger.warning("GDC API request failed", extra={"error": str(exc)})
            return {}

        for hit in data.get("data", {}).get("hits", []):
            file_id = hit.get("file_id", "")
            cases = hit.get("cases", [])
            if cases:
                barcode = cases[0].get("submitter_id", "")
                if barcode:
                    uuid_to_barcode[file_id] = barcode[:12]

    logger.info(
        "GDC UUID→barcode mapping complete",
        extra={"resolved": len(uuid_to_barcode), "total": len(uuids)},
    )

    # Build filename -> barcode map
    result: dict[str, str] = {}
    for uid, fname in uuid_to_filename.items():
        barcode = uuid_to_barcode.get(uid)
        if barcode:
            result[fname] = barcode

    # Cache for next run
    try:
        with open(cache_path, "w") as fh:
            json.dump(result, fh, indent=2)
    except OSError:
        pass

    return result


def load_gdc_clinical(clinical_path: Path) -> pd.DataFrame:
    """Load and clean TCGA-PAAD clinical data from GDC download.

    Handles the GDC aggregated clinical.tsv format with dot-separated
    column names (e.g., 'demographic.vital_status', 'diagnoses.ajcc_pathologic_stage').

    Also handles JSON and simpler CSV formats.

    Args:
        clinical_path: Path to clinical TSV/CSV/JSON file downloaded from GDC.

    Returns:
        DataFrame with columns: patient_id, os_days, os_event, age,
        stage_numeric, ajcc_stage. One row per unique patient.
    """
    suffix = clinical_path.suffix.lower()
    if suffix == ".json":
        df = _load_gdc_json_clinical(clinical_path)
    elif suffix in {".tsv", ".csv"}:
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(clinical_path, sep=sep, low_memory=False)
    else:
        msg = f"Unsupported clinical file format: {suffix}"
        raise ValueError(msg)

    logger.info("Raw clinical loaded", extra={"rows": len(df), "cols": len(df.columns)})

    # Detect GDC dot-separated format
    if any("." in c for c in df.columns[:5]):
        return _standardize_gdc_dotted_clinical(df)
    return _standardize_clinical(df)


def _load_gdc_json_clinical(path: Path) -> pd.DataFrame:
    """Parse GDC JSON clinical export into a flat DataFrame."""
    import json

    with path.open() as f:
        records = json.load(f)

    rows = []
    for rec in records:
        row: dict = {"submitter_id": rec.get("submitter_id")}
        demo = rec.get("demographic", {})
        row["vital_status"] = demo.get("vital_status")
        row["days_to_death"] = demo.get("days_to_death")
        row["days_to_last_follow_up"] = demo.get("days_to_last_follow_up")
        row["age_at_diagnosis"] = demo.get("age_at_diagnosis")

        diagnoses = rec.get("diagnoses", [{}])
        if diagnoses:
            row["ajcc_pathologic_stage"] = diagnoses[0].get("ajcc_pathologic_stage")
            row["days_to_last_follow_up"] = (
                row["days_to_last_follow_up"]
                or diagnoses[0].get("days_to_last_follow_up")
            )
        rows.append(row)

    return pd.DataFrame(rows)


def _standardize_gdc_dotted_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Handle GDC aggregated clinical.tsv with dot-separated column names.

    The GDC clinical TSV has multiple rows per patient (one per treatment/diagnosis).
    We deduplicate to one row per patient, preferring the primary diagnosis row.

    Key columns:
      cases.submitter_id, demographic.vital_status, demographic.days_to_death,
      diagnoses.days_to_last_follow_up, diagnoses.age_at_diagnosis,
      diagnoses.ajcc_pathologic_stage
    """
    def _clean_val(val: object) -> object:
        """Replace GDC missing sentinels with NaN."""
        if isinstance(val, str) and val.strip() in GDC_MISSING:
            return np.nan
        return val

    df = df.map(_clean_val)

    # Map to standard column names
    col_map = {
        "cases.submitter_id": "patient_id",
        "demographic.vital_status": "vital_status",
        "demographic.days_to_death": "days_to_death",
        "diagnoses.days_to_last_follow_up": "days_to_last_follow_up",
        "diagnoses.age_at_diagnosis": "age_at_diagnosis",
        "diagnoses.ajcc_pathologic_stage": "ajcc_pathologic_stage",
    }

    # Keep only columns we need
    available = {k: v for k, v in col_map.items() if k in df.columns}
    subset = df[list(available.keys())].rename(columns=available)

    # Convert numeric columns
    for col in ["days_to_death", "days_to_last_follow_up", "age_at_diagnosis"]:
        if col in subset.columns:
            subset[col] = pd.to_numeric(subset[col], errors="coerce")

    # Deduplicate: keep one row per patient with the most complete primary diagnosis data
    # Prefer rows where ajcc_pathologic_stage is not NaN
    subset["_has_stage"] = subset["ajcc_pathologic_stage"].notna().astype(int)
    subset = subset.sort_values(["patient_id", "_has_stage"], ascending=[True, False])
    subset = subset.drop_duplicates(subset="patient_id", keep="first")
    subset = subset.drop(columns=["_has_stage"])

    # Build output
    out = pd.DataFrame()
    out["patient_id"] = subset["patient_id"].astype(str)

    # OS days: days_to_death for dead patients, days_to_last_follow_up for alive
    out["os_days"] = subset["days_to_death"].fillna(subset["days_to_last_follow_up"])

    # OS event: 1 = dead, 0 = alive/censored
    out["os_event"] = (
        subset["vital_status"].astype(str).str.lower().str.strip()
        .map({"dead": 1.0, "alive": 0.0})
    )

    # Age in years (GDC stores age_at_diagnosis in days)
    age_raw = subset["age_at_diagnosis"]
    if age_raw.median(skipna=True) > 365:
        out["age"] = (age_raw / 365.25).round(1)
    else:
        out["age"] = age_raw

    # AJCC stage
    out["ajcc_stage"] = subset["ajcc_pathologic_stage"].astype(str)
    out["stage_numeric"] = out["ajcc_stage"].map(_stage_to_numeric)

    # Drop rows missing survival endpoint
    pre = len(out)
    out = out.dropna(subset=["os_days", "os_event"])
    out = out[out["os_days"] > 0]

    logger.info(
        "GDC clinical standardized",
        extra={"kept": len(out), "dropped": pre - len(out), "unique_patients": out["patient_id"].nunique()},
    )
    return out.reset_index(drop=True)


def _standardize_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Map various non-GDC column naming conventions to our standard schema."""
    col_map = _detect_column_mapping(df)
    out = pd.DataFrame()

    out["patient_id"] = df[col_map["patient_id"]].astype(str)
    out["os_days"] = _compute_os_days(df, col_map)
    out["os_event"] = _compute_os_event(df, col_map)

    age_raw = pd.to_numeric(df[col_map["age"]], errors="coerce")
    if age_raw.median() > 365:
        out["age"] = (age_raw / 365.25).round(1)
    else:
        out["age"] = age_raw

    out["ajcc_stage"] = df[col_map["stage"]].astype(str)
    out["stage_numeric"] = out["ajcc_stage"].map(_stage_to_numeric)

    pre = len(out)
    out = out.dropna(subset=["os_days", "os_event"])
    out = out[out["os_days"] > 0]
    logger.info(
        "Clinical standardized",
        extra={"kept": len(out), "dropped": pre - len(out)},
    )
    return out.reset_index(drop=True)


def _detect_column_mapping(df: pd.DataFrame) -> dict[str, str]:
    """Auto-detect which columns map to our schema."""
    cols = set(df.columns.str.lower())
    mapping: dict[str, str] = {}

    for candidate in ["submitter_id", "bcr_patient_barcode", "patient_id", "case_id"]:
        if candidate in df.columns:
            mapping["patient_id"] = candidate
            break
        if candidate in cols:
            match = _find_col_case_insensitive(df, candidate)
            if match:
                mapping["patient_id"] = match
                break
    if "patient_id" not in mapping:
        msg = f"Cannot find patient ID column in: {list(df.columns)}"
        raise ValueError(msg)

    for candidate in ["vital_status", "vitalstatus"]:
        match = _find_col_case_insensitive(df, candidate)
        if match:
            mapping["vital_status"] = match
            break

    for candidate in ["days_to_death", "daystodeath"]:
        match = _find_col_case_insensitive(df, candidate)
        if match:
            mapping["days_to_death"] = match
            break

    for candidate in ["days_to_last_follow_up", "days_to_last_followup", "daystolastfollowup"]:
        match = _find_col_case_insensitive(df, candidate)
        if match:
            mapping["days_to_last_follow_up"] = match
            break

    for candidate in ["age_at_diagnosis", "age_at_initial_pathologic_diagnosis", "age"]:
        match = _find_col_case_insensitive(df, candidate)
        if match:
            mapping["age"] = match
            break
    if "age" not in mapping:
        msg = f"Cannot find age column in: {list(df.columns)}"
        raise ValueError(msg)

    for candidate in [
        "ajcc_pathologic_stage", "ajcc_staging_edition",
        "pathologic_stage", "clinical_stage", "tumor_stage",
    ]:
        match = _find_col_case_insensitive(df, candidate)
        if match:
            mapping["stage"] = match
            break
    if "stage" not in mapping:
        msg = f"Cannot find stage column in: {list(df.columns)}"
        raise ValueError(msg)

    logger.info("Column mapping detected", extra={"mapping": mapping})
    return mapping


def _find_col_case_insensitive(df: pd.DataFrame, target: str) -> str | None:
    """Find a column name case-insensitively."""
    lower_map = {c.lower(): c for c in df.columns}
    return lower_map.get(target.lower())


def _compute_os_days(df: pd.DataFrame, col_map: dict[str, str]) -> pd.Series:
    """Compute overall survival days from death/follow-up columns."""
    dtd = pd.to_numeric(
        df.get(col_map.get("days_to_death", ""), pd.Series(dtype=float)),
        errors="coerce",
    )
    dtf = pd.to_numeric(
        df.get(col_map.get("days_to_last_follow_up", ""), pd.Series(dtype=float)),
        errors="coerce",
    )
    return dtd.fillna(dtf)


def _compute_os_event(df: pd.DataFrame, col_map: dict[str, str]) -> pd.Series:
    """Compute event indicator: 1 = dead, 0 = alive/censored."""
    if "vital_status" not in col_map:
        return pd.Series(np.nan, index=df.index)

    vs = df[col_map["vital_status"]].astype(str).str.lower().str.strip()
    return vs.map({"dead": 1, "alive": 0, "1": 1, "0": 0}).astype(float)


def _stage_to_numeric(stage: str) -> float | None:
    """Map AJCC stage string to numeric ordinal.

    Args:
        stage: AJCC stage string (e.g., 'Stage IIA', 'Stage IV').

    Returns:
        Numeric ordinal 1-4, or None if unparseable.
    """
    if not isinstance(stage, str):
        return None
    s = stage.upper().strip().replace("STAGE ", "")
    stage_map = {
        "I": 1, "IA": 1, "IB": 1,
        "II": 2, "IIA": 2, "IIB": 2,
        "III": 3, "IIIA": 3, "IIIB": 3,
        "IV": 4, "IVA": 4, "IVB": 4,
    }
    return stage_map.get(s)


def load_expression_matrix(expression_path: Path) -> pd.DataFrame:
    """Load a gene expression matrix (genes x samples or samples x genes).

    Handles TSV/CSV files. Detects orientation and transposes if needed
    so that rows = samples, columns = genes.

    Args:
        expression_path: Path to expression TSV/CSV file.

    Returns:
        DataFrame with rows = samples, columns = gene symbols.
    """
    suffix = expression_path.suffix.lower()
    sep = "\t" if suffix == ".tsv" else ","
    df = pd.read_csv(expression_path, sep=sep, index_col=0)

    # Heuristic: if more columns than rows, it's genes-as-rows (standard GDC format)
    if df.shape[1] > df.shape[0]:
        logger.info(
            "Expression appears genes-as-rows; transposing",
            extra={"genes": df.shape[0], "samples": df.shape[1]},
        )
        df = df.T

    logger.info(
        "Expression matrix loaded",
        extra={"samples": df.shape[0], "genes": df.shape[1]},
    )
    return df


def build_mutation_matrix_from_maf_dir(
    maf_dir: Path,
    min_frequency: float = 0.02,
    exclude_silent: bool = True,
) -> pd.DataFrame:
    """Build a binary gene-level mutation matrix from per-sample MAF.gz files.

    Args:
        maf_dir: Directory containing UUID subdirectories, each with a .maf.gz file.
        min_frequency: Minimum mutation frequency across patients to include a gene.
            Genes mutated in fewer than this fraction of patients are dropped.
        exclude_silent: If True, exclude Silent/Intron/IGR/3'UTR/5'UTR variants.

    Returns:
        Binary DataFrame of shape (n_patients, n_genes). 1 = gene mutated in patient.
        Index = patient barcodes (12-char TCGA IDs). Columns = gene symbols.
    """
    maf_files = sorted(maf_dir.rglob("*.maf.gz"))
    if not maf_files:
        msg = f"No MAF.gz files found in {maf_dir}"
        raise FileNotFoundError(msg)

    logger.info("Building mutation matrix from MAF files", extra={"n_files": len(maf_files)})

    silent_classes = {
        "Silent", "Intron", "IGR", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
        "RNA", "lincRNA",
    }

    records: list[tuple[str, str]] = []
    for maf_path in maf_files:
        with gzip.open(maf_path, "rt") as f:
            df = pd.read_csv(
                f, sep="\t", comment="#", low_memory=False,
                usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"],
            )

        if exclude_silent:
            df = df[~df["Variant_Classification"].isin(silent_classes)]

        df["patient_id"] = df["Tumor_Sample_Barcode"].str[:12]
        for _, row in df[["patient_id", "Hugo_Symbol"]].drop_duplicates().iterrows():
            records.append((row["patient_id"], row["Hugo_Symbol"]))

    if not records:
        msg = "No mutations found after filtering"
        raise ValueError(msg)

    mutations_df = pd.DataFrame(records, columns=["patient_id", "gene"])
    mutations_df = mutations_df.drop_duplicates()

    # Pivot to binary matrix
    mutations_df["mutated"] = 1
    matrix = mutations_df.pivot_table(
        index="patient_id", columns="gene", values="mutated", fill_value=0,
    ).astype(int)

    # Filter by minimum frequency
    n_patients = len(matrix)
    gene_freq = matrix.sum(axis=0) / n_patients
    kept_genes = gene_freq[gene_freq >= min_frequency].index
    matrix = matrix[kept_genes]

    logger.info(
        "Mutation matrix built",
        extra={
            "n_patients": n_patients,
            "n_genes_total": len(gene_freq),
            "n_genes_kept": len(kept_genes),
            "min_frequency": min_frequency,
        },
    )
    return matrix


def compute_tmb_from_maf_dir(maf_dir: Path) -> pd.Series:
    """Compute tumor mutation burden from a directory of per-sample MAF.gz files.

    Args:
        maf_dir: Directory containing UUID subdirectories, each with a .maf.gz file.

    Returns:
        Series indexed by patient barcode (first 12 chars of Tumor_Sample_Barcode)
        with log2(mutation_count + 1) values.
    """
    maf_files = sorted(maf_dir.rglob("*.maf.gz"))
    if not maf_files:
        msg = f"No MAF.gz files found in {maf_dir}"
        raise FileNotFoundError(msg)

    logger.info("Computing TMB from MAF files", extra={"n_files": len(maf_files)})

    all_counts: dict[str, int] = {}
    for maf_path in maf_files:
        with gzip.open(maf_path, "rt") as f:
            df = pd.read_csv(f, sep="\t", comment="#", low_memory=False, usecols=["Tumor_Sample_Barcode"])

        if "Tumor_Sample_Barcode" not in df.columns:
            continue

        # Extract patient-level barcode (first 12 chars: TCGA-XX-XXXX)
        barcodes = df["Tumor_Sample_Barcode"].str[:12]
        for patient_id, count in barcodes.value_counts().items():
            all_counts[patient_id] = all_counts.get(patient_id, 0) + count

    counts = pd.Series(all_counts)
    tmb_log = np.log2(counts + 1)
    tmb_log.name = "tmb_log"
    tmb_log.index.name = "patient_id"

    logger.info(
        "TMB computed",
        extra={"n_patients": len(tmb_log), "median_mutations": int(counts.median())},
    )
    return tmb_log


def compute_tmb(maf_path: Path) -> pd.Series:
    """Compute tumor mutation burden from a single concatenated MAF file.

    Args:
        maf_path: Path to MAF file (tab-separated, standard GDC format).

    Returns:
        Series indexed by patient barcode with log2(mutation_count + 1) values.
    """
    maf = pd.read_csv(maf_path, sep="\t", comment="#", low_memory=False)
    barcode_col = "Tumor_Sample_Barcode"
    if barcode_col not in maf.columns:
        msg = f"MAF missing {barcode_col} column"
        raise ValueError(msg)

    maf["patient_id"] = maf[barcode_col].str[:12]
    counts = maf.groupby("patient_id").size()
    tmb_log = np.log2(counts + 1)
    tmb_log.name = "tmb_log"
    logger.info(
        "TMB computed",
        extra={"n_patients": len(tmb_log), "median_mutations": int(counts.median())},
    )
    return tmb_log
