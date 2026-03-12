"""H13h: Supplementary analyses addressing reviewer concerns.

Three analyses:
  1. RAS mutation interaction model (drug_sens ~ SLC2A8 + KRAS + SLC2A8:KRAS)
  2. Classification performance (ROC AUC, top-quartile sensitivity difference)
  3. TCGA-PAAD descriptive (SLC2A8 distribution, KRAS co-occurrence in tumors)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

from pdac.h13_trehalose_vulnerability.data import (
    get_ras_mutant_model_ids,
    load_crispr_dependencies,
    load_expression,
    load_model_metadata,
    load_mutations,
)
from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection
from pdac.h13_trehalose_vulnerability.h13b_prism import (
    _find_compound_col,
    load_prism_compound_info,
    load_prism_sensitivity,
)
from pdac.h13_trehalose_vulnerability.scoring import (
    compute_vulnerability_score,
    stratify_by_target_and_ras,
)
from pdac.shared.logging import get_logger

logger = get_logger(__name__)

DEPMAP_DIR = Path("data/raw/depmap")
TCGA_DIR = Path("data/processed")
RESULTS_DIR = Path("results/h13h")


def _load_depmap() -> dict:
    """Load DepMap data."""
    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    crispr = load_crispr_dependencies(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)

    gene_sets = GeneSetCollection()
    ras_mutant_ids = get_ras_mutant_model_ids(mutations, gene_sets.ras_mutations)

    common_ids = sorted(
        set(expression.index) & set(crispr.index) & set(metadata.index)
    )
    expression = expression.loc[common_ids]
    crispr = crispr.loc[common_ids]

    scores_df = compute_vulnerability_score(
        expression, gene_sets.target, gene_sets.scoring_sets(),
    )
    quadrant_labels = stratify_by_target_and_ras(scores_df, ras_mutant_ids)
    target_ids = quadrant_labels[
        quadrant_labels == "SLC2A8-high / RAS-mut"
    ].index.tolist()

    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)
    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])

    return {
        "expression": expression,
        "crispr": crispr,
        "metadata": metadata,
        "mutations": mutations,
        "sensitivity": sensitivity,
        "compound_info": compound_info,
        "name_col": name_col,
        "id_col": id_col,
        "target_ids": target_ids,
        "ras_mutant_ids": ras_mutant_ids,
        "common_ids": common_ids,
    }


def _resolve_drug(drug_name: str, data: dict) -> pd.Series | None:
    """Look up drug sensitivity vector by name."""
    ci = data["compound_info"]
    mask = ci[data["name_col"]].str.upper() == drug_name.upper()
    if mask.sum() == 0:
        return None
    for brd_id in ci.loc[mask, data["id_col"]].tolist():
        if brd_id in data["sensitivity"].columns:
            return data["sensitivity"][brd_id].dropna()
    return None


# =====================================================================
# Analysis 1: RAS mutation interaction model
# =====================================================================

def analysis_1_ras_interaction(data: dict) -> None:
    """Test drug_sensitivity ~ SLC2A8 + KRAS_mut + SLC2A8:KRAS_mut."""
    logger.info("=" * 70)
    logger.info("ANALYSIS 1: RAS mutation interaction model")
    logger.info("=" * 70)

    expr = data["expression"]
    ras_ids = set(data["ras_mutant_ids"])

    drugs = {
        "Afuresertib": "GSK2110183",
        "Ipatasertib": "GDC-0068",
        "Capivasertib": "AZD5363",
    }

    for label, name in drugs.items():
        sens = _resolve_drug(name, data)
        if sens is None:
            logger.info("  %s: not found", label)
            continue

        shared = sorted(set(expr.index) & set(sens.index))
        slc2a8 = expr.loc[shared, "SLC2A8"].values
        drug_val = sens[shared].values
        kras_mut = np.array([1.0 if s in ras_ids else 0.0 for s in shared])

        # Standardize continuous predictors
        slc2a8_z = (slc2a8 - slc2a8.mean()) / slc2a8.std()

        # Build design matrix: intercept, SLC2A8, KRAS, SLC2A8*KRAS
        n = len(shared)
        x_mat = np.column_stack([
            np.ones(n),
            slc2a8_z,
            kras_mut,
            slc2a8_z * kras_mut,
        ])

        # OLS
        beta, residuals, rank, sv = np.linalg.lstsq(x_mat, drug_val, rcond=None)
        y_hat = x_mat @ beta
        ss_res = np.sum((drug_val - y_hat) ** 2)
        df_res = n - x_mat.shape[1]
        mse = ss_res / df_res

        # Standard errors and p-values
        cov = mse * np.linalg.inv(x_mat.T @ x_mat)
        se = np.sqrt(np.diag(cov))
        t_stat = beta / se
        p_vals = 2 * stats.t.sf(np.abs(t_stat), df_res)

        names = ["Intercept", "SLC2A8", "KRAS_mut", "SLC2A8:KRAS"]
        logger.info("  %s (n=%d):", label, n)
        for nm, b, s, t, p in zip(names, beta, se, t_stat, p_vals):
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            logger.info("    %-20s  beta=%+.4f  SE=%.4f  t=%.3f  p=%.4f %s",
                        nm, b, s, t, p, sig)

        # Also test SLC2A8 alone vs SLC2A8+KRAS (F-test for KRAS contribution)
        x_reduced = np.column_stack([np.ones(n), slc2a8_z])
        beta_r, _, _, _ = np.linalg.lstsq(x_reduced, drug_val, rcond=None)
        ss_reduced = np.sum((drug_val - x_reduced @ beta_r) ** 2)
        f_stat = ((ss_reduced - ss_res) / 2) / mse
        p_f = 1 - stats.f.cdf(f_stat, 2, df_res)
        logger.info("    F-test (KRAS terms): F=%.3f, p=%.4f", f_stat, p_f)
        logger.info("")


# =====================================================================
# Analysis 2: Classification performance (ROC AUC)
# =====================================================================

def analysis_2_classification(data: dict) -> None:
    """ROC AUC for SLC2A8 predicting AKT drug sensitivity."""
    logger.info("=" * 70)
    logger.info("ANALYSIS 2: Classification performance (ROC AUC)")
    logger.info("=" * 70)

    expr = data["expression"]

    drugs = {
        "Afuresertib": "GSK2110183",
        "Ipatasertib": "GDC-0068",
        "Capivasertib": "AZD5363",
    }

    for label, name in drugs.items():
        sens = _resolve_drug(name, data)
        if sens is None:
            logger.info("  %s: not found", label)
            continue

        shared = sorted(set(expr.index) & set(sens.index))
        slc2a8 = expr.loc[shared, "SLC2A8"].values
        drug_val = sens[shared].values

        # Binary: sensitive = below median (more negative = more killing)
        is_sensitive = (drug_val < np.median(drug_val)).astype(int)

        # ROC AUC (higher SLC2A8 -> more sensitive -> more negative drug value)
        # SLC2A8 predicts sensitivity (negative correlation), so use SLC2A8 directly
        auc = roc_auc_score(is_sensitive, slc2a8)
        fpr, tpr, _ = roc_curve(is_sensitive, slc2a8)

        logger.info("  %s:", label)
        logger.info("    ROC AUC (median split): %.3f", auc)

        # Top-quartile vs bottom-quartile sensitivity difference
        q75 = np.percentile(slc2a8, 75)
        q25 = np.percentile(slc2a8, 25)
        top_q = drug_val[slc2a8 >= q75]
        bot_q = drug_val[slc2a8 <= q25]
        d_q = (top_q.mean() - bot_q.mean()) / np.sqrt(
            (top_q.std() ** 2 + bot_q.std() ** 2) / 2
        )
        _, p_q = stats.ttest_ind(top_q, bot_q)
        logger.info("    Top-quartile mean sensitivity: %.4f (n=%d)", top_q.mean(), len(top_q))
        logger.info("    Bottom-quartile mean sensitivity: %.4f (n=%d)", bot_q.mean(), len(bot_q))
        logger.info("    Quartile difference: Cohen's d=%.3f, p=%.4f", d_q, p_q)

        # Also: 25th percentile threshold (more stringent)
        is_sensitive_25 = (drug_val < np.percentile(drug_val, 25)).astype(int)
        auc_25 = roc_auc_score(is_sensitive_25, slc2a8)
        logger.info("    ROC AUC (25th pct split): %.3f", auc_25)
        logger.info("")

    # Save ROC data for figure generation
    roc_data = {}
    for label, name in drugs.items():
        sens = _resolve_drug(name, data)
        if sens is None:
            continue
        shared = sorted(set(expr.index) & set(sens.index))
        slc2a8 = expr.loc[shared, "SLC2A8"].values
        drug_val = sens[shared].values
        is_sensitive = (drug_val < np.median(drug_val)).astype(int)
        fpr, tpr, _ = roc_curve(is_sensitive, slc2a8)
        auc = roc_auc_score(is_sensitive, slc2a8)
        roc_data[label] = {"fpr": fpr.tolist(), "tpr": tpr.tolist(), "auc": float(auc)}

    import json
    roc_path = RESULTS_DIR / "roc_data.json"
    with open(roc_path, "w") as f:
        json.dump(roc_data, f, indent=2)
    logger.info("  ROC data saved to %s", roc_path)


# =====================================================================
# Analysis 3: TCGA-PAAD descriptive
# =====================================================================

def analysis_3_tcga(data: dict) -> None:
    """TCGA-PAAD SLC2A8 expression distribution and KRAS co-occurrence."""
    logger.info("=" * 70)
    logger.info("ANALYSIS 3: TCGA-PAAD descriptive analysis")
    logger.info("=" * 70)

    tcga_expr_path = TCGA_DIR / "tcga_paad_expression_tpm.parquet"
    if not tcga_expr_path.exists():
        logger.info("  TCGA expression parquet not found at %s", tcga_expr_path)
        return

    tcga_expr = pd.read_parquet(tcga_expr_path)
    logger.info("  TCGA-PAAD expression matrix: %d patients x %d genes",
                *tcga_expr.shape)

    # Check if log-transformed
    if tcga_expr.max().max() > 100:
        logger.info("  Applying log2(TPM+1) transformation")
        tcga_expr = np.log2(tcga_expr + 1)

    if "SLC2A8" not in tcga_expr.columns:
        logger.info("  SLC2A8 not found in TCGA expression columns")
        return

    slc2a8_tcga = tcga_expr["SLC2A8"].dropna()
    logger.info("  SLC2A8 expression: n=%d, mean=%.3f, median=%.3f, std=%.3f",
                len(slc2a8_tcga), slc2a8_tcga.mean(), slc2a8_tcga.median(),
                slc2a8_tcga.std())

    # Distribution stats
    q25 = slc2a8_tcga.quantile(0.25)
    q75 = slc2a8_tcga.quantile(0.75)
    logger.info("  SLC2A8 IQR: [%.3f, %.3f]", q25, q75)

    # Compare to DepMap PDAC cell lines
    depmap_expr = data["expression"]
    meta = data["metadata"]
    lineage_col = "OncotreeLineage"
    if lineage_col in meta.columns:
        pdac_ids = meta.index[meta[lineage_col].str.contains("Pancrea", na=False)]
        pdac_expr = depmap_expr.loc[depmap_expr.index.isin(pdac_ids), "SLC2A8"].dropna()
        logger.info("  DepMap PDAC SLC2A8: n=%d, mean=%.3f, median=%.3f",
                    len(pdac_expr), pdac_expr.mean(), pdac_expr.median())

        # Test distribution similarity
        stat_ks, p_ks = stats.ks_2samp(slc2a8_tcga.values, pdac_expr.values)
        logger.info("  KS test (TCGA vs DepMap PDAC): stat=%.3f, p=%.4f", stat_ks, p_ks)

    # KRAS mutation co-occurrence (if available)
    if "KRAS" in tcga_expr.columns:
        kras_expr = tcga_expr["KRAS"].dropna()
        shared = sorted(set(slc2a8_tcga.index) & set(kras_expr.index))
        rho, p = stats.spearmanr(slc2a8_tcga[shared], kras_expr[shared])
        logger.info("  SLC2A8 vs KRAS expression: rho=%.3f, p=%.4f (n=%d)",
                    rho, p, len(shared))

    # AKT1/MTOR co-expression in TCGA tumors
    for gene in ["AKT1", "MTOR", "FASN", "VPS35", "SCD"]:
        if gene in tcga_expr.columns:
            gene_expr = tcga_expr[gene].dropna()
            shared = sorted(set(slc2a8_tcga.index) & set(gene_expr.index))
            if len(shared) > 20:
                rho, p = stats.spearmanr(slc2a8_tcga[shared], gene_expr[shared])
                logger.info("  SLC2A8 vs %s (TCGA): rho=%.3f, p=%.4f (n=%d)",
                            gene, rho, p, len(shared))

    # Save TCGA distribution for figure
    import json
    tcga_data = {
        "slc2a8_values": slc2a8_tcga.values.tolist(),
        "n_patients": len(slc2a8_tcga),
        "mean": float(slc2a8_tcga.mean()),
        "median": float(slc2a8_tcga.median()),
    }
    tcga_path = RESULTS_DIR / "tcga_descriptive.json"
    with open(tcga_path, "w") as f:
        json.dump(tcga_data, f, indent=2)
    logger.info("  TCGA data saved to %s", tcga_path)


# =====================================================================
# Main
# =====================================================================

def main() -> None:
    """Run all reviewer analyses."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("H13h: Supplementary reviewer analyses")
    logger.info("=" * 70)

    data = _load_depmap()

    analysis_1_ras_interaction(data)
    analysis_2_classification(data)
    analysis_3_tcga(data)

    logger.info("=" * 70)
    logger.info("All reviewer analyses complete. Results in %s", RESULTS_DIR)


if __name__ == "__main__":
    main()
