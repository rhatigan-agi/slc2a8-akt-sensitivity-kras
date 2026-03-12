"""H13e Sanity Checks: PRISM orientation, within-group correlation, lineage control.

Three checks before trusting the H13d Analysis 4 paradox:
  1. Confirm PRISM sign convention (lower = more killing)
  2. Correlation within target group only (does paradox disappear?)
  3. Lineage-controlled regression (drug_response ~ SLC2A8 + lineage)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

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

# PI3K/AKT core genes (excluding PTEN)
PI3K_CORE_GENES: list[str] = [
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
    "AKT1", "AKT2", "AKT3", "MTOR", "RPTOR", "RICTOR",
    "RPS6KB1", "RPS6KB2", "EIF4EBP1",
]

KEY_DRUGS: list[str] = [
    "ALPELISIB", "GSK2110183", "MK-2206", "GDC-0077",
    "AZD8835", "TASELISIB", "GDC-0068", "AZD5363",
]


def _load_all() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame,
    list[str], list[str], pd.Series, pd.Series,
]:
    """Load data and return expression, sensitivity, compound_info,
    target_ids, background_ids, pathway_score, quadrant_labels."""
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

    scores_df = compute_vulnerability_score(
        expression, gene_sets.target, gene_sets.scoring_sets(),
    )
    quadrant_labels = stratify_by_target_and_ras(scores_df, ras_mutant_ids)

    target_ids = quadrant_labels[
        quadrant_labels == "SLC2A8-high / RAS-mut"
    ].index.tolist()
    background_ids = [i for i in common_ids if i not in target_ids]

    # PI3K pathway z-score
    pi3k_genes = [g for g in PI3K_CORE_GENES if g in expression.columns]
    subset = expression[pi3k_genes]
    zscored = subset.apply(stats.zscore, nan_policy="omit")
    pathway_score = zscored.mean(axis=1)
    pathway_score.name = "pi3k_zscore"

    return (
        expression, metadata, scores_df,
        target_ids, background_ids, pathway_score, quadrant_labels,
    )


def _resolve_drug_sensitivity(
    drug_name: str,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> pd.Series | None:
    """Look up a drug's sensitivity vector by name."""
    mask = compound_info[name_col].str.upper() == drug_name.upper()
    if mask.sum() == 0:
        return None
    brd_ids = compound_info.loc[mask, id_col].tolist()
    for brd_id in brd_ids:
        if brd_id in sensitivity.columns:
            return sensitivity[brd_id].dropna()
    return None


def _sig(p: float) -> str:
    """Significance stars."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.1:
        return "."
    return ""


def check_1_prism_orientation(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> None:
    """Verify PRISM sign convention: lower values should = more killing."""
    border = "=" * 90
    logger.info(border)
    logger.info("CHECK 1: PRISM SIGN CONVENTION")
    logger.info(border)

    # Use a well-known cytotoxic drug as positive control
    # Bortezomib (proteasome inhibitor) should kill most cells
    control_drugs = ["BORTEZOMIB", "STAUROSPORINE", "DOXORUBICIN", "PACLITAXEL"]

    for drug in control_drugs:
        sens = _resolve_drug_sensitivity(
            drug, sensitivity, compound_info, name_col, id_col,
        )
        if sens is not None:
            logger.info(
                f"  {drug}: mean={sens.mean():.4f}, median={sens.median():.4f}, "
                f"min={sens.min():.4f}, max={sens.max():.4f}, n={len(sens)}"
            )
            pct_negative = (sens < 0).mean() * 100
            logger.info(
                f"    {pct_negative:.1f}% of values are negative (expected: most, if lower=more killing)"
            )
        else:
            logger.info(f"  {drug}: NOT FOUND")

    # Also check: overall distribution of all PRISM values
    all_vals = sensitivity.values.flatten()
    all_vals = all_vals[~np.isnan(all_vals)]
    logger.info(f"  Overall PRISM distribution:")
    logger.info(f"    mean={np.mean(all_vals):.4f}, median={np.median(all_vals):.4f}")
    logger.info(f"    10th pctl={np.percentile(all_vals, 10):.4f}, 90th={np.percentile(all_vals, 90):.4f}")
    logger.info(f"    % negative={100 * (all_vals < 0).mean():.1f}%")
    logger.info(f"    Interpretation: {'CONFIRMED lower=more killing' if np.mean(all_vals) < 0 else 'CHECK — mean is positive, orientation may differ'}")


def check_2_within_group_correlation(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    expression: pd.DataFrame,
    pathway_score: pd.Series,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Test correlation WITHIN target group vs WITHIN background."""
    border = "=" * 90
    logger.info(border)
    logger.info("CHECK 2: WITHIN-GROUP CORRELATION (does paradox disappear?)")
    logger.info(border)

    logger.info(f"  {'Drug':<20} {'All (rho)':>10} {'Target (rho)':>14} {'BG (rho)':>10} {'Target p':>12} {'Target n':>10}")
    logger.info("  " + "-" * 90)

    for drug_name in KEY_DRUGS:
        sens = _resolve_drug_sensitivity(
            drug_name, sensitivity, compound_info, name_col, id_col,
        )
        if sens is None:
            logger.info(f"  {drug_name:<20} NOT FOUND")
            continue

        # All cells
        shared_all = sorted(set(pathway_score.dropna().index) & set(sens.index))
        if len(shared_all) < 30:
            continue
        rho_all, _ = stats.spearmanr(pathway_score[shared_all], sens[shared_all])

        # Within target only
        shared_target = sorted(
            set(target_ids) & set(pathway_score.dropna().index) & set(sens.index)
        )
        rho_tgt, p_tgt = (np.nan, np.nan)
        if len(shared_target) >= 10:
            rho_tgt, p_tgt = stats.spearmanr(
                pathway_score[shared_target], sens[shared_target],
            )

        # Within background only
        shared_bg = sorted(
            set(background_ids) & set(pathway_score.dropna().index) & set(sens.index)
        )
        rho_bg, _ = (np.nan, np.nan)
        if len(shared_bg) >= 10:
            rho_bg, _ = stats.spearmanr(
                pathway_score[shared_bg], sens[shared_bg],
            )

        logger.info(
            f"  {drug_name:<20} {rho_all:>+10.4f} {rho_tgt:>+14.4f} {rho_bg:>+10.4f} "
            f"{p_tgt:>12.4e} {len(shared_target):>10} {_sig(p_tgt)}"
        )

    # Also test with SLC2A8 expression directly instead of PI3K pathway score
    logger.info("")
    logger.info("  SLC2A8 EXPRESSION vs DRUG SENSITIVITY (continuous):")
    if "SLC2A8" not in expression.columns:
        logger.info("  SLC2A8 not in expression matrix")
        return

    slc2a8 = expression["SLC2A8"]

    logger.info(f"  {'Drug':<20} {'All (rho)':>10} {'Target (rho)':>14} {'BG (rho)':>10} {'All p':>12}")
    logger.info("  " + "-" * 80)

    for drug_name in KEY_DRUGS:
        sens = _resolve_drug_sensitivity(
            drug_name, sensitivity, compound_info, name_col, id_col,
        )
        if sens is None:
            continue

        shared = sorted(set(slc2a8.dropna().index) & set(sens.index))
        if len(shared) < 30:
            continue
        rho_all, p_all = stats.spearmanr(slc2a8[shared], sens[shared])

        shared_tgt = sorted(set(target_ids) & set(slc2a8.dropna().index) & set(sens.index))
        rho_tgt = np.nan
        if len(shared_tgt) >= 10:
            rho_tgt, _ = stats.spearmanr(slc2a8[shared_tgt], sens[shared_tgt])

        shared_bg = sorted(set(background_ids) & set(slc2a8.dropna().index) & set(sens.index))
        rho_bg = np.nan
        if len(shared_bg) >= 10:
            rho_bg, _ = stats.spearmanr(slc2a8[shared_bg], sens[shared_bg])

        logger.info(
            f"  {drug_name:<20} {rho_all:>+10.4f} {rho_tgt:>+14.4f} {rho_bg:>+10.4f} "
            f"{p_all:>12.4e} {_sig(p_all)}"
        )


def check_3_lineage_controlled(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> None:
    """Regression: drug_response ~ SLC2A8 + lineage to control for lineage effects."""
    border = "=" * 90
    logger.info(border)
    logger.info("CHECK 3: LINEAGE-CONTROLLED REGRESSION (drug ~ SLC2A8 + lineage)")
    logger.info(border)

    if "SLC2A8" not in expression.columns:
        logger.info("  SLC2A8 not in expression matrix")
        return

    # Get top lineages for dummy encoding
    lineage_col = "OncotreeLineage"
    lineage_map = metadata[lineage_col].reindex(expression.index).dropna()
    top_lineages = lineage_map.value_counts().head(15).index.tolist()

    logger.info(f"  Top lineages controlled for: {top_lineages[:8]}...")
    logger.info("")

    # For each drug: partial correlation of SLC2A8 controlling for lineage
    logger.info(f"  {'Drug':<20} {'Raw rho':>10} {'Partial r':>12} {'Partial p':>12} {'Interpretation'}")
    logger.info("  " + "-" * 80)

    slc2a8 = expression["SLC2A8"]

    for drug_name in KEY_DRUGS:
        sens = _resolve_drug_sensitivity(
            drug_name, sensitivity, compound_info, name_col, id_col,
        )
        if sens is None:
            logger.info(f"  {drug_name:<20} NOT FOUND")
            continue

        # Build regression DataFrame
        shared = sorted(
            set(slc2a8.dropna().index) & set(sens.index) & set(lineage_map.index)
        )
        if len(shared) < 50:
            continue

        df = pd.DataFrame({
            "drug": sens[shared].values,
            "slc2a8": slc2a8[shared].values,
        }, index=shared)

        # Add lineage dummies
        lin = lineage_map[shared]
        for lineage in top_lineages:
            df[f"lin_{lineage[:10]}"] = (lin == lineage).astype(float).values

        # Raw correlation
        rho_raw, _ = stats.spearmanr(df["slc2a8"], df["drug"])

        # Partial correlation: regress both drug and SLC2A8 on lineage dummies,
        # then correlate residuals
        lin_cols = [c for c in df.columns if c.startswith("lin_")]
        if not lin_cols:
            continue

        from numpy.linalg import lstsq  # noqa: E402

        x_lin = df[lin_cols].values
        # Add intercept
        x_lin_int = np.column_stack([np.ones(len(x_lin)), x_lin])

        # Residualize drug
        coef_drug, _, _, _ = lstsq(x_lin_int, df["drug"].values, rcond=None)
        resid_drug = df["drug"].values - x_lin_int @ coef_drug

        # Residualize SLC2A8
        coef_slc, _, _, _ = lstsq(x_lin_int, df["slc2a8"].values, rcond=None)
        resid_slc = df["slc2a8"].values - x_lin_int @ coef_slc

        # Partial correlation
        r_partial, p_partial = stats.spearmanr(resid_slc, resid_drug)

        if r_partial < 0:
            interp = "SLC2A8 → SENSITIVITY (lineage-controlled)"
        else:
            interp = "SLC2A8 → resistance (lineage-controlled)"

        logger.info(
            f"  {drug_name:<20} {rho_raw:>+10.4f} {r_partial:>+12.4f} {p_partial:>12.4e} "
            f"{interp} {_sig(p_partial)}"
        )

    # Also check: is PDAC lineage itself predictive?
    logger.info("")
    logger.info("  PDAC LINEAGE ALONE vs DRUG SENSITIVITY:")
    is_pdac = (lineage_map == "Pancreas").astype(float)

    for drug_name in KEY_DRUGS:
        sens = _resolve_drug_sensitivity(
            drug_name, sensitivity, compound_info, name_col, id_col,
        )
        if sens is None:
            continue
        shared = sorted(set(is_pdac.index) & set(sens.index))
        if len(shared) < 50:
            continue
        pdac_sens = sens[shared][is_pdac[shared] == 1]
        other_sens = sens[shared][is_pdac[shared] == 0]
        if len(pdac_sens) < 5:
            continue
        t, p = stats.ttest_ind(pdac_sens, other_sens)
        d = (pdac_sens.mean() - other_sens.mean()) / np.sqrt(
            (pdac_sens.std() ** 2 + other_sens.std() ** 2) / 2
        )
        logger.info(
            f"  {drug_name:<20} PDAC mean={pdac_sens.mean():.4f}, "
            f"other={other_sens.mean():.4f}, d={d:+.3f}, p={p:.4e} {_sig(p)}"
        )


def main() -> None:
    """Run all three sanity checks."""
    border = "=" * 90
    logger.info(border)
    logger.info("H13e SANITY CHECKS")
    logger.info(border)

    logger.info("Loading data")
    (
        expression, metadata, scores_df,
        target_ids, background_ids, pathway_score, quadrant_labels,
    ) = _load_all()

    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    name_col = _find_compound_col(
        compound_info, ["Drug.Name", "name", "Name", "drug_name"],
    )
    id_col = _find_compound_col(
        compound_info, ["IDs", "broad_id", "column_name", "compound_id"],
    )

    logger.info("Data loaded")

    check_1_prism_orientation(sensitivity, compound_info, name_col, id_col)
    check_2_within_group_correlation(
        sensitivity, compound_info, name_col, id_col,
        expression, pathway_score, target_ids, background_ids,
    )
    check_3_lineage_controlled(
        sensitivity, compound_info, name_col, id_col,
        expression, metadata,
    )

    logger.info(border)
    logger.info("H13e SANITY CHECKS COMPLETE")
    logger.info(border)


if __name__ == "__main__":
    main()
