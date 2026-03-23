"""H13i: PIK3CA/PTEN mutation confound check for AKT inhibitor sensitivity signal.

Davies et al. (2012) showed PIK3CA/PTEN mutations predict AZD5363 sensitivity and
RAS mutations predict resistance (pan-cancer). This script tests whether H13's
SLC2A8 → AKT drug sensitivity signal survives after adding PIK3CA and PTEN mutation
status as covariates.

Three models per drug:
  M1: drug ~ SLC2A8 + lineage              (existing H13e result)
  M2: drug ~ SLC2A8 + lineage + PIK3CA_mut (primary confound check)
  M3: drug ~ SLC2A8 + lineage + PIK3CA_mut + PTEN_mut (full check)

Also checks: collinearity between SLC2A8 expression and PIK3CA/PTEN mutation status.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pdac.h13_trehalose_vulnerability.data import (
    _find_column,
    load_expression,
    load_model_metadata,
    load_mutations,
)
from pdac.h13_trehalose_vulnerability.h13b_prism import (
    _find_compound_col,
    load_prism_compound_info,
    load_prism_sensitivity,
)
from pdac.shared.logging import get_logger

logger = get_logger(__name__)

DEPMAP_DIR = Path("data/raw/depmap")
RESULTS_DIR = Path("results/h13i")

AKT_DRUGS: dict[str, str] = {
    "Capivasertib": "AZD5363",
    "Afuresertib": "GSK2110183",
    "Ipatasertib": "GDC-0068",
}

TOP_LINEAGES = 15


def _get_mutant_ids(mutations: pd.DataFrame, gene: str) -> set[str]:
    """Return set of ModelIDs with any somatic mutation in the given gene."""
    hugo_col = _find_column(mutations, ["HugoSymbol", "Hugo_Symbol"])
    model_col = _find_column(mutations, ["ModelID", "DepMap_ID"])
    mask = mutations[hugo_col] == gene
    return set(mutations.loc[mask, model_col].unique())


def _residualize(x: np.ndarray, confounders: np.ndarray) -> np.ndarray:
    """Return residuals of x after projecting out confounders (with intercept)."""
    x_int = np.column_stack([np.ones(len(confounders)), confounders])
    coef, _, _, _ = np.linalg.lstsq(x_int, x, rcond=None)
    return x - x_int @ coef


def _ols_summary(
    y: np.ndarray,
    x_mat: np.ndarray,
    param_names: list[str],
) -> dict[str, dict[str, float]]:
    """OLS fit; return dict of {name: {beta, se, t, p}} for each param."""
    beta, _, _, _ = np.linalg.lstsq(x_mat, y, rcond=None)
    y_hat = x_mat @ beta
    n, k = x_mat.shape
    ss_res = float(np.sum((y - y_hat) ** 2))
    mse = ss_res / (n - k)
    cov = mse * np.linalg.inv(x_mat.T @ x_mat)
    se = np.sqrt(np.diag(cov))
    t_stat = beta / se
    p_vals = 2 * stats.t.sf(np.abs(t_stat), df=n - k)
    return {
        name: {"beta": float(b), "se": float(s), "t": float(t), "p": float(p)}
        for name, b, s, t, p in zip(param_names, beta, se, t_stat, p_vals)
    }


def _resolve_drug(
    name: str,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> pd.Series | None:
    mask = compound_info[name_col].str.upper() == name.upper()
    if mask.sum() == 0:
        return None
    for brd_id in compound_info.loc[mask, id_col].tolist():
        if brd_id in sensitivity.columns:
            return sensitivity[brd_id].dropna()
    return None


def _sig(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.1:
        return "."
    return ""


def check_collinearity(
    expression: pd.DataFrame,
    pik3ca_ids: set[str],
    pten_ids: set[str],
) -> None:
    """Check correlation between SLC2A8 expression and PIK3CA/PTEN mutation status."""
    logger.info("=" * 70)
    logger.info("COLLINEARITY CHECK: SLC2A8 vs PIK3CA/PTEN mutation status")
    logger.info("=" * 70)

    slc2a8 = expression["SLC2A8"].dropna()
    shared = slc2a8.index.tolist()

    pik3ca_vec = np.array([1.0 if s in pik3ca_ids else 0.0 for s in shared])
    pten_vec = np.array([1.0 if s in pten_ids else 0.0 for s in shared])

    rho_pik3ca, p_pik3ca = stats.spearmanr(slc2a8.values, pik3ca_vec)
    rho_pten, p_pten = stats.spearmanr(slc2a8.values, pten_vec)

    n_pik3ca = int(pik3ca_vec.sum())
    n_pten = int(pten_vec.sum())
    n_total = len(shared)

    logger.info(
        "  PIK3CA-mutant: %d / %d lines (%.1f%%)",
        n_pik3ca, n_total, 100 * n_pik3ca / n_total,
    )
    logger.info(
        "  PTEN-mutant:   %d / %d lines (%.1f%%)",
        n_pten, n_total, 100 * n_pten / n_total,
    )
    logger.info("")
    logger.info("  SLC2A8 expression vs PIK3CA_mut: rho=%.3f, p=%.4f %s",
                rho_pik3ca, p_pik3ca, _sig(p_pik3ca))
    logger.info("  SLC2A8 expression vs PTEN_mut:   rho=%.3f, p=%.4f %s",
                rho_pten, p_pten, _sig(p_pten))

    # Point-biserial: mean SLC2A8 in mutant vs wildtype
    for gene, vec, n_mut in [("PIK3CA", pik3ca_vec, n_pik3ca), ("PTEN", pten_vec, n_pten)]:
        mut_expr = slc2a8.values[vec == 1]
        wt_expr = slc2a8.values[vec == 0]
        if n_mut >= 5:
            t, p = stats.ttest_ind(mut_expr, wt_expr)
            d = (mut_expr.mean() - wt_expr.mean()) / np.sqrt(
                (mut_expr.std() ** 2 + wt_expr.std() ** 2) / 2
            )
            logger.info(
                "  SLC2A8 in %s-mut: mean=%.3f vs WT: mean=%.3f  d=%+.3f  p=%.4f %s",
                gene, mut_expr.mean(), wt_expr.mean(), d, p, _sig(p),
            )


def check_pik3ca_confound(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    pik3ca_ids: set[str],
    pten_ids: set[str],
) -> None:
    """For each AKT drug, fit M1/M2/M3 and report SLC2A8 beta/p across models."""
    logger.info("")
    logger.info("=" * 70)
    logger.info("CONFOUND CHECK: drug ~ SLC2A8 + lineage [+ PIK3CA_mut] [+ PTEN_mut]")
    logger.info("=" * 70)

    lineage_col = "OncotreeLineage"
    lineage_map = metadata[lineage_col].reindex(expression.index).dropna()
    top_lineages = lineage_map.value_counts().head(TOP_LINEAGES).index.tolist()

    slc2a8 = expression["SLC2A8"].dropna()

    for label, drug_name in AKT_DRUGS.items():
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            logger.info("  %s: NOT FOUND in PRISM", label)
            continue

        shared = sorted(
            set(slc2a8.index) & set(sens.index) & set(lineage_map.index)
        )
        if len(shared) < 50:
            logger.info("  %s: insufficient overlap (n=%d)", label, len(shared))
            continue

        n = len(shared)
        slc2a8_vals = slc2a8[shared].values
        drug_vals = sens[shared].values
        lin = lineage_map[shared]
        pik3ca_vec = np.array([1.0 if s in pik3ca_ids else 0.0 for s in shared])
        pten_vec = np.array([1.0 if s in pten_ids else 0.0 for s in shared])

        # Build lineage dummy matrix
        lin_dummies = np.column_stack([
            (lin == l).astype(float).values for l in top_lineages
        ])

        # -------------------------------------------------------
        # Partial correlation approach (matches H13e methodology)
        # -------------------------------------------------------
        def partial_r_slc2a8(extra_covariates: np.ndarray | None = None) -> tuple[float, float, int]:
            """Residualize drug and SLC2A8 on lineage [+ extra_covariates]."""
            if extra_covariates is not None:
                conf = np.column_stack([lin_dummies, extra_covariates])
            else:
                conf = lin_dummies
            resid_drug = _residualize(drug_vals, conf)
            resid_slc = _residualize(slc2a8_vals, conf)
            r, p = stats.spearmanr(resid_slc, resid_drug)
            return float(r), float(p), n

        r_m1, p_m1, _ = partial_r_slc2a8(None)
        r_m2, p_m2, _ = partial_r_slc2a8(pik3ca_vec.reshape(-1, 1))
        r_m3, p_m3, _ = partial_r_slc2a8(
            np.column_stack([pik3ca_vec, pten_vec])
        )

        # Raw Spearman for reference
        rho_raw, p_raw = stats.spearmanr(slc2a8_vals, drug_vals)

        logger.info("")
        logger.info("  %s (n=%d):", label, n)
        logger.info("    Raw Spearman:          rho=%+.3f  p=%.4f %s",
                    rho_raw, p_raw, _sig(p_raw))
        logger.info("    M1 (+lineage):         r=%+.3f  p=%.4f %s  [H13e result]",
                    r_m1, p_m1, _sig(p_m1))
        logger.info("    M2 (+lineage+PIK3CA):  r=%+.3f  p=%.4f %s",
                    r_m2, p_m2, _sig(p_m2))
        logger.info("    M3 (+lineage+PIK3CA+PTEN): r=%+.3f  p=%.4f %s",
                    r_m3, p_m3, _sig(p_m3))

        # Interpretation
        survived = p_m3 < 0.05 and r_m3 < 0
        attenuation = abs(r_m1 - r_m3) / (abs(r_m1) + 1e-9)
        logger.info(
            "    Attenuation M1->M3: %.1f%%  |  Signal survives M3: %s",
            attenuation * 100,
            "YES" if survived else "NO",
        )

        # Also: are PIK3CA/PTEN mutants enriched in SLC2A8-high cells?
        slc2a8_median = np.median(slc2a8_vals)
        slc2a8_high = set(np.array(shared)[slc2a8_vals >= slc2a8_median])
        pik3ca_in_high = sum(1 for s in shared if s in pik3ca_ids and s in slc2a8_high)
        pik3ca_total = sum(1 for s in shared if s in pik3ca_ids)
        pten_in_high = sum(1 for s in shared if s in pten_ids and s in slc2a8_high)
        pten_total = sum(1 for s in shared if s in pten_ids)

        if pik3ca_total > 0:
            logger.info(
                "    PIK3CA-mut in SLC2A8-high half: %d / %d (%.0f%%)",
                pik3ca_in_high, pik3ca_total, 100 * pik3ca_in_high / pik3ca_total,
            )
        if pten_total > 0:
            logger.info(
                "    PTEN-mut in SLC2A8-high half:   %d / %d (%.0f%%)",
                pten_in_high, pten_total, 100 * pten_in_high / pten_total,
            )


def main() -> None:
    """Run PIK3CA/PTEN confound check."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("H13i: PIK3CA/PTEN confound check")
    logger.info("=" * 70)

    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)
    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])

    pik3ca_ids = _get_mutant_ids(mutations, "PIK3CA")
    pten_ids = _get_mutant_ids(mutations, "PTEN")

    logger.info("PIK3CA-mutant lines found: %d", len(pik3ca_ids))
    logger.info("PTEN-mutant lines found:   %d", len(pten_ids))

    check_collinearity(expression, pik3ca_ids, pten_ids)
    check_pik3ca_confound(
        expression, metadata, sensitivity, compound_info,
        name_col, id_col, pik3ca_ids, pten_ids,
    )

    logger.info("")
    logger.info("=" * 70)
    logger.info("H13i complete. Results in %s", RESULTS_DIR)


if __name__ == "__main__":
    main()
