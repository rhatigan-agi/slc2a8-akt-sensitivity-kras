"""H13g Convergence Analysis: Do AKT drug sensitivity and lipid gene dependency
co-occur at the individual cell-line level?

Previous analyses showed that SLC2A8-high/RAS-mut cells have BOTH:
  - Higher AKT drug sensitivity (H13c/e: afuresertib partial r=-0.144)
  - Stronger SCD/FASN CRISPR dependency (H13f: SCD d=-0.316, FASN d=-0.173)

But these were group-level comparisons. This script tests whether the same
individual cell lines that are most AKT-sensitive are also the most SCD/FASN-
dependent -- the cell-level convergence that would support a combination
therapy rationale.

Analyses:
  1. Cell-line-level correlation: SCD/FASN dependency vs AKT drug sensitivity
     (pan-cancer and within SLC2A8-high/RAS-mut target group)
  2. SCD inhibitor PRISM compound-level data (which compounds, what sensitivity?)
  3. Interaction model: drug_sensitivity ~ SLC2A8 + SCD_dep + SLC2A8*SCD_dep
  4. Quadrant analysis: classify lines as AKT-sensitive AND/OR SCD-dependent
  5. Lineage-controlled convergence: does the cell-level correlation survive
     lineage correction?
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
RESULTS_DIR = Path("results/h13g")

AKT_DRUGS: list[str] = ["GSK2110183", "GDC-0068", "AZD5363"]
LIPID_DEP_GENES: list[str] = ["SCD", "FASN", "ACLY", "ACACA", "SREBF1"]


def _load_all() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame,
    list[str], list[str],
]:
    """Load expression, CRISPR, metadata; compute target/background groups."""
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
    background_ids = [i for i in common_ids if i not in target_ids]

    return expression, crispr, metadata, target_ids, background_ids


def _resolve_drug(
    drug_name: str,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> pd.Series | None:
    """Resolve drug name to sensitivity vector."""
    mask = compound_info[name_col].str.upper() == drug_name.upper()
    if mask.sum() == 0:
        return None
    for brd_id in compound_info.loc[mask, id_col].tolist():
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


def _header(title: str) -> None:
    """Print section header."""
    border = "=" * 90
    logger.info(border)
    logger.info(title)
    logger.info(border)


# =========================================================================
# ANALYSIS 1: Cell-line-level correlation of CRISPR dependency vs drug sens
# =========================================================================

def analysis_1_cellline_convergence(
    crispr: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Do individual cell lines with strong SCD/FASN dependency also show
    higher AKT drug sensitivity?"""
    _header("ANALYSIS 1: CELL-LINE-LEVEL CONVERGENCE (CRISPR dep vs drug sens)")

    logger.info("  For each (lipid gene, AKT drug) pair:")
    logger.info("  Spearman rho between CRISPR dependency and drug sensitivity")
    logger.info("  (both negative = more essential / more killing)")
    logger.info("  Positive rho = convergence (same lines are vulnerable to both)")
    logger.info("")

    for scope_name, scope_ids in [
        ("PAN-CANCER (all lines)", None),
        ("TARGET ONLY (SLC2A8-hi/RAS-mut)", target_ids),
        ("BACKGROUND ONLY", background_ids),
    ]:
        logger.info(f"  --- {scope_name} ---")
        logger.info(
            f"  {'Gene':<10} {'Drug':<20} {'rho':>8} {'p':>12} {'n':>6}"
        )
        logger.info("  " + "-" * 65)

        for gene in LIPID_DEP_GENES:
            if gene not in crispr.columns:
                continue
            for drug_name in AKT_DRUGS:
                sens = _resolve_drug(
                    drug_name, sensitivity, compound_info, name_col, id_col,
                )
                if sens is None:
                    continue

                dep = crispr[gene]
                shared = sorted(set(dep.dropna().index) & set(sens.index))
                if scope_ids is not None:
                    shared = [s for s in shared if s in scope_ids]
                if len(shared) < 10:
                    logger.info(
                        f"  {gene:<10} {drug_name:<20} {'n/a':>8} {'n/a':>12} "
                        f"{len(shared):>6}"
                    )
                    continue

                rho, p = stats.spearmanr(dep[shared], sens[shared])
                logger.info(
                    f"  {gene:<10} {drug_name:<20} {rho:>+8.4f} {p:>12.2e} "
                    f"{len(shared):>6} {_sig(p)}"
                )

        logger.info("")


# =========================================================================
# ANALYSIS 2: SCD inhibitor PRISM compound-level detail
# =========================================================================

def analysis_2_scd_inhibitor_detail(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Identify and report all SCD-related compounds in PRISM."""
    _header("ANALYSIS 2: SCD INHIBITOR PRISM DETAIL")

    # Find MOA column
    moa_col = _find_compound_col(compound_info, ["moa", "MOA", "mechanism"])
    if moa_col is None:
        # Try broader search
        for col in compound_info.columns:
            if "moa" in col.lower() or "mechanism" in col.lower():
                moa_col = col
                break

    if moa_col is None:
        logger.info("  No MOA column found in compound info. Searching by name.")
        # Fall back to name-based search
        scd_keywords = ["SCD", "stearoyl", "desaturase"]
        mask = compound_info[name_col].str.upper().apply(
            lambda x: any(kw.upper() in str(x).upper() for kw in scd_keywords)
        )
    else:
        logger.info(f"  MOA column: {moa_col}")
        scd_keywords = ["SCD", "stearoyl", "desaturase"]
        lipid_keywords = [
            "SCD", "stearoyl", "desaturase", "FASN", "fatty acid synthase",
            "lipid", "ACLY", "citrate lyase", "ACC ", "acetyl-CoA carboxylase",
        ]
        mask_scd = compound_info[moa_col].fillna("").str.upper().apply(
            lambda x: any(kw.upper() in x for kw in scd_keywords)
        )
        mask_lipid = compound_info[moa_col].fillna("").str.upper().apply(
            lambda x: any(kw.upper() in x for kw in lipid_keywords)
        )
        mask = mask_scd | mask_lipid

    scd_compounds = compound_info[mask]
    logger.info(f"  Found {len(scd_compounds)} lipid metabolism compounds")
    logger.info("")

    if len(scd_compounds) == 0:
        logger.info("  No lipid metabolism compounds found by MOA/name search")
        return

    logger.info(
        f"  {'Name':<30} {'MOA':<35} {'t_mean':>8} {'bg_mean':>8} "
        f"{'diff':>8} {'p':>12}"
    )
    logger.info("  " + "-" * 110)

    for _, row in scd_compounds.iterrows():
        brd_id = row[id_col]
        drug_name_val = row[name_col]
        moa_val = row.get(moa_col, "N/A") if moa_col else "N/A"

        if brd_id not in sensitivity.columns:
            continue

        sens = sensitivity[brd_id].dropna()
        t_sens = sens[sens.index.isin(target_ids)]
        b_sens = sens[sens.index.isin(background_ids)]

        if len(t_sens) < 5 or len(b_sens) < 5:
            continue

        _, p = stats.ttest_ind(t_sens, b_sens)
        diff = t_sens.mean() - b_sens.mean()

        logger.info(
            f"  {str(drug_name_val)[:30]:<30} {str(moa_val)[:35]:<35} "
            f"{t_sens.mean():>+8.4f} {b_sens.mean():>+8.4f} "
            f"{diff:>+8.4f} {p:>12.2e} {_sig(p)}"
        )

    logger.info("")
    logger.info("  (Negative mean = more killing. Negative diff = target more sensitive)")


# =========================================================================
# ANALYSIS 3: Interaction model
# =========================================================================

def analysis_3_interaction_model(
    expression: pd.DataFrame,
    crispr: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> None:
    """Test interaction: drug_sens ~ SLC2A8 + SCD_dep + SLC2A8 * SCD_dep."""
    _header("ANALYSIS 3: INTERACTION MODEL (SLC2A8 x SCD dependency)")

    if "SLC2A8" not in expression.columns or "SCD" not in crispr.columns:
        logger.info("  Missing SLC2A8 expression or SCD CRISPR data")
        return

    slc2a8 = expression["SLC2A8"]
    scd_dep = crispr["SCD"]

    for drug_name in AKT_DRUGS:
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            logger.info(f"  {drug_name}: NOT FOUND")
            continue

        shared = sorted(
            set(slc2a8.dropna().index)
            & set(scd_dep.dropna().index)
            & set(sens.index)
        )
        if len(shared) < 50:
            logger.info(f"  {drug_name}: insufficient overlap (n={len(shared)})")
            continue

        # Build regression DataFrame
        df = pd.DataFrame({
            "drug": sens[shared].values,
            "slc2a8": slc2a8[shared].values,
            "scd_dep": scd_dep[shared].values,
        }, index=shared)

        # Standardize predictors for interpretability
        df["slc2a8_z"] = (df["slc2a8"] - df["slc2a8"].mean()) / df["slc2a8"].std()
        df["scd_dep_z"] = (df["scd_dep"] - df["scd_dep"].mean()) / df["scd_dep"].std()
        df["interaction"] = df["slc2a8_z"] * df["scd_dep_z"]

        # OLS regression
        from numpy.linalg import lstsq

        x = np.column_stack([
            np.ones(len(df)),
            df["slc2a8_z"].values,
            df["scd_dep_z"].values,
            df["interaction"].values,
        ])
        y = df["drug"].values

        coef, residuals, rank, sv = lstsq(x, y, rcond=None)
        y_hat = x @ coef
        resid = y - y_hat
        n_obs = len(y)
        n_params = x.shape[1]
        mse = np.sum(resid ** 2) / (n_obs - n_params)
        se = np.sqrt(np.diag(mse * np.linalg.inv(x.T @ x)))
        t_vals = coef / se
        p_vals = 2 * stats.t.sf(np.abs(t_vals), df=n_obs - n_params)

        r2 = 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)

        logger.info(f"  {drug_name} (n={n_obs}, R2={r2:.4f}):")
        labels = ["intercept", "SLC2A8 (z)", "SCD_dep (z)", "SLC2A8 x SCD_dep"]
        for label, c, s, t, p in zip(labels, coef, se, t_vals, p_vals):
            logger.info(
                f"    {label:<25} coef={c:>+8.5f}  se={s:.5f}  "
                f"t={t:>+7.3f}  p={p:.4e} {_sig(p)}"
            )
        logger.info("")

    # Repeat with FASN
    logger.info("  --- Same model with FASN dependency ---")
    if "FASN" not in crispr.columns:
        logger.info("  FASN not in CRISPR data")
        return

    fasn_dep = crispr["FASN"]

    for drug_name in AKT_DRUGS:
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            continue

        shared = sorted(
            set(slc2a8.dropna().index)
            & set(fasn_dep.dropna().index)
            & set(sens.index)
        )
        if len(shared) < 50:
            continue

        df = pd.DataFrame({
            "drug": sens[shared].values,
            "slc2a8": slc2a8[shared].values,
            "fasn_dep": fasn_dep[shared].values,
        }, index=shared)

        df["slc2a8_z"] = (df["slc2a8"] - df["slc2a8"].mean()) / df["slc2a8"].std()
        df["fasn_dep_z"] = (df["fasn_dep"] - df["fasn_dep"].mean()) / df["fasn_dep"].std()
        df["interaction"] = df["slc2a8_z"] * df["fasn_dep_z"]

        from numpy.linalg import lstsq

        x = np.column_stack([
            np.ones(len(df)),
            df["slc2a8_z"].values,
            df["fasn_dep_z"].values,
            df["interaction"].values,
        ])
        y = df["drug"].values

        coef, _, _, _ = lstsq(x, y, rcond=None)
        y_hat = x @ coef
        resid = y - y_hat
        n_obs = len(y)
        n_params = x.shape[1]
        mse = np.sum(resid ** 2) / (n_obs - n_params)
        se = np.sqrt(np.diag(mse * np.linalg.inv(x.T @ x)))
        t_vals = coef / se
        p_vals = 2 * stats.t.sf(np.abs(t_vals), df=n_obs - n_params)

        r2 = 1 - np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2)

        logger.info(f"  {drug_name} (n={n_obs}, R2={r2:.4f}):")
        labels = ["intercept", "SLC2A8 (z)", "FASN_dep (z)", "SLC2A8 x FASN_dep"]
        for label, c, s, t, p in zip(labels, coef, se, t_vals, p_vals):
            logger.info(
                f"    {label:<25} coef={c:>+8.5f}  se={s:.5f}  "
                f"t={t:>+7.3f}  p={p:.4e} {_sig(p)}"
            )
        logger.info("")


# =========================================================================
# ANALYSIS 4: Quadrant classification (AKT-sensitive AND/OR SCD-dependent)
# =========================================================================

def analysis_4_quadrant_classification(
    crispr: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Classify lines into quadrants of AKT sensitivity x SCD dependency."""
    _header("ANALYSIS 4: QUADRANT CLASSIFICATION (AKT sens x SCD dep)")

    if "SCD" not in crispr.columns:
        logger.info("  SCD not in CRISPR data")
        return

    # Use afuresertib (strongest AKT signal)
    sens = _resolve_drug("GSK2110183", sensitivity, compound_info, name_col, id_col)
    if sens is None:
        logger.info("  GSK2110183 (afuresertib) not found")
        return

    scd_dep = crispr["SCD"]
    shared = sorted(set(scd_dep.dropna().index) & set(sens.index))

    df = pd.DataFrame({
        "akt_sens": sens[shared],
        "scd_dep": scd_dep[shared],
    })
    df["is_target"] = df.index.isin(target_ids)

    # Define thresholds: median split
    akt_thresh = df["akt_sens"].median()
    scd_thresh = df["scd_dep"].median()

    logger.info(f"  Thresholds: AKT sensitivity median={akt_thresh:.4f}, "
                f"SCD dependency median={scd_thresh:.4f}")
    logger.info(f"  (AKT: lower=more sensitive; SCD: lower=more dependent)")
    logger.info("")

    df["akt_sensitive"] = df["akt_sens"] < akt_thresh
    df["scd_dependent"] = df["scd_dep"] < scd_thresh

    # 4 quadrants
    q_labels = {
        (True, True): "AKT-sens + SCD-dep (DUAL VULN)",
        (True, False): "AKT-sens only",
        (False, True): "SCD-dep only",
        (False, False): "Neither",
    }

    logger.info(
        f"  {'Quadrant':<35} {'Total':>8} {'Target':>8} {'Target%':>8} "
        f"{'OR':>8} {'p':>12}"
    )
    logger.info("  " + "-" * 85)

    # Build contingency data for enrichment
    total_target = df["is_target"].sum()
    total_bg = (~df["is_target"]).sum()

    for (akt_s, scd_d), label in q_labels.items():
        in_q = df["akt_sensitive"].eq(akt_s) & df["scd_dependent"].eq(scd_d)
        n_total = in_q.sum()
        n_target = (in_q & df["is_target"]).sum()
        pct_target = 100 * n_target / n_total if n_total > 0 else 0

        # Fisher exact test: is target enriched in this quadrant?
        n_bg_in_q = n_total - n_target
        n_target_out = total_target - n_target
        n_bg_out = total_bg - n_bg_in_q

        contingency = [[n_target, n_bg_in_q], [n_target_out, n_bg_out]]
        odds_ratio, p_fisher = stats.fisher_exact(contingency)

        logger.info(
            f"  {label:<35} {n_total:>8} {n_target:>8} {pct_target:>7.1f}% "
            f"{odds_ratio:>8.2f} {p_fisher:>12.2e} {_sig(p_fisher)}"
        )

    logger.info("")

    # Same analysis with stricter thresholds (25th percentile)
    logger.info("  --- Stricter thresholds (25th percentile) ---")
    akt_q25 = df["akt_sens"].quantile(0.25)
    scd_q25 = df["scd_dep"].quantile(0.25)
    logger.info(f"  Thresholds: AKT 25th pctl={akt_q25:.4f}, SCD 25th pctl={scd_q25:.4f}")

    df["akt_sens_strict"] = df["akt_sens"] < akt_q25
    df["scd_dep_strict"] = df["scd_dep"] < scd_q25

    q_labels_strict = {
        (True, True): "AKT-sens + SCD-dep (DUAL VULN)",
        (True, False): "AKT-sens only",
        (False, True): "SCD-dep only",
        (False, False): "Neither",
    }

    logger.info(
        f"  {'Quadrant':<35} {'Total':>8} {'Target':>8} {'Target%':>8} "
        f"{'OR':>8} {'p':>12}"
    )
    logger.info("  " + "-" * 85)

    for (akt_s, scd_d), label in q_labels_strict.items():
        in_q = df["akt_sens_strict"].eq(akt_s) & df["scd_dep_strict"].eq(scd_d)
        n_total = in_q.sum()
        n_target = (in_q & df["is_target"]).sum()
        pct_target = 100 * n_target / n_total if n_total > 0 else 0

        n_bg_in_q = n_total - n_target
        n_target_out = total_target - n_target
        n_bg_out = total_bg - n_bg_in_q

        contingency = [[n_target, n_bg_in_q], [n_target_out, n_bg_out]]
        odds_ratio, p_fisher = stats.fisher_exact(contingency)

        logger.info(
            f"  {label:<35} {n_total:>8} {n_target:>8} {pct_target:>7.1f}% "
            f"{odds_ratio:>8.2f} {p_fisher:>12.2e} {_sig(p_fisher)}"
        )


# =========================================================================
# ANALYSIS 5: Lineage-controlled convergence
# =========================================================================

def analysis_5_lineage_controlled_convergence(
    expression: pd.DataFrame,
    crispr: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    metadata: pd.DataFrame,
) -> None:
    """Does the SCD dep <-> AKT sensitivity correlation survive lineage control?"""
    _header("ANALYSIS 5: LINEAGE-CONTROLLED CONVERGENCE")

    if "SCD" not in crispr.columns:
        logger.info("  SCD not in CRISPR data")
        return

    from numpy.linalg import lstsq

    lineage_col = "OncotreeLineage"
    lineage_map = metadata[lineage_col].reindex(crispr.index).dropna()
    top_lineages = lineage_map.value_counts().head(15).index.tolist()

    logger.info(f"  Controlling for top 15 lineages")
    logger.info("")

    for dep_gene in ["SCD", "FASN"]:
        if dep_gene not in crispr.columns:
            continue

        logger.info(f"  --- {dep_gene} dependency vs AKT drug sensitivity ---")
        logger.info(
            f"  {'Drug':<20} {'Raw rho':>10} {'Partial r':>12} {'Partial p':>12} "
            f"{'n':>6}"
        )
        logger.info("  " + "-" * 70)

        dep = crispr[dep_gene]

        for drug_name in AKT_DRUGS:
            sens = _resolve_drug(
                drug_name, sensitivity, compound_info, name_col, id_col,
            )
            if sens is None:
                continue

            shared = sorted(
                set(dep.dropna().index) & set(sens.index) & set(lineage_map.index)
            )
            if len(shared) < 50:
                continue

            # Raw correlation
            rho_raw, _ = stats.spearmanr(dep[shared], sens[shared])

            # Lineage residualization
            df_reg = pd.DataFrame({
                "dep": dep[shared].values,
                "drug": sens[shared].values,
            }, index=shared)

            lin = lineage_map[shared]
            for lineage in top_lineages:
                df_reg[f"lin_{lineage[:10]}"] = (lin == lineage).astype(float).values

            lin_cols = [c for c in df_reg.columns if c.startswith("lin_")]
            x_lin = np.column_stack([
                np.ones(len(df_reg)),
                df_reg[lin_cols].values,
            ])

            coef_dep, _, _, _ = lstsq(x_lin, df_reg["dep"].values, rcond=None)
            resid_dep = df_reg["dep"].values - x_lin @ coef_dep

            coef_drug, _, _, _ = lstsq(x_lin, df_reg["drug"].values, rcond=None)
            resid_drug = df_reg["drug"].values - x_lin @ coef_drug

            r_partial, p_partial = stats.spearmanr(resid_dep, resid_drug)

            logger.info(
                f"  {drug_name:<20} {rho_raw:>+10.4f} {r_partial:>+12.4f} "
                f"{p_partial:>12.2e} {len(shared):>6} {_sig(p_partial)}"
            )

        logger.info("")

    # Also: SLC2A8 expression vs SCD dependency (lineage-controlled)
    logger.info("  --- SLC2A8 expression vs lipid gene dependency (lineage-controlled) ---")
    if "SLC2A8" not in expression.columns:
        logger.info("  SLC2A8 not in expression")
        return

    slc2a8 = expression["SLC2A8"]

    logger.info(
        f"  {'Gene':<10} {'Raw rho':>10} {'Partial r':>12} {'Partial p':>12} {'n':>6}"
    )
    logger.info("  " + "-" * 60)

    for dep_gene in LIPID_DEP_GENES:
        if dep_gene not in crispr.columns:
            continue

        dep = crispr[dep_gene]
        shared = sorted(
            set(slc2a8.dropna().index) & set(dep.dropna().index)
            & set(lineage_map.index)
        )
        if len(shared) < 50:
            continue

        rho_raw, _ = stats.spearmanr(slc2a8[shared], dep[shared])

        df_reg = pd.DataFrame({
            "slc2a8": slc2a8[shared].values,
            "dep": dep[shared].values,
        }, index=shared)

        lin = lineage_map[shared]
        for lineage in top_lineages:
            df_reg[f"lin_{lineage[:10]}"] = (lin == lineage).astype(float).values

        lin_cols = [c for c in df_reg.columns if c.startswith("lin_")]
        x_lin = np.column_stack([np.ones(len(df_reg)), df_reg[lin_cols].values])

        coef_s, _, _, _ = lstsq(x_lin, df_reg["slc2a8"].values, rcond=None)
        resid_s = df_reg["slc2a8"].values - x_lin @ coef_s

        coef_d, _, _, _ = lstsq(x_lin, df_reg["dep"].values, rcond=None)
        resid_d = df_reg["dep"].values - x_lin @ coef_d

        r_partial, p_partial = stats.spearmanr(resid_s, resid_d)

        logger.info(
            f"  {dep_gene:<10} {rho_raw:>+10.4f} {r_partial:>+12.4f} "
            f"{p_partial:>12.2e} {len(shared):>6} {_sig(p_partial)}"
        )


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    """Run H13g convergence analysis."""
    _header("H13g: CONVERGENCE ANALYSIS — AKT DRUG SENSITIVITY x LIPID DEPENDENCY")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading data")
    expression, crispr, metadata, target_ids, background_ids = _load_all()

    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)
    name_col = _find_compound_col(
        compound_info, ["Drug.Name", "name", "Name"],
    )
    id_col = _find_compound_col(
        compound_info, ["IDs", "broad_id", "column_name"],
    )
    logger.info(
        f"Data loaded: {len(expression)} cell lines, "
        f"{len(target_ids)} target, {len(background_ids)} background"
    )

    analysis_1_cellline_convergence(
        crispr, sensitivity, compound_info, name_col, id_col,
        target_ids, background_ids,
    )
    analysis_2_scd_inhibitor_detail(
        sensitivity, compound_info, name_col, id_col,
        target_ids, background_ids,
    )
    analysis_3_interaction_model(
        expression, crispr, sensitivity, compound_info, name_col, id_col,
    )
    analysis_4_quadrant_classification(
        crispr, sensitivity, compound_info, name_col, id_col,
        target_ids, background_ids,
    )
    analysis_5_lineage_controlled_convergence(
        expression, crispr, sensitivity, compound_info, name_col, id_col,
        metadata,
    )

    _header("H13g CONVERGENCE ANALYSIS COMPLETE")


if __name__ == "__main__":
    main()
