"""H13d PI3K/AKT Pathway Activity Validation.

Tests whether SLC2A8-high / RAS-mutant cells show elevated PI3K/AKT pathway
expression, which would explain their preferential sensitivity to PI3K/AKT
inhibitors observed in H13c.

Analyses:
  1. PI3K/AKT pathway gene expression (target vs background)
  2. Correlation of individual pathway genes with vulnerability score
  3. PDAC-specific vs pan-cancer comparison
  4. PI3K pathway score vs drug sensitivity (closing the loop)
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from pdac.h13_trehalose_vulnerability.data import (
    get_pdac_model_ids,
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
RESULTS_DIR = Path("results/h13d")

# ---------------------------------------------------------------------------
# PI3K/AKT pathway genes to test
# ---------------------------------------------------------------------------

# Core signalling cascade members
PI3K_AKT_CORE: dict[str, str] = {
    "PIK3CA": "PI3K catalytic subunit alpha (class I)",
    "PIK3CB": "PI3K catalytic subunit beta",
    "PIK3CD": "PI3K catalytic subunit delta",
    "PIK3R1": "PI3K regulatory subunit 1",
    "PIK3R2": "PI3K regulatory subunit 2",
    "AKT1": "AKT serine/threonine kinase 1",
    "AKT2": "AKT serine/threonine kinase 2",
    "AKT3": "AKT serine/threonine kinase 3",
    "MTOR": "mechanistic target of rapamycin",
    "RPTOR": "regulatory associated protein of mTOR (mTORC1)",
    "RICTOR": "RPTOR independent companion of mTOR (mTORC2)",
    "RPS6KB1": "ribosomal protein S6 kinase B1 (p70S6K)",
    "RPS6KB2": "ribosomal protein S6 kinase B2",
    "EIF4EBP1": "eukaryotic translation initiation factor 4E binding protein 1",
    "PTEN": "phosphatase and tensin homolog (negative regulator)",
}

# Downstream readouts / activity markers
PI3K_AKT_READOUTS: dict[str, str] = {
    "FOXO3": "forkhead box O3 (AKT substrate)",
    "GSK3B": "glycogen synthase kinase 3 beta (AKT substrate)",
    "BAD": "BCL2 associated agonist of cell death (AKT substrate)",
    "TSC2": "TSC complex subunit 2 / tuberin (AKT → mTOR relay)",
    "SGK1": "serum/glucocorticoid regulated kinase 1 (PI3K target)",
    "PDPK1": "3-phosphoinositide dependent protein kinase 1 (PDK1)",
}

# Combine for full analysis
ALL_PI3K_GENES: dict[str, str] = {**PI3K_AKT_CORE, **PI3K_AKT_READOUTS}


def _load_data() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame,
    list[str], list[str], pd.DataFrame,
]:
    """Load DepMap data and compute target/background groups."""
    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    crispr = load_crispr_dependencies(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)

    gene_sets = GeneSetCollection()
    ras_mutant_ids = get_ras_mutant_model_ids(mutations, gene_sets.ras_mutations)

    common_ids = sorted(
        set(expression.index) & set(crispr.index) & set(metadata.index)
    )
    expression_common = expression.loc[common_ids]
    crispr_common = crispr.loc[common_ids]

    scores_df = compute_vulnerability_score(
        expression_common, gene_sets.target, gene_sets.scoring_sets(),
    )
    quadrant_labels = stratify_by_target_and_ras(scores_df, ras_mutant_ids)

    target_ids = quadrant_labels[
        quadrant_labels == "SLC2A8-high / RAS-mut"
    ].index.tolist()
    background_ids = [i for i in common_ids if i not in target_ids]

    return expression_common, crispr_common, scores_df, target_ids, background_ids, metadata


# =========================================================================
# ANALYSIS 1: PI3K/AKT pathway expression — target vs background
# =========================================================================


def analysis_1_pathway_expression(
    expression: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> pd.DataFrame:
    """Compare PI3K/AKT gene expression between target and background."""
    header = "ANALYSIS 1: PI3K/AKT PATHWAY EXPRESSION (TARGET vs BACKGROUND)"
    logger.info(header)
    _print_header(header)

    _print(f"  Target group (SLC2A8-high / RAS-mut): {len(target_ids)}")
    _print(f"  Background: {len(background_ids)}")
    _print("")

    results = []
    for gene, desc in ALL_PI3K_GENES.items():
        if gene not in expression.columns:
            _print(f"  {gene}: NOT IN EXPRESSION MATRIX")
            continue

        target_vals = expression.loc[
            expression.index.isin(target_ids), gene
        ].dropna()
        bg_vals = expression.loc[
            expression.index.isin(background_ids), gene
        ].dropna()

        if len(target_vals) < 5 or len(bg_vals) < 5:
            continue

        t_stat, p_val = stats.ttest_ind(target_vals, bg_vals, alternative="two-sided")
        cohens_d = (target_vals.mean() - bg_vals.mean()) / np.sqrt(
            (target_vals.std() ** 2 + bg_vals.std() ** 2) / 2
        )

        # Direction: positive d means target has HIGHER expression
        direction = "HIGHER" if cohens_d > 0 else "lower"

        results.append({
            "gene": gene,
            "description": desc,
            "target_mean": target_vals.mean(),
            "bg_mean": bg_vals.mean(),
            "diff": target_vals.mean() - bg_vals.mean(),
            "cohens_d": cohens_d,
            "p_value": p_val,
            "direction": direction,
            "target_n": len(target_vals),
            "bg_n": len(bg_vals),
        })

    results_df = pd.DataFrame(results).sort_values("p_value")

    _print(f"  {'Gene':<12} {'Description':<50} {'d':>8} {'p-value':>12} {'Dir':>8}")
    _print("  " + "-" * 100)
    for _, row in results_df.iterrows():
        sig = _sig_stars(row["p_value"])
        _print(
            f"  {row['gene']:<12} {row['description']:<50} "
            f"{row['cohens_d']:>8.3f} {row['p_value']:>12.4e} {row['direction']:>8} {sig}"
        )

    # Composite PI3K pathway z-score
    _print("")
    _print("  COMPOSITE PI3K/AKT PATHWAY Z-SCORE:")
    pi3k_genes_available = [g for g in PI3K_AKT_CORE if g in expression.columns and g != "PTEN"]
    if pi3k_genes_available:
        subset = expression[pi3k_genes_available]
        zscored = subset.apply(stats.zscore, nan_policy="omit")
        pathway_score = zscored.mean(axis=1)

        target_scores = pathway_score[pathway_score.index.isin(target_ids)].dropna()
        bg_scores = pathway_score[pathway_score.index.isin(background_ids)].dropna()
        t_stat, p_val = stats.ttest_ind(target_scores, bg_scores)
        d = (target_scores.mean() - bg_scores.mean()) / np.sqrt(
            (target_scores.std() ** 2 + bg_scores.std() ** 2) / 2
        )
        _print(f"    Genes used: {pi3k_genes_available}")
        _print(f"    Target mean z: {target_scores.mean():.4f}")
        _print(f"    Background mean z: {bg_scores.mean():.4f}")
        _print(f"    Cohen's d: {d:.4f}")
        _print(f"    p-value: {p_val:.4e} {_sig_stars(p_val)}")
        _print(f"    Direction: {'TARGET HIGHER' if d > 0 else 'target lower'}")

    return results_df


# =========================================================================
# ANALYSIS 2: Correlation with vulnerability score
# =========================================================================


def analysis_2_correlation(
    expression: pd.DataFrame,
    scores_df: pd.DataFrame,
) -> pd.DataFrame:
    """Correlate PI3K/AKT gene expression with vulnerability score."""
    header = "ANALYSIS 2: PI3K/AKT GENE EXPRESSION vs VULNERABILITY SCORE"
    _print_header(header)

    common = sorted(set(expression.index) & set(scores_df.index))
    vuln = scores_df.loc[common, "vulnerability_score"]

    results = []
    for gene in ALL_PI3K_GENES:
        if gene not in expression.columns:
            continue
        vals = expression.loc[common, gene].dropna()
        shared = sorted(set(vals.index) & set(vuln.dropna().index))
        if len(shared) < 30:
            continue

        rho, p_val = stats.spearmanr(vals[shared], vuln[shared])
        results.append({
            "gene": gene,
            "spearman_rho": rho,
            "p_value": p_val,
            "n": len(shared),
        })

    results_df = pd.DataFrame(results).sort_values("p_value")

    _print(f"  {'Gene':<12} {'Spearman rho':>14} {'p-value':>12} {'n':>6}")
    _print("  " + "-" * 60)
    for _, row in results_df.iterrows():
        sig = _sig_stars(row["p_value"])
        _print(
            f"  {row['gene']:<12} {row['spearman_rho']:>14.4f} "
            f"{row['p_value']:>12.4e} {row['n']:>6} {sig}"
        )

    # SLC2A8 alone vs pathway genes
    _print("")
    _print("  SLC2A8 EXPRESSION vs PI3K PATHWAY GENES (direct correlation):")
    if "SLC2A8" in expression.columns:
        slc2a8 = expression.loc[common, "SLC2A8"].dropna()
        for gene in ["AKT1", "AKT2", "PIK3CA", "PIK3CB", "MTOR", "RPS6KB1", "PTEN"]:
            if gene not in expression.columns:
                continue
            gvals = expression.loc[common, gene].dropna()
            shared = sorted(set(slc2a8.index) & set(gvals.index))
            if len(shared) < 30:
                continue
            rho, p_val = stats.spearmanr(slc2a8[shared], gvals[shared])
            sig = _sig_stars(p_val)
            _print(f"    SLC2A8 vs {gene:<10}: rho={rho:+.4f}, p={p_val:.4e} {sig}")

    return results_df


# =========================================================================
# ANALYSIS 3: PDAC-specific pathway activity check
# =========================================================================


def analysis_3_pdac_specific(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Check pathway activity specifically within PDAC cell lines."""
    header = "ANALYSIS 3: PDAC-SPECIFIC PI3K/AKT PATHWAY CHECK"
    _print_header(header)

    pdac_ids = get_pdac_model_ids(metadata)
    pdac_target = [i for i in target_ids if i in pdac_ids]
    pdac_bg = [i for i in background_ids if i in pdac_ids]

    _print(f"  PDAC target (SLC2A8-high / RAS-mut): {len(pdac_target)}")
    _print(f"  PDAC background: {len(pdac_bg)}")
    _print("")

    key_genes = ["AKT1", "AKT2", "PIK3CA", "PIK3CB", "MTOR", "RPS6KB1",
                 "PTEN", "PDPK1", "RICTOR", "RPTOR"]

    for gene in key_genes:
        if gene not in expression.columns:
            continue

        t_vals = expression.loc[
            expression.index.isin(pdac_target), gene
        ].dropna()
        b_vals = expression.loc[
            expression.index.isin(pdac_bg), gene
        ].dropna()

        if len(t_vals) < 3 or len(b_vals) < 3:
            _print(f"  {gene}: insufficient PDAC samples")
            continue

        t_stat, p_val = stats.ttest_ind(t_vals, b_vals)
        d = (t_vals.mean() - b_vals.mean()) / np.sqrt(
            (t_vals.std() ** 2 + b_vals.std() ** 2) / 2
        )
        direction = "HIGHER" if d > 0 else "lower"
        sig = _sig_stars(p_val)
        _print(
            f"  {gene:<10}: target={t_vals.mean():.3f}, bg={b_vals.mean():.3f}, "
            f"d={d:+.3f} ({direction}), p={p_val:.4e} {sig}"
        )


# =========================================================================
# ANALYSIS 4: PI3K pathway score vs drug sensitivity (causal link)
# =========================================================================


def analysis_4_pathway_vs_drugs(
    expression: pd.DataFrame,
    scores_df: pd.DataFrame,
) -> None:
    """Test whether PI3K pathway activity predicts drug sensitivity directly."""
    header = "ANALYSIS 4: PI3K PATHWAY SCORE → DRUG SENSITIVITY (CAUSAL LINK)"
    _print_header(header)

    # Compute PI3K pathway z-score (excluding PTEN, which is a negative regulator)
    pi3k_genes = [g for g in PI3K_AKT_CORE if g in expression.columns and g != "PTEN"]
    subset = expression[pi3k_genes]
    zscored = subset.apply(stats.zscore, nan_policy="omit")
    pathway_score = zscored.mean(axis=1)
    pathway_score.name = "pi3k_pathway_zscore"

    # Load PRISM
    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    # Key PI3K/AKT drugs from H13c top hits
    key_drugs = ["ALPELISIB", "GSK2110183", "MK-2206", "GDC-0077",
                 "AZD8835", "TASELISIB", "GDC-0068", "AZD5363"]

    # Resolve column names from PRISM compound info
    name_col = _find_compound_col(
        compound_info, ["Drug.Name", "name", "Name", "drug_name", "compound_name"],
    )
    id_col = _find_compound_col(
        compound_info, ["IDs", "broad_id", "column_name", "compound_id"],
    )

    _print(f"  Compound info columns: name={name_col}, id={id_col}")
    _print("  Testing: does higher PI3K pathway expression → greater drug sensitivity?")
    _print("")
    _print(f"  {'Drug':<20} {'Spearman rho':>14} {'p-value':>12} {'n':>6} {'Interpretation'}")
    _print("  " + "-" * 80)

    for drug_name in key_drugs:
        # Find compound IDs for this drug
        mask = compound_info[name_col].str.upper() == drug_name.upper()
        if mask.sum() == 0:
            _print(f"  {drug_name:<20} NOT FOUND IN PRISM")
            continue

        brd_ids = compound_info.loc[mask, id_col].tolist()

        for brd_id in brd_ids:
            if brd_id not in sensitivity.columns:
                continue

            drug_sens = sensitivity[brd_id].dropna()
            shared = sorted(set(pathway_score.dropna().index) & set(drug_sens.index))
            if len(shared) < 30:
                continue

            rho, p_val = stats.spearmanr(pathway_score[shared], drug_sens[shared])
            # Negative rho means higher pathway → more sensitivity (more negative AUC)
            interp = "HIGHER PI3K → MORE SENSITIVE" if rho < 0 else "no predictive link"
            sig = _sig_stars(p_val)
            _print(
                f"  {drug_name:<20} {rho:>14.4f} {p_val:>12.4e} {len(shared):>6} {interp} {sig}"
            )
            break  # one compound ID per drug is enough

    # Also test: vulnerability score vs same drugs
    _print("")
    _print("  VULNERABILITY SCORE vs SAME DRUGS:")
    vuln = scores_df["vulnerability_score"]

    for drug_name in key_drugs:
        mask = compound_info[name_col].str.upper() == drug_name.upper()
        if mask.sum() == 0:
            continue

        brd_ids = compound_info.loc[mask, id_col].tolist()
        for brd_id in brd_ids:
            if brd_id not in sensitivity.columns:
                continue
            drug_sens = sensitivity[brd_id].dropna()
            shared = sorted(set(vuln.dropna().index) & set(drug_sens.index))
            if len(shared) < 30:
                continue
            rho, p_val = stats.spearmanr(vuln[shared], drug_sens[shared])
            sig = _sig_stars(p_val)
            _print(
                f"  {drug_name:<20} rho={rho:+.4f}, p={p_val:.4e} {sig}"
            )
            break


# =========================================================================
# ANALYSIS 5: PTEN loss check
# =========================================================================


def analysis_5_pten_loss(
    expression: pd.DataFrame,
    crispr: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Check if target cells have lower PTEN (removing the brake on PI3K)."""
    header = "ANALYSIS 5: PTEN STATUS IN TARGET vs BACKGROUND"
    _print_header(header)

    if "PTEN" not in expression.columns:
        _print("  PTEN not in expression matrix")
        return

    t_pten = expression.loc[expression.index.isin(target_ids), "PTEN"].dropna()
    b_pten = expression.loc[expression.index.isin(background_ids), "PTEN"].dropna()

    _print(f"  Target PTEN expression: mean={t_pten.mean():.4f}, median={t_pten.median():.4f}")
    _print(f"  Background PTEN expression: mean={b_pten.mean():.4f}, median={b_pten.median():.4f}")

    t_stat, p_val = stats.ttest_ind(t_pten, b_pten)
    d = (t_pten.mean() - b_pten.mean()) / np.sqrt(
        (t_pten.std() ** 2 + b_pten.std() ** 2) / 2
    )
    _print(f"  Cohen's d: {d:+.4f} ({'target LOWER PTEN' if d < 0 else 'target higher PTEN'})")
    _print(f"  p-value: {p_val:.4e} {_sig_stars(p_val)}")

    # Fraction with very low PTEN (< 25th percentile)
    pten_25 = expression["PTEN"].quantile(0.25)
    t_low = (t_pten < pten_25).mean()
    b_low = (b_pten < pten_25).mean()
    _print(f"  Fraction PTEN-low (<25th pctl): target={t_low:.3f}, bg={b_low:.3f}")

    # PTEN dependency (should be non-essential if PTEN is a tumor suppressor)
    if "PTEN" in crispr.columns:
        _print("")
        t_dep = crispr.loc[crispr.index.isin(target_ids), "PTEN"].dropna()
        b_dep = crispr.loc[crispr.index.isin(background_ids), "PTEN"].dropna()
        _print(f"  PTEN CRISPR dependency:")
        _print(f"    Target: mean={t_dep.mean():.4f} (positive = non-essential)")
        _print(f"    Background: mean={b_dep.mean():.4f}")


# =========================================================================
# Utilities
# =========================================================================


def _print(msg: str) -> None:
    """Print to stdout for script output."""
    logger.info(msg.strip() if msg.strip() else "---")


def _print_header(title: str) -> None:
    """Print a section header."""
    border = "=" * 90
    logger.info(border)
    logger.info(title)
    logger.info(border)


def _sig_stars(p: float) -> str:
    """Return significance stars."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.1:
        return "."
    return ""


# =========================================================================
# Main
# =========================================================================


def main() -> None:
    """Run all PI3K/AKT pathway validation analyses."""
    border = "=" * 90
    logger.info(border)
    logger.info("H13d: PI3K/AKT PATHWAY ACTIVITY VALIDATION")
    logger.info(border)

    logger.info("Loading DepMap data")
    expression, crispr, scores_df, target_ids, background_ids, metadata = _load_data()
    logger.info("Data loaded")

    # Run all analyses
    expr_results = analysis_1_pathway_expression(expression, target_ids, background_ids)
    corr_results = analysis_2_correlation(expression, scores_df)
    analysis_3_pdac_specific(expression, metadata, target_ids, background_ids)
    analysis_4_pathway_vs_drugs(expression, scores_df)
    analysis_5_pten_loss(expression, crispr, target_ids, background_ids)

    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    expr_results.to_csv(RESULTS_DIR / "pi3k_expression_target_vs_bg.csv", index=False)
    corr_results.to_csv(RESULTS_DIR / "pi3k_vuln_correlation.csv", index=False)

    logger.info(border)
    logger.info("H13d PI3K/AKT VALIDATION COMPLETE")
    logger.info(border)


if __name__ == "__main__":
    main()
