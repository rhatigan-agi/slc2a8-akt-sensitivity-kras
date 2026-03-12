"""H13c Deep Dive: Unbiased PRISM screen, MOA enrichment, per-gene driver analysis.

Extends H13b results with:
  1. SAR405 / PIK-93 (VPS34 inhibitor) status check
  2. Unbiased genome-wide PRISM screen (ALL compounds)
  3. MOA enrichment analysis (Fisher exact test)
  4. Autophagy inducer negative control
  5. Per-gene driver decomposition (which initiation gene drives the signal?)
  6. Pan-GI pooled analysis (reframe from PDAC-specific to GI vulnerability)
  7. PIK3C3 individual gene dependency deep dive
"""

from collections import Counter
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
from pdac.h13_trehalose_vulnerability.enrichment import check_dependency_enrichment
from pdac.h13_trehalose_vulnerability.gene_sets import (
    GeneSetCollection,
    RAS_PATHWAY_MUTATIONS,
    TARGET_GENE,
)
from pdac.h13_trehalose_vulnerability.h13b_gene_sets import (
    AUTOPHAGY_EXECUTION,
    AUTOPHAGY_INITIATION,
)
from pdac.h13_trehalose_vulnerability.h13b_lineage import get_lineage_model_ids
from pdac.h13_trehalose_vulnerability.h13b_prism import (
    _find_compound_col,
    compute_drug_sensitivity,
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
RESULTS_DIR = Path("results/h13c")


def _load_data() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame,
    list[str], list[str], list[str], list[str],
]:
    """Load all DepMap data and compute target/background groups."""
    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    crispr = load_crispr_dependencies(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)

    gene_sets = GeneSetCollection()
    pdac_ids = get_pdac_model_ids(metadata)
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

    return (
        expression, crispr, metadata, mutations,
        target_ids, background_ids, ras_mutant_ids, common_ids,
    )


# =========================================================================
# ANALYSIS 1: SAR405 / PIK-93 status check
# =========================================================================


def check_vps34_inhibitors(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Check whether VPS34-specific inhibitors (SAR405, PIK-93) were tested."""
    print("\n" + "=" * 90)
    print("ANALYSIS 1: VPS34 INHIBITOR STATUS CHECK")
    print("=" * 90)

    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])
    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])

    vps34_drugs = ["SAR405", "PIK-93", "PIK-III", "VPS34-IN"]
    for drug in vps34_drugs:
        mask = compound_info[name_col].str.contains(drug, case=False, na=False)
        matches = compound_info[mask]

        if matches.empty:
            print(f"\n  {drug}: NOT in compound list")
            continue

        for _, row in matches.iterrows():
            cid = str(row[id_col])
            print(f"\n  {drug} (compound ID: {cid})")
            print(f"    MOA: {row.get('MOA', 'N/A')}")
            print(f"    Target: {row.get('repurposing_target', 'N/A')}")

            # Check if in sensitivity matrix
            if cid in sensitivity.columns:
                t_vals = sensitivity.loc[
                    [i for i in target_ids if i in sensitivity.index], cid
                ].dropna()
                b_vals = sensitivity.loc[
                    [i for i in background_ids if i in sensitivity.index], cid
                ].dropna()

                print(f"    In sensitivity matrix: YES")
                print(f"    Target lines with data: {len(t_vals)}")
                print(f"    Background lines with data: {len(b_vals)}")

                if len(t_vals) >= 3 and len(b_vals) >= 3:
                    u_stat, p_val = stats.mannwhitneyu(
                        t_vals.values, b_vals.values, alternative="less",
                    )
                    mean_diff = float(np.mean(t_vals.values) - np.mean(b_vals.values))
                    print(f"    p-value (one-sided): {p_val:.4e}")
                    print(f"    mean_diff: {mean_diff:.4f}")
                    print(f"    target mean sensitivity: {np.mean(t_vals.values):.4f}")
                    print(f"    background mean sensitivity: {np.mean(b_vals.values):.4f}")
                else:
                    print(f"    Insufficient data for test")
            else:
                # Try partial match
                partial = [c for c in sensitivity.columns if cid.split("-")[0] in c]
                if partial:
                    print(f"    Exact ID not found, partial matches: {partial[:3]}")
                else:
                    print(f"    In sensitivity matrix: NO")


# =========================================================================
# ANALYSIS 2: Unbiased genome-wide PRISM screen
# =========================================================================


def unbiased_prism_screen(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
    top_n: int = 50,
) -> pd.DataFrame:
    """Test ALL PRISM compounds, not just autophagy-related ones.

    Args:
        sensitivity: Full PRISM matrix.
        compound_info: Compound metadata.
        target_ids: SLC2A8-high / RAS-mut cell line IDs.
        background_ids: All other cell line IDs.
        top_n: Number of top compounds to report.

    Returns:
        DataFrame with all tested compounds ranked by p-value.
    """
    print("\n" + "=" * 90)
    print("ANALYSIS 2: UNBIASED GENOME-WIDE PRISM SCREEN")
    print("=" * 90)

    all_compounds = list(sensitivity.columns)
    print(f"\n  Total compounds in PRISM matrix: {len(all_compounds)}")

    result_df = compute_drug_sensitivity(
        sensitivity, target_ids, background_ids, all_compounds,
    )

    if result_df.empty:
        print("  No compounds tested successfully")
        return result_df

    # Resolve names
    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])
    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    moa_col = _find_compound_col(compound_info, ["MOA", "moa", "mechanism_of_action"])
    target_col = _find_compound_col(compound_info, ["repurposing_target", "target"])

    # Build lookup
    id_to_name: dict[str, str] = {}
    id_to_moa: dict[str, str] = {}
    id_to_target: dict[str, str] = {}
    for _, row in compound_info.iterrows():
        cid = str(row.get(id_col, ""))
        id_to_name[cid] = str(row.get(name_col, "UNKNOWN"))
        id_to_moa[cid] = str(row.get(moa_col, "UNKNOWN"))
        id_to_target[cid] = str(row.get(target_col, "UNKNOWN"))

    result_df["drug_name"] = result_df["compound"].map(id_to_name).fillna("UNKNOWN")
    result_df["moa"] = result_df["compound"].map(id_to_moa).fillna("UNKNOWN")
    result_df["drug_target"] = result_df["compound"].map(id_to_target).fillna("UNKNOWN")

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    _, fdr_pvals, _, _ = multipletests(result_df["p_value"].values, method="fdr_bh")
    result_df["fdr_q"] = fdr_pvals

    n_tested = len(result_df)
    n_nom_sig = (result_df["p_value"] < 0.05).sum()
    n_fdr_sig = (result_df["fdr_q"] < 0.1).sum()

    print(f"  Compounds tested: {n_tested}")
    print(f"  Nominal p < 0.05: {n_nom_sig} ({100*n_nom_sig/n_tested:.1f}%)")
    print(f"  FDR q < 0.10: {n_fdr_sig}")
    print(f"\n  Expected by chance at p<0.05: {n_tested * 0.05:.0f}")
    print(f"  Observed: {n_nom_sig}")
    print(f"  Enrichment ratio: {n_nom_sig / (n_tested * 0.05):.2f}x")

    print(f"\n  TOP {top_n} COMPOUNDS (by p-value):")
    print(f"  {'Rank':<5} {'Drug':<30} {'MOA':<40} {'p-value':<12} {'FDR q':<12} {'mean_diff':<10}")
    print("  " + "-" * 105)
    for i, (_, row) in enumerate(result_df.head(top_n).iterrows()):
        drug = str(row["drug_name"])[:28]
        moa = str(row["moa"])[:38]
        sig = "***" if row["p_value"] < 0.01 else "**" if row["p_value"] < 0.05 else "*" if row["p_value"] < 0.1 else ""
        print(
            f"  {i+1:<5} {drug:<30} {moa:<40} "
            f"{row['p_value']:<12.4e} {row['fdr_q']:<12.4e} {row['mean_diff']:<10.4f} {sig}"
        )

    return result_df


# =========================================================================
# ANALYSIS 3: MOA enrichment (Fisher exact test)
# =========================================================================


def moa_enrichment_analysis(
    result_df: pd.DataFrame,
    p_threshold: float = 0.05,
    min_moa_count: int = 3,
) -> None:
    """Test whether specific MOA classes are enriched in significant compounds.

    For each MOA class, test:
      H0: proportion of significant compounds is same as overall
      H1: MOA class is enriched among significant hits

    Uses Fisher's exact test.
    """
    print("\n" + "=" * 90)
    print("ANALYSIS 3: MOA ENRICHMENT (Fisher exact test)")
    print("=" * 90)

    sig_mask = result_df["p_value"] < p_threshold
    n_sig = sig_mask.sum()
    n_total = len(result_df)
    baseline_rate = n_sig / n_total

    print(f"\n  Overall: {n_sig}/{n_total} compounds significant (p < {p_threshold})")
    print(f"  Baseline rate: {baseline_rate:.3f}")

    # Parse MOA — compounds can have multiple MOA terms
    moa_counts: Counter[str] = Counter()
    moa_sig_counts: Counter[str] = Counter()

    for _, row in result_df.iterrows():
        moa_str = str(row["moa"]).upper()
        if moa_str in ("NAN", "UNKNOWN", ""):
            continue

        # Split on comma to get individual MOA terms
        terms = [t.strip() for t in moa_str.split(",")]
        is_sig = row["p_value"] < p_threshold

        for term in terms:
            if len(term) > 2:
                moa_counts[term] += 1
                if is_sig:
                    moa_sig_counts[term] += 1

    # Filter to MOA classes with enough compounds
    eligible_moas = {m: c for m, c in moa_counts.items() if c >= min_moa_count}
    print(f"\n  MOA classes with >= {min_moa_count} compounds: {len(eligible_moas)}")

    # Fisher exact test per MOA
    enrichment_results: list[dict[str, object]] = []
    for moa_term, count in sorted(eligible_moas.items(), key=lambda x: -x[1]):
        sig_in_moa = moa_sig_counts.get(moa_term, 0)
        not_sig_in_moa = count - sig_in_moa
        sig_not_in_moa = n_sig - sig_in_moa
        not_sig_not_in_moa = (n_total - n_sig) - not_sig_in_moa

        # 2x2 contingency table
        table = np.array([
            [sig_in_moa, not_sig_in_moa],
            [sig_not_in_moa, not_sig_not_in_moa],
        ])

        odds_ratio, fisher_p = stats.fisher_exact(table, alternative="greater")

        enrichment_results.append({
            "moa": moa_term,
            "n_total": count,
            "n_sig": sig_in_moa,
            "rate": sig_in_moa / count if count > 0 else 0,
            "odds_ratio": odds_ratio,
            "fisher_p": fisher_p,
        })

    enrich_df = pd.DataFrame(enrichment_results).sort_values("fisher_p")

    print(f"\n  TOP MOA CLASSES (by Fisher enrichment p-value):")
    print(f"  {'MOA':<45} {'Sig/Total':<12} {'Rate':<8} {'OR':<8} {'p-value':<12}")
    print("  " + "-" * 85)
    for _, row in enrich_df.head(25).iterrows():
        sig_marker = "***" if row["fisher_p"] < 0.01 else "**" if row["fisher_p"] < 0.05 else "*" if row["fisher_p"] < 0.1 else ""
        print(
            f"  {str(row['moa']):<45} "
            f"{row['n_sig']}/{row['n_total']:<9} "
            f"{row['rate']:<8.3f} {row['odds_ratio']:<8.2f} "
            f"{row['fisher_p']:<12.4e} {sig_marker}"
        )

    # Highlight PI3K/mTOR specifically
    print("\n  PI3K/mTOR SPECIFIC:")
    pi3k_mtor_terms = [t for t in eligible_moas if "PI3K" in t or "MTOR" in t or "AKT" in t]
    for term in pi3k_mtor_terms:
        row = enrich_df[enrich_df["moa"] == term]
        if not row.empty:
            r = row.iloc[0]
            print(
                f"    {term}: {r['n_sig']}/{r['n_total']} sig, "
                f"OR={r['odds_ratio']:.2f}, p={r['fisher_p']:.4e}"
            )


# =========================================================================
# ANALYSIS 4: Autophagy inducer negative control
# =========================================================================


def autophagy_inducer_control(
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Check if autophagy INDUCERS show opposite effect (less sensitive)."""
    print("\n" + "=" * 90)
    print("ANALYSIS 4: AUTOPHAGY INDUCER NEGATIVE CONTROL")
    print("=" * 90)

    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])
    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    moa_col = _find_compound_col(compound_info, ["MOA", "moa"])

    # Find autophagy inducers
    inducer_mask = compound_info[moa_col].str.contains("AUTOPHAGY INDUCER", case=False, na=False)
    inducers = compound_info[inducer_mask]

    # Also include rapamycin/everolimus (mTOR inhibitors are autophagy inducers)
    rapa_mask = compound_info[name_col].str.contains(
        "rapamycin|everolimus|temsirolimus", case=False, na=False,
    )
    rapa_compounds = compound_info[rapa_mask]

    print(f"\n  Autophagy inducers found: {len(inducers)}")
    print(f"  Rapalogs found: {len(rapa_compounds)}")

    all_inducer_compounds = pd.concat([inducers, rapa_compounds]).drop_duplicates(subset=[id_col])

    for _, row in all_inducer_compounds.iterrows():
        cid = str(row[id_col])
        drug = str(row[name_col])
        moa = str(row[moa_col])

        if cid not in sensitivity.columns:
            print(f"\n  {drug} ({cid}): NOT in sensitivity matrix")
            continue

        t_vals = sensitivity.loc[
            [i for i in target_ids if i in sensitivity.index], cid
        ].dropna()
        b_vals = sensitivity.loc[
            [i for i in background_ids if i in sensitivity.index], cid
        ].dropna()

        if len(t_vals) < 3 or len(b_vals) < 3:
            print(f"\n  {drug}: insufficient data (target={len(t_vals)}, bg={len(b_vals)})")
            continue

        # Two-sided test (inducers could go either way)
        _, p_two = stats.mannwhitneyu(t_vals.values, b_vals.values, alternative="two-sided")
        # One-sided: target LESS sensitive (more resistant = positive diff)
        _, p_greater = stats.mannwhitneyu(t_vals.values, b_vals.values, alternative="greater")
        mean_diff = float(np.mean(t_vals.values) - np.mean(b_vals.values))

        direction = "MORE SENSITIVE" if mean_diff < 0 else "LESS SENSITIVE (protective)"
        print(f"\n  {drug} ({moa})")
        print(f"    mean_diff: {mean_diff:.4f} ({direction})")
        print(f"    p (two-sided): {p_two:.4e}")
        print(f"    p (target less sensitive): {p_greater:.4e}")
        print(f"    target n={len(t_vals)}, background n={len(b_vals)}")


# =========================================================================
# ANALYSIS 5: Per-gene driver decomposition
# =========================================================================


def per_gene_driver_analysis(
    crispr: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Decompose which initiation genes drive the dependency signal."""
    print("\n" + "=" * 90)
    print("ANALYSIS 5: PER-GENE DRIVER DECOMPOSITION")
    print("=" * 90)

    all_genes = AUTOPHAGY_INITIATION + AUTOPHAGY_EXECUTION

    print(f"\n  {'Gene':<12} {'Set':<15} {'Target mean':<14} {'BG mean':<14} {'Diff':<10} {'p-value':<14} {'Effect size':<12}")
    print("  " + "-" * 95)

    for gene in all_genes:
        if gene not in crispr.columns:
            print(f"  {gene:<12} NOT IN DATA")
            continue

        set_label = "INITIATION" if gene in AUTOPHAGY_INITIATION else "EXECUTION"

        t_vals = crispr.loc[[i for i in target_ids if i in crispr.index], gene].dropna().values
        b_vals = crispr.loc[[i for i in background_ids if i in crispr.index], gene].dropna().values

        if len(t_vals) < 3 or len(b_vals) < 3:
            continue

        _, p_val = stats.mannwhitneyu(t_vals, b_vals, alternative="less")
        mean_diff = float(np.mean(t_vals) - np.mean(b_vals))

        # Cohen's d effect size
        pooled_std = np.sqrt(
            ((len(t_vals) - 1) * np.var(t_vals, ddof=1) +
             (len(b_vals) - 1) * np.var(b_vals, ddof=1)) /
            (len(t_vals) + len(b_vals) - 2)
        )
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0.0

        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(
            f"  {gene:<12} {set_label:<15} "
            f"{np.mean(t_vals):<14.4f} {np.mean(b_vals):<14.4f} "
            f"{mean_diff:<10.4f} {p_val:<14.4e} d={cohens_d:<10.3f} {sig}"
        )

    # Also test: is PIK3C3 the single strongest driver?
    print("\n  PIK3C3 (VPS34) DETAILED:")
    if "PIK3C3" in crispr.columns:
        t_vals = crispr.loc[[i for i in target_ids if i in crispr.index], "PIK3C3"].dropna()
        b_vals = crispr.loc[[i for i in background_ids if i in crispr.index], "PIK3C3"].dropna()

        print(f"    Target distribution: mean={t_vals.mean():.4f}, std={t_vals.std():.4f}, n={len(t_vals)}")
        print(f"    Background distribution: mean={b_vals.mean():.4f}, std={b_vals.std():.4f}, n={len(b_vals)}")

        # Percentile rank of target mean within background distribution
        pct_rank = (b_vals < t_vals.mean()).mean() * 100
        print(f"    Target mean falls at {pct_rank:.1f}th percentile of background")

        # How many target lines have PIK3C3 dependency < -0.5 (strong dependency)?
        strong_dep_target = (t_vals < -0.5).sum()
        strong_dep_bg = (b_vals < -0.5).sum()
        print(f"    Strong dependency (< -0.5): target {strong_dep_target}/{len(t_vals)} "
              f"({100*strong_dep_target/len(t_vals):.1f}%), "
              f"bg {strong_dep_bg}/{len(b_vals)} ({100*strong_dep_bg/len(b_vals):.1f}%)")


# =========================================================================
# ANALYSIS 6: Pan-GI pooled analysis
# =========================================================================


def pan_gi_pooled_analysis(
    crispr: pd.DataFrame,
    expression: pd.DataFrame,
    mutations: pd.DataFrame,
    metadata: pd.DataFrame,
    ras_mutant_ids: list[str],
) -> None:
    """Pool all GI lineages and test initiation dependency in the combined cohort."""
    print("\n" + "=" * 90)
    print("ANALYSIS 6: PAN-GI POOLED ANALYSIS")
    print("=" * 90)

    gi_lineages = {
        "Pancreas": "Pancrea",
        "Bowel": "Bowel",
        "Esophagus/Stomach": "Esophagus/Stomach",
        "Liver": "Liver",
        "Biliary": "Biliary",
    }

    all_gi_ids: list[str] = []
    for name, pattern in gi_lineages.items():
        ids = get_lineage_model_ids(metadata, pattern)
        available = [i for i in ids if i in expression.index and i in crispr.index]
        all_gi_ids.extend(available)
        print(f"  {name}: {len(available)} cell lines")

    all_gi_ids = sorted(set(all_gi_ids))
    print(f"\n  Total GI cell lines: {len(all_gi_ids)}")

    # Stratify by SLC2A8 within GI cohort
    if TARGET_GENE not in expression.columns:
        print("  SLC2A8 not in expression columns")
        return

    gi_expr = expression.loc[all_gi_ids, TARGET_GENE]
    median_val = gi_expr.median()
    high_ids = gi_expr[gi_expr >= median_val].index.tolist()

    # Target: SLC2A8-high AND RAS-mutant within GI
    gi_target = [i for i in high_ids if i in ras_mutant_ids]
    gi_background = [i for i in all_gi_ids if i not in gi_target]

    print(f"  GI SLC2A8-high / RAS-mut (target): {len(gi_target)}")
    print(f"  GI background: {len(gi_background)}")

    # Test initiation enrichment in pan-GI
    if len(gi_target) >= 3:
        result = check_dependency_enrichment(
            crispr, gi_target, gi_background,
            AUTOPHAGY_INITIATION, "autophagy_initiation_panGI",
        )
        print(f"\n  INITIATION enrichment (pan-GI): p = {result.set_level_p:.4e}, "
              f"mean_diff = {result.set_level_statistic:.4f}")

        exec_result = check_dependency_enrichment(
            crispr, gi_target, gi_background,
            AUTOPHAGY_EXECUTION, "autophagy_execution_panGI",
        )
        print(f"  EXECUTION control (pan-GI):     p = {exec_result.set_level_p:.4e}, "
              f"mean_diff = {exec_result.set_level_statistic:.4f}")

        # Per-gene within GI
        print(f"\n  Per-gene initiation dependency (pan-GI):")
        for gene in AUTOPHAGY_INITIATION:
            if gene not in crispr.columns:
                continue
            t_vals = crispr.loc[[i for i in gi_target if i in crispr.index], gene].dropna().values
            b_vals = crispr.loc[[i for i in gi_background if i in crispr.index], gene].dropna().values
            if len(t_vals) >= 3 and len(b_vals) >= 3:
                _, p_val = stats.mannwhitneyu(t_vals, b_vals, alternative="less")
                mean_diff = float(np.mean(t_vals) - np.mean(b_vals))
                sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
                print(f"    {gene:<10}: p = {p_val:.4e}, diff = {mean_diff:.4f} {sig}")

    # Also test: PRISM drug sensitivity in pan-GI only
    try:
        sensitivity = load_prism_sensitivity(DEPMAP_DIR)
        gi_target_prism = [i for i in gi_target if i in sensitivity.index]
        gi_bg_prism = [i for i in gi_background if i in sensitivity.index]

        if len(gi_target_prism) >= 3 and len(gi_bg_prism) >= 3:
            # Test all compounds
            all_compounds = list(sensitivity.columns)
            gi_drug_results = compute_drug_sensitivity(
                sensitivity, gi_target_prism, gi_bg_prism, all_compounds,
            )

            if not gi_drug_results.empty:
                compound_info = load_prism_compound_info(DEPMAP_DIR)
                id_col = _find_compound_col(compound_info, ["IDs", "broad_id"])
                name_col = _find_compound_col(compound_info, ["Drug.Name", "name"])
                moa_col = _find_compound_col(compound_info, ["MOA", "moa"])

                id_to_name = {}
                id_to_moa = {}
                for _, row in compound_info.iterrows():
                    cid = str(row.get(id_col, ""))
                    id_to_name[cid] = str(row.get(name_col, ""))
                    id_to_moa[cid] = str(row.get(moa_col, ""))

                gi_drug_results["drug_name"] = gi_drug_results["compound"].map(id_to_name).fillna("UNKNOWN")
                gi_drug_results["moa"] = gi_drug_results["compound"].map(id_to_moa).fillna("UNKNOWN")

                n_sig = (gi_drug_results["p_value"] < 0.05).sum()
                print(f"\n  Pan-GI PRISM: {n_sig}/{len(gi_drug_results)} compounds significant (p < 0.05)")
                print(f"  Expected by chance: {len(gi_drug_results) * 0.05:.0f}")

                print(f"\n  Top 20 pan-GI PRISM hits:")
                for i, (_, row) in enumerate(gi_drug_results.head(20).iterrows()):
                    drug = str(row["drug_name"])[:25]
                    moa = str(row["moa"])[:35]
                    sig = "***" if row["p_value"] < 0.01 else "**" if row["p_value"] < 0.05 else "*" if row["p_value"] < 0.1 else ""
                    print(f"    {i+1:<4} {drug:<27} {moa:<37} p={row['p_value']:.4e} diff={row['mean_diff']:.3f} {sig}")
    except Exception as exc:
        logger.warning("Pan-GI PRISM analysis failed", extra={"error": str(exc)})


# =========================================================================
# ANALYSIS 7: PIK3C3 deep dive — is it the single driver?
# =========================================================================


def pik3c3_deep_dive(
    crispr: pd.DataFrame,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Deep dive into PIK3C3/VPS34: expression-dependency correlation, lineage pattern."""
    print("\n" + "=" * 90)
    print("ANALYSIS 7: PIK3C3 (VPS34) DEEP DIVE")
    print("=" * 90)

    if "PIK3C3" not in crispr.columns or "PIK3C3" not in expression.columns:
        print("  PIK3C3 not in data")
        return

    # Expression-dependency correlation across ALL cell lines
    common = sorted(set(expression.index) & set(crispr.index))
    expr_vals = expression.loc[common, "PIK3C3"].dropna()
    dep_vals = crispr.loc[common, "PIK3C3"].dropna()
    common2 = sorted(set(expr_vals.index) & set(dep_vals.index))

    rho, p_val = stats.spearmanr(expr_vals.loc[common2].values, dep_vals.loc[common2].values)
    print(f"\n  PIK3C3 expression vs dependency (all cell lines):")
    print(f"    Spearman rho = {rho:.4f}, p = {p_val:.4e}")
    print(f"    n = {len(common2)}")

    # Same for SLC2A8 expression vs PIK3C3 dependency
    if TARGET_GENE in expression.columns:
        slc2a8_vals = expression.loc[common, TARGET_GENE].dropna()
        common3 = sorted(set(slc2a8_vals.index) & set(dep_vals.index))
        rho2, p_val2 = stats.spearmanr(slc2a8_vals.loc[common3].values, dep_vals.loc[common3].values)
        print(f"\n  SLC2A8 expression vs PIK3C3 dependency:")
        print(f"    Spearman rho = {rho2:.4f}, p = {p_val2:.4e}")
        print(f"    n = {len(common3)}")

    # PIK3C3 dependency by lineage (top 10 most dependent)
    dep_by_lineage: dict[str, float] = {}
    if "OncotreeLineage" in metadata.columns:
        for lineage in metadata["OncotreeLineage"].dropna().unique():
            lineage_ids = metadata[
                metadata["OncotreeLineage"] == lineage
            ].index.tolist()
            lineage_ids = [i for i in lineage_ids if i in crispr.index]
            if len(lineage_ids) >= 5:
                mean_dep = float(crispr.loc[lineage_ids, "PIK3C3"].mean())
                dep_by_lineage[lineage] = mean_dep

        sorted_lineages = sorted(dep_by_lineage.items(), key=lambda x: x[1])
        print(f"\n  PIK3C3 mean dependency by lineage (most dependent first):")
        for lineage, dep in sorted_lineages[:15]:
            marker = " <-- GI" if any(gi in lineage.lower() for gi in ["pancrea", "bowel", "esophag", "stomach", "liver", "biliar"]) else ""
            print(f"    {lineage:<30}: {dep:.4f}{marker}")

    # Remove PIK3C3 from initiation set and re-test
    print(f"\n  INITIATION WITHOUT PIK3C3:")
    reduced_init = [g for g in AUTOPHAGY_INITIATION if g != "PIK3C3"]
    result_without = check_dependency_enrichment(
        crispr, target_ids, background_ids,
        reduced_init, "initiation_no_PIK3C3",
    )
    print(f"    Genes: {reduced_init}")
    print(f"    p = {result_without.set_level_p:.4e}, mean_diff = {result_without.set_level_statistic:.4f}")

    # PIK3C3 alone
    print(f"\n  PIK3C3 ALONE:")
    result_alone = check_dependency_enrichment(
        crispr, target_ids, background_ids,
        ["PIK3C3"], "PIK3C3_alone",
    )
    print(f"    p = {result_alone.set_level_p:.4e}, mean_diff = {result_alone.set_level_statistic:.4f}")


# =========================================================================
# MAIN
# =========================================================================


def main() -> None:
    """Run all H13c deep-dive analyses."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("H13c DEEP DIVE: UNBIASED PRISM + MOA ENRICHMENT + DRIVER DECOMPOSITION")
    print("=" * 90)

    # Load data
    logger.info("Loading DepMap data")
    (
        expression, crispr, metadata, mutations,
        target_ids, background_ids, ras_mutant_ids, common_ids,
    ) = _load_data()

    logger.info(
        "Data loaded",
        extra={
            "n_common": len(common_ids),
            "n_target": len(target_ids),
            "n_background": len(background_ids),
        },
    )

    # Load PRISM data
    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    # Run all analyses
    check_vps34_inhibitors(sensitivity, compound_info, target_ids, background_ids)
    unbiased_results = unbiased_prism_screen(sensitivity, compound_info, target_ids, background_ids)
    moa_enrichment_analysis(unbiased_results)
    autophagy_inducer_control(sensitivity, compound_info, target_ids, background_ids)
    per_gene_driver_analysis(crispr, target_ids, background_ids)
    pan_gi_pooled_analysis(crispr, expression, mutations, metadata, ras_mutant_ids)
    pik3c3_deep_dive(crispr, expression, metadata, target_ids, background_ids)

    # Save unbiased results
    if not unbiased_results.empty:
        out_path = RESULTS_DIR / "unbiased_prism_results.csv"
        unbiased_results.to_csv(out_path, index=False)
        logger.info("Unbiased results saved", extra={"path": str(out_path)})

    print("\n" + "=" * 90)
    print("H13c DEEP DIVE COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
