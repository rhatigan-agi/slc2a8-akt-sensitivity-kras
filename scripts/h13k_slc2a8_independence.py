"""H13k: SLC2A8 independence validation — is the AKT signal real without PIK3CA?

Three definitive tests:
  1. Continuous correlation: SLC2A8 expression vs ALL 6,790 PRISM drugs
     (no mutation stratifier). Where do AKT inhibitors rank?
  2. PIK3CA-WT only: SLC2A8 vs AKT drugs restricted to PIK3CA-wildtype lines.
  3. KRAS-mut / PIK3CA-WT: SLC2A8 vs AKT drugs in the clinically relevant
     subgroup (PDAC-like: KRAS-mutant but PIK3CA-wildtype).

Also includes:
  4. PIK3CA-mut only: does SLC2A8 add information beyond PIK3CA status?
  5. Lineage-corrected versions of tests 2-4.

Output: results/h13k/h13k_slc2a8_independence.txt
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
RESULTS_DIR = Path("results/h13k")

AKT_DRUGS: dict[str, str] = {
    "Capivasertib": "AZD5363",
    "Afuresertib": "GSK2110183",
    "Ipatasertib": "GDC-0068",
}

PI3K_DRUGS: dict[str, str] = {
    "Alpelisib": "ALPELISIB",
    "GDC-0077": "GDC-0077",
    "Taselisib": "TASELISIB",
}

MEK_DRUGS: dict[str, str] = {
    "Trametinib": "TRAMETINIB",
    "Selumetinib": "SELUMETINIB",
    "Cobimetinib": "COBIMETINIB",
}

TOP_LINEAGES = 15


def _get_mutant_ids(mutations: pd.DataFrame, gene: str) -> set[str]:
    """Return set of ModelIDs with any somatic mutation in the given gene."""
    hugo_col = _find_column(mutations, ["HugoSymbol", "Hugo_Symbol"])
    model_col = _find_column(mutations, ["ModelID", "DepMap_ID"])
    mask = mutations[hugo_col] == gene
    return set(mutations.loc[mask, model_col].unique())


def _resolve_drug(
    name: str,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> pd.Series | None:
    """Resolve drug name to sensitivity vector."""
    mask = compound_info[name_col].str.upper() == name.upper()
    if mask.sum() == 0:
        return None
    for brd_id in compound_info.loc[mask, id_col].tolist():
        if brd_id in sensitivity.columns:
            return sensitivity[brd_id].dropna()
    return None


def _residualize(x: np.ndarray, confounders: np.ndarray) -> np.ndarray:
    """Return residuals of x after projecting out confounders (with intercept)."""
    c = np.column_stack([np.ones(len(confounders)), confounders])
    coef, _, _, _ = np.linalg.lstsq(c, x, rcond=None)
    return x - c @ coef


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
    return "ns"


# =====================================================================
# Analysis 1: Continuous SLC2A8 x ALL drugs (no mutation info)
# =====================================================================

def analysis_1_continuous_ranking(
    expression: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> list[str]:
    """Rank all PRISM drugs by Spearman correlation with SLC2A8 expression."""
    lines: list[str] = []
    lines.append("=" * 70)
    lines.append("ANALYSIS 1: Continuous SLC2A8 x drug correlation (NO mutation info)")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Method: Spearman rank correlation of SLC2A8 expression vs each")
    lines.append("PRISM drug across ALL cell lines. No stratifier. No mutation data.")
    lines.append("")

    slc2a8 = expression["SLC2A8"].dropna()
    shared_lines = sorted(set(slc2a8.index) & set(sensitivity.index))
    slc2a8_vals = slc2a8[shared_lines]

    lines.append(f"Cell lines with SLC2A8 + PRISM data: {len(shared_lines)}")

    # Build drug name lookup
    id_to_name: dict[str, str] = {}
    id_to_moa: dict[str, str] = {}
    moa_col = _find_compound_col(compound_info, ["MOA", "moa", "mechanism_of_action"])
    for _, row in compound_info.iterrows():
        cid = str(row.get(id_col, ""))
        id_to_name[cid] = str(row.get(name_col, "UNKNOWN"))
        id_to_moa[cid] = str(row.get(moa_col, "UNKNOWN"))

    results: list[dict[str, object]] = []
    for col in sensitivity.columns:
        drug_vals = sensitivity.loc[shared_lines, col].dropna()
        overlap = sorted(set(slc2a8_vals.index) & set(drug_vals.index))
        if len(overlap) < 50:
            continue
        rho, p = stats.spearmanr(slc2a8_vals[overlap].values, drug_vals[overlap].values)
        results.append({
            "compound": col,
            "drug_name": id_to_name.get(col, "UNKNOWN"),
            "moa": id_to_moa.get(col, "UNKNOWN"),
            "rho": rho,
            "p": p,
            "n": len(overlap),
        })

    result_df = pd.DataFrame(results).sort_values("rho")
    n_tested = len(result_df)

    lines.append(f"Drugs tested: {n_tested}")
    lines.append("")

    # Where do AKT inhibitors rank?
    lines.append("--- AKT INHIBITOR RANKING (most negative rho = most sensitive) ---")
    for label, drug_name in AKT_DRUGS.items():
        match = result_df[result_df["drug_name"].str.upper() == drug_name.upper()]
        if match.empty:
            # Try partial match
            match = result_df[result_df["drug_name"].str.upper().str.contains(
                drug_name.upper()[:6]
            )]
        if not match.empty:
            row = match.iloc[0]
            rank = int((result_df["rho"] <= float(row["rho"])).sum())
            pctile = 100 * rank / n_tested
            lines.append(
                f"  {label}: rho={float(row['rho']):+.4f}, p={float(row['p']):.2e} "
                f"{_sig(float(row['p']))}, n={int(row['n'])}, "
                f"rank={rank}/{n_tested} (top {pctile:.1f}%)"
            )
        else:
            lines.append(f"  {label}: NOT FOUND")

    lines.append("")
    lines.append("--- PI3K INHIBITOR RANKING (control: expect lineage-confounded) ---")
    for label, drug_name in PI3K_DRUGS.items():
        match = result_df[result_df["drug_name"].str.upper() == drug_name.upper()]
        if not match.empty:
            row = match.iloc[0]
            rank = int((result_df["rho"] <= float(row["rho"])).sum())
            pctile = 100 * rank / n_tested
            lines.append(
                f"  {label}: rho={float(row['rho']):+.4f}, p={float(row['p']):.2e} "
                f"{_sig(float(row['p']))}, n={int(row['n'])}, "
                f"rank={rank}/{n_tested} (top {pctile:.1f}%)"
            )
        else:
            lines.append(f"  {label}: NOT FOUND")

    lines.append("")
    lines.append("--- MEK INHIBITOR RANKING (control: expect RAS-dependent) ---")
    for label, drug_name in MEK_DRUGS.items():
        match = result_df[result_df["drug_name"].str.upper() == drug_name.upper()]
        if not match.empty:
            row = match.iloc[0]
            rank = int((result_df["rho"] <= float(row["rho"])).sum())
            pctile = 100 * rank / n_tested
            lines.append(
                f"  {label}: rho={float(row['rho']):+.4f}, p={float(row['p']):.2e} "
                f"{_sig(float(row['p']))}, n={int(row['n'])}, "
                f"rank={rank}/{n_tested} (top {pctile:.1f}%)"
            )
        else:
            lines.append(f"  {label}: NOT FOUND")

    # Top 20 most negatively correlated drugs
    lines.append("")
    lines.append("--- TOP 30 MOST NEGATIVELY CORRELATED DRUGS WITH SLC2A8 ---")
    lines.append(f"{'Rank':<6} {'Drug':<30} {'MOA':<40} {'rho':<10} {'p':<12} {'n':<6}")
    lines.append("-" * 104)
    for i, (_, row) in enumerate(result_df.head(30).iterrows()):
        lines.append(
            f"{i+1:<6} {str(row['drug_name'])[:28]:<30} "
            f"{str(row['moa'])[:38]:<40} "
            f"{float(row['rho']):+.4f}    {float(row['p']):<12.2e} {int(row['n']):<6}"
        )

    # MOA enrichment among top 5% most negatively correlated
    lines.append("")
    top_pct = 0.05
    n_top = int(n_tested * top_pct)
    top_drugs = result_df.head(n_top)
    lines.append(f"--- MOA ENRICHMENT IN TOP {top_pct*100:.0f}% NEGATIVE CORRELATIONS (n={n_top}) ---")

    from collections import Counter
    top_moas: Counter[str] = Counter()
    all_moas: Counter[str] = Counter()
    for _, row in result_df.iterrows():
        moa_str = str(row["moa"]).upper()
        if moa_str in ("NAN", "UNKNOWN", ""):
            continue
        for term in [t.strip() for t in moa_str.split(",")]:
            if len(term) > 2:
                all_moas[term] += 1
    for _, row in top_drugs.iterrows():
        moa_str = str(row["moa"]).upper()
        if moa_str in ("NAN", "UNKNOWN", ""):
            continue
        for term in [t.strip() for t in moa_str.split(",")]:
            if len(term) > 2:
                top_moas[term] += 1

    # Fisher exact for key MOAs
    for moa_label in ["AKT INHIBITOR", "PI3K INHIBITOR", "MEK INHIBITOR",
                       "MTOR INHIBITOR", "RAF INHIBITOR"]:
        in_top = top_moas.get(moa_label, 0)
        total = all_moas.get(moa_label, 0)
        if total < 3:
            lines.append(f"  {moa_label}: {in_top}/{total} in top (too few to test)")
            continue
        not_in_top = total - in_top
        other_in_top = n_top - in_top  # approximate
        other_not_in_top = n_tested - n_top - not_in_top
        table = np.array([[in_top, not_in_top], [other_in_top, other_not_in_top]])
        odds_r, fisher_p = stats.fisher_exact(table, alternative="greater")
        lines.append(
            f"  {moa_label}: {in_top}/{total} in top {top_pct*100:.0f}%, "
            f"OR={odds_r:.2f}, p={fisher_p:.2e} {_sig(fisher_p)}"
        )

    return lines


# =====================================================================
# Analysis 2-4: Subgroup analyses
# =====================================================================

def analysis_subgroup(
    expression: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    metadata: pd.DataFrame,
    name_col: str,
    id_col: str,
    subgroup_ids: list[str],
    subgroup_label: str,
    lineage_correct: bool = True,
) -> list[str]:
    """Correlate SLC2A8 with AKT/PI3K/MEK drugs within a subgroup."""
    lines: list[str] = []
    lines.append("")
    lines.append(f"--- {subgroup_label} (n={len(subgroup_ids)}) ---")

    slc2a8 = expression["SLC2A8"].dropna()
    shared_base = sorted(set(subgroup_ids) & set(slc2a8.index) & set(sensitivity.index))

    if len(shared_base) < 20:
        lines.append(f"  Insufficient overlap: n={len(shared_base)}")
        return lines

    lines.append(f"  Lines with SLC2A8 + PRISM: {len(shared_base)}")

    # Lineage info for correction
    lineage_col = "OncotreeLineage"
    lineage_map = metadata[lineage_col].reindex(expression.index).dropna()
    top_lineages = lineage_map.value_counts().head(TOP_LINEAGES).index.tolist()

    all_drugs = {**AKT_DRUGS, **PI3K_DRUGS, **MEK_DRUGS}

    for label, drug_name in all_drugs.items():
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            lines.append(f"  {label}: NOT FOUND")
            continue

        shared = sorted(set(shared_base) & set(sens.index))
        if len(shared) < 15:
            lines.append(f"  {label}: n={len(shared)} (too few)")
            continue

        slc_vals = slc2a8[shared].values
        drug_vals = sens[shared].values

        # Raw correlation
        rho, p = stats.spearmanr(slc_vals, drug_vals)
        result_str = f"  {label:<15} raw: rho={rho:+.4f}, p={p:.4f} {_sig(p)}, n={len(shared)}"

        # Lineage-corrected
        if lineage_correct:
            lin = lineage_map.reindex(shared).fillna("Other")
            shared_with_lin = [s for s in shared if s in lineage_map.index]
            if len(shared_with_lin) >= 30:
                lin_vals = lineage_map[shared_with_lin]
                lin_dummies = np.column_stack([
                    (lin_vals == l).astype(float).values for l in top_lineages
                ])
                slc_sub = slc2a8[shared_with_lin].values
                drug_sub = sens[shared_with_lin].values
                resid_drug = _residualize(drug_sub, lin_dummies)
                resid_slc = _residualize(slc_sub, lin_dummies)
                r_partial, p_partial = stats.spearmanr(resid_slc, resid_drug)
                result_str += f"  |  lineage-adj: r={r_partial:+.4f}, p={p_partial:.4f} {_sig(p_partial)}"
            else:
                result_str += "  |  lineage-adj: too few"

        lines.append(result_str)

    return lines


def main() -> None:
    """Run SLC2A8 independence validation."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("H13k: SLC2A8 independence validation")

    # Load data
    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)
    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])

    # Mutation subgroups
    pik3ca_ids = _get_mutant_ids(mutations, "PIK3CA")
    kras_ids = _get_mutant_ids(mutations, "KRAS")

    all_ids = set(expression.index) & set(sensitivity.index)
    pik3ca_wt_ids = sorted(all_ids - pik3ca_ids)
    pik3ca_mut_ids = sorted(all_ids & pik3ca_ids)
    kras_mut_pik3ca_wt_ids = sorted((all_ids & kras_ids) - pik3ca_ids)
    kras_mut_ids = sorted(all_ids & kras_ids)

    logger.info(
        "Subgroup sizes",
        extra={
            "all": len(all_ids),
            "pik3ca_wt": len(pik3ca_wt_ids),
            "pik3ca_mut": len(pik3ca_mut_ids),
            "kras_mut_pik3ca_wt": len(kras_mut_pik3ca_wt_ids),
            "kras_mut": len(kras_mut_ids),
        },
    )

    output: list[str] = []
    output.append("H13k: SLC2A8 Independence Validation")
    output.append("=" * 70)
    output.append("")
    output.append("Central question: Does SLC2A8 expression predict AKT inhibitor")
    output.append("sensitivity INDEPENDENTLY of PIK3CA mutation status?")
    output.append("")
    output.append(f"Total lines with expression + PRISM: {len(all_ids)}")
    output.append(f"PIK3CA-wildtype:                     {len(pik3ca_wt_ids)}")
    output.append(f"PIK3CA-mutant:                       {len(pik3ca_mut_ids)}")
    output.append(f"KRAS-mutant (any PIK3CA):             {len(kras_mut_ids)}")
    output.append(f"KRAS-mutant / PIK3CA-wildtype:        {len(kras_mut_pik3ca_wt_ids)}")
    output.append("")

    # Analysis 1: Continuous ranking
    output.extend(
        analysis_1_continuous_ranking(
            expression, sensitivity, compound_info, name_col, id_col,
        )
    )

    # Analysis 2-5: Subgroup analyses
    output.append("")
    output.append("")
    output.append("=" * 70)
    output.append("ANALYSIS 2-5: Subgroup SLC2A8 x drug correlations")
    output.append("=" * 70)
    output.append("")
    output.append("For each subgroup: raw Spearman correlation of SLC2A8 expression")
    output.append("vs drug sensitivity, plus lineage-corrected partial correlation.")
    output.append("")
    output.append("KEY: If SLC2A8 predicts AKT sensitivity in PIK3CA-WT lines,")
    output.append("the signal is INDEPENDENT of PIK3CA mutation.")

    for subgroup_ids, label in [
        (sorted(all_ids), "ALL LINES"),
        (pik3ca_wt_ids, "PIK3CA-WILDTYPE ONLY (critical independence test)"),
        (kras_mut_pik3ca_wt_ids, "KRAS-MUTANT / PIK3CA-WILDTYPE (clinical subgroup)"),
        (pik3ca_mut_ids, "PIK3CA-MUTANT ONLY (does SLC2A8 add info?)"),
        (kras_mut_ids, "KRAS-MUTANT (any PIK3CA status)"),
    ]:
        output.extend(
            analysis_subgroup(
                expression, sensitivity, compound_info, metadata,
                name_col, id_col, subgroup_ids, label,
            )
        )

    # Summary
    output.append("")
    output.append("")
    output.append("=" * 70)
    output.append("INTERPRETATION GUIDE")
    output.append("=" * 70)
    output.append("")
    output.append("If AKT drugs show significant negative correlation with SLC2A8 in")
    output.append("PIK3CA-WT lines: SLC2A8 is an INDEPENDENT biomarker (signal is real).")
    output.append("")
    output.append("If AKT drugs lose significance in PIK3CA-WT but retain it in")
    output.append("PIK3CA-MUT: the signal is MEDIATED by PIK3CA (not independent).")
    output.append("")
    output.append("If AKT drugs rank in top 5-10% of continuous correlation: SLC2A8")
    output.append("has a genuine association with AKT sensitivity independent of any")
    output.append("mutation stratifier.")

    output_text = "\n".join(output)
    output_path = RESULTS_DIR / "h13k_slc2a8_independence.txt"
    output_path.write_text(output_text)
    logger.info("Results written to %s", output_path)
    print(output_text)


if __name__ == "__main__":
    main()
