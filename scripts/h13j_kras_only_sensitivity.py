"""H13j: KRAS-only sensitivity analysis — exclude PIK3CA from stratifier.

Re-runs the PRISM MOA enrichment using KRAS/NRAS/HRAS/BRAF only
(excluding PIK3CA) as the RAS-pathway mutation stratifier. Compares
results to the original 5-gene stratification to confirm the AKT
inhibitor signal is not driven by PIK3CA inclusion.

Output: results/h13j/h13j_kras_only_sensitivity.txt
"""

from collections import Counter
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
from pdac.h13_trehalose_vulnerability.gene_sets import (
    GeneSetCollection,
    RAS_PATHWAY_MUTATIONS,
)
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
RESULTS_DIR = Path("results/h13j")

# The 4-gene list excluding PIK3CA
KRAS_ONLY_GENES: list[str] = ["KRAS", "NRAS", "HRAS", "BRAF"]


def _moa_enrichment(
    result_df: pd.DataFrame,
    p_threshold: float = 0.05,
    min_moa_count: int = 3,
) -> pd.DataFrame:
    """Compute MOA enrichment via Fisher exact test. Returns DataFrame."""
    sig_mask = result_df["p_value"] < p_threshold
    n_sig = int(sig_mask.sum())
    n_total = len(result_df)

    moa_counts: Counter[str] = Counter()
    moa_sig_counts: Counter[str] = Counter()

    for _, row in result_df.iterrows():
        moa_str = str(row["moa"]).upper()
        if moa_str in ("NAN", "UNKNOWN", ""):
            continue
        terms = [t.strip() for t in moa_str.split(",")]
        is_sig = row["p_value"] < p_threshold
        for term in terms:
            if len(term) > 2:
                moa_counts[term] += 1
                if is_sig:
                    moa_sig_counts[term] += 1

    eligible_moas = {m: c for m, c in moa_counts.items() if c >= min_moa_count}

    enrichment_results: list[dict[str, object]] = []
    for moa_term, count in sorted(eligible_moas.items(), key=lambda x: -x[1]):
        sig_in_moa = moa_sig_counts.get(moa_term, 0)
        not_sig_in_moa = count - sig_in_moa
        sig_not_in_moa = n_sig - sig_in_moa
        not_sig_not_in_moa = (n_total - n_sig) - not_sig_in_moa

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

    return pd.DataFrame(enrichment_results).sort_values("fisher_p")


def _run_screen(
    expression: pd.DataFrame,
    mutations: pd.DataFrame,
    crispr: pd.DataFrame,
    metadata: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    ras_genes: list[str],
    label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, int, int]:
    """Run PRISM screen with a given RAS gene list. Return (screen_df, moa_df, n_target, n_bg)."""
    gene_sets = GeneSetCollection()
    common_ids = sorted(
        set(expression.index) & set(crispr.index) & set(metadata.index)
    )
    expr_common = expression.loc[common_ids]

    ras_mutant_ids = get_ras_mutant_model_ids(mutations, ras_genes)

    scores_df = compute_vulnerability_score(
        expr_common, gene_sets.target, gene_sets.scoring_sets(),
    )
    quadrant_labels = stratify_by_target_and_ras(scores_df, ras_mutant_ids)
    target_ids = quadrant_labels[
        quadrant_labels == "SLC2A8-high / RAS-mut"
    ].index.tolist()
    background_ids = [i for i in common_ids if i not in target_ids]

    logger.info(
        "%s stratification",
        label,
        extra={"n_target": len(target_ids), "n_background": len(background_ids)},
    )

    all_compounds = list(sensitivity.columns)
    result_df = compute_drug_sensitivity(
        sensitivity, target_ids, background_ids, all_compounds,
    )

    # Annotate with MOA
    id_col = _find_compound_col(compound_info, ["IDs", "broad_id", "column_name"])
    name_col = _find_compound_col(compound_info, ["Drug.Name", "name", "Name"])
    moa_col = _find_compound_col(compound_info, ["MOA", "moa", "mechanism_of_action"])

    id_to_name: dict[str, str] = {}
    id_to_moa: dict[str, str] = {}
    for _, row in compound_info.iterrows():
        cid = str(row.get(id_col, ""))
        id_to_name[cid] = str(row.get(name_col, "UNKNOWN"))
        id_to_moa[cid] = str(row.get(moa_col, "UNKNOWN"))

    result_df["drug_name"] = result_df["compound"].map(id_to_name).fillna("UNKNOWN")
    result_df["moa"] = result_df["compound"].map(id_to_moa).fillna("UNKNOWN")

    moa_df = _moa_enrichment(result_df)

    return result_df, moa_df, len(target_ids), len(background_ids)


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


KEY_MOAS = [
    "AKT INHIBITOR",
    "PI3K INHIBITOR",
    "MEK INHIBITOR",
    "MTOR INHIBITOR",
    "SCD INHIBITOR",
]


def main() -> None:
    """Run KRAS-only vs 5-gene sensitivity comparison."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("H13j: KRAS-only sensitivity analysis")
    logger.info("=" * 70)

    # Load data
    expression = load_expression(DEPMAP_DIR)
    mutations = load_mutations(DEPMAP_DIR)
    crispr = load_crispr_dependencies(DEPMAP_DIR)
    metadata = load_model_metadata(DEPMAP_DIR)
    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)

    # Run both stratifications
    logger.info("Running 5-gene stratification (KRAS/NRAS/HRAS/BRAF/PIK3CA)...")
    screen_5g, moa_5g, n_t_5g, n_b_5g = _run_screen(
        expression, mutations, crispr, metadata, sensitivity, compound_info,
        ras_genes=list(RAS_PATHWAY_MUTATIONS),
        label="5-gene",
    )

    logger.info("Running 4-gene stratification (KRAS/NRAS/HRAS/BRAF only)...")
    screen_4g, moa_4g, n_t_4g, n_b_4g = _run_screen(
        expression, mutations, crispr, metadata, sensitivity, compound_info,
        ras_genes=KRAS_ONLY_GENES,
        label="4-gene (no PIK3CA)",
    )

    # Build comparison output
    lines: list[str] = []
    lines.append("H13j: KRAS-Only Sensitivity Analysis")
    lines.append("=" * 70)
    lines.append("")
    lines.append("Question: Does the AKT inhibitor signal persist when PIK3CA")
    lines.append("is excluded from the RAS-pathway mutation stratifier?")
    lines.append("")
    lines.append("Stratification comparison:")
    lines.append(f"  5-gene (KRAS/NRAS/HRAS/BRAF/PIK3CA): target={n_t_5g}, background={n_b_5g}")
    lines.append(f"  4-gene (KRAS/NRAS/HRAS/BRAF only):   target={n_t_4g}, background={n_b_4g}")
    lines.append("")
    lines.append("-" * 70)
    lines.append("MOA Enrichment Comparison (Fisher exact test, p < 0.05 threshold)")
    lines.append("-" * 70)
    lines.append("")
    lines.append(f"{'MOA':<25} {'5-gene OR':<12} {'5-gene p':<14} {'4-gene OR':<12} {'4-gene p':<14} {'Status'}")
    lines.append("-" * 85)

    for moa_key in KEY_MOAS:
        row_5g = moa_5g[moa_5g["moa"] == moa_key]
        row_4g = moa_4g[moa_4g["moa"] == moa_key]

        or_5g = float(row_5g["odds_ratio"].iloc[0]) if not row_5g.empty else float("nan")
        p_5g = float(row_5g["fisher_p"].iloc[0]) if not row_5g.empty else float("nan")
        n_sig_5g = int(row_5g["n_sig"].iloc[0]) if not row_5g.empty else 0
        n_tot_5g = int(row_5g["n_total"].iloc[0]) if not row_5g.empty else 0

        or_4g = float(row_4g["odds_ratio"].iloc[0]) if not row_4g.empty else float("nan")
        p_4g = float(row_4g["fisher_p"].iloc[0]) if not row_4g.empty else float("nan")
        n_sig_4g = int(row_4g["n_sig"].iloc[0]) if not row_4g.empty else 0
        n_tot_4g = int(row_4g["n_total"].iloc[0]) if not row_4g.empty else 0

        status = "SURVIVES" if p_4g < 0.05 else "ATTENUATED" if p_4g < 0.1 else "LOST"
        lines.append(
            f"{moa_key:<25} "
            f"{or_5g:<12.2f} {p_5g:<14.2e} "
            f"{or_4g:<12.2f} {p_4g:<14.2e} "
            f"{status}"
        )
        lines.append(
            f"  {'':25} ({n_sig_5g}/{n_tot_5g} sig)     "
            f"({n_sig_4g}/{n_tot_4g} sig)"
        )

    lines.append("")
    lines.append("-" * 70)
    lines.append("Full top-20 MOA classes (4-gene stratification):")
    lines.append("-" * 70)
    lines.append(f"{'MOA':<45} {'Sig/Tot':<10} {'OR':<10} {'p':<14}")
    for _, row in moa_4g.head(20).iterrows():
        lines.append(
            f"{str(row['moa']):<45} "
            f"{row['n_sig']}/{row['n_total']:<7} "
            f"{row['odds_ratio']:<10.2f} "
            f"{row['fisher_p']:<14.2e} {_sig(float(row['fisher_p']))}"
        )

    lines.append("")
    lines.append("=" * 70)
    lines.append("CONCLUSION:")
    # We'll write a placeholder; the actual conclusion depends on results
    akt_row = moa_4g[moa_4g["moa"] == "AKT INHIBITOR"]
    if not akt_row.empty:
        akt_p = float(akt_row["fisher_p"].iloc[0])
        akt_or = float(akt_row["odds_ratio"].iloc[0])
        if akt_p < 0.05:
            lines.append(
                f"AKT inhibitor enrichment SURVIVES (OR={akt_or:.2f}, p={akt_p:.2e})"
            )
            lines.append(
                "The signal is not dependent on PIK3CA inclusion in the stratifier."
            )
        else:
            lines.append(
                f"AKT inhibitor enrichment attenuated (OR={akt_or:.2f}, p={akt_p:.2e})"
            )
    else:
        lines.append("AKT INHIBITOR MOA class not found in results.")

    output_text = "\n".join(lines)
    output_path = RESULTS_DIR / "h13j_kras_only_sensitivity.txt"
    output_path.write_text(output_text)
    logger.info("Results written to %s", output_path)

    # Also print to stdout for immediate visibility
    print(output_text)


if __name__ == "__main__":
    main()
