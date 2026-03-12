"""H13f Metabolic Exploration: Why does SLC2A8 predict AKT inhibitor sensitivity?

SLC2A8 (GLUT8) is a lysosomal/endosomal glucose transporter. Cells with high
SLC2A8 may be in an altered metabolic state that creates AKT dependence.

Analyses:
  1. Metabolic MOA enrichment in PRISM hits (filter for metabolism-related MOAs)
  2. Metabolic gene co-expression with SLC2A8 (glycolysis, glucose transport,
     lipid synthesis, nutrient scavenging)
  3. SLC2A8 vs other GLUT/SLC family members (metabolic phenotype fingerprint)
  4. Metabolic gene expression → AKT drug sensitivity (does metabolism explain
     the AKT link better than PI3K pathway expression?)
  5. Lysosomal/endosomal trafficking genes (SLC2A8 localizes to lysosomes — does
     endosomal machinery co-express and predict AKT sensitivity?)
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

# ---------------------------------------------------------------------------
# Gene sets for metabolic exploration
# ---------------------------------------------------------------------------

# Glycolysis / glucose metabolism
GLYCOLYSIS_GENES: dict[str, str] = {
    "HK1": "hexokinase 1",
    "HK2": "hexokinase 2 (AKT-regulated)",
    "GPI": "glucose-6-phosphate isomerase",
    "PFKP": "phosphofructokinase, platelet",
    "PFKL": "phosphofructokinase, liver",
    "ALDOA": "aldolase, fructose-bisphosphate A",
    "GAPDH": "glyceraldehyde-3-P dehydrogenase",
    "PKM": "pyruvate kinase M1/2",
    "LDHA": "lactate dehydrogenase A",
    "LDHB": "lactate dehydrogenase B",
    "ENO1": "enolase 1",
    "PDK1": "pyruvate dehydrogenase kinase 1",
}

# Glucose transporters (SLC2A family)
GLUCOSE_TRANSPORTERS: dict[str, str] = {
    "SLC2A1": "GLUT1 — main glucose transporter",
    "SLC2A2": "GLUT2 — liver/beta-cell",
    "SLC2A3": "GLUT3 — neuronal, high affinity",
    "SLC2A4": "GLUT4 — insulin-responsive (AKT-regulated!)",
    "SLC2A5": "GLUT5 — fructose",
    "SLC2A6": "GLUT6 — lysosomal",
    "SLC2A8": "GLUT8 — lysosomal/endosomal (our target)",
    "SLC2A9": "GLUT9 — urate",
    "SLC2A10": "GLUT10",
    "SLC2A11": "GLUT11",
    "SLC2A12": "GLUT12",
}

# Lysosomal/endosomal trafficking (where SLC2A8 lives)
LYSO_TRAFFICKING: dict[str, str] = {
    "LAMP1": "lysosomal membrane protein 1",
    "LAMP2": "lysosomal membrane protein 2",
    "RAB7A": "Rab7 — late endosome/lysosome",
    "RAB5A": "Rab5 — early endosome",
    "VPS35": "retromer complex",
    "VPS26A": "retromer complex",
    "VPS29": "retromer complex",
    "SNX1": "sorting nexin 1",
    "MCOLN1": "mucolipin-1 — lysosomal Ca2+ channel",
    "ATP6V0A1": "V-ATPase — lysosomal acidification",
    "ATP6V1A": "V-ATPase subunit A",
    "CTSD": "cathepsin D",
    "TFEB": "master TF for lysosome biogenesis",
    "TMEM55B": "lysosomal membrane protein",
}

# Nutrient scavenging / macropinocytosis / amino acid sensing
NUTRIENT_SENSING: dict[str, str] = {
    "SLC7A5": "LAT1 — amino acid transporter (mTOR activator)",
    "SLC3A2": "4F2hc — LAT1 partner",
    "SLC1A5": "ASCT2 — glutamine transporter",
    "RRAGC": "Rag GTPase C — mTOR amino acid sensing",
    "RRAGA": "Rag GTPase A — mTOR nutrient sensing",
    "LAMTOR1": "Ragulator complex — lysosomal mTOR",
    "LAMTOR2": "Ragulator complex",
    "FLCN": "folliculin — Rag GTPase GAP",
    "FNIP1": "folliculin interacting protein",
    "CASTOR1": "arginine sensor for mTORC1",
    "SESN2": "sestrin 2 — leucine sensor",
    "SLC38A9": "lysosomal arginine transporter → mTOR",
}

# Lipid synthesis (stearoyl-CoA desaturase was enriched in H13c MOA)
LIPID_GENES: dict[str, str] = {
    "SCD": "stearoyl-CoA desaturase (SCD1)",
    "SCD5": "stearoyl-CoA desaturase 5",
    "FASN": "fatty acid synthase",
    "ACLY": "ATP citrate lyase",
    "ACACA": "acetyl-CoA carboxylase alpha",
    "HMGCR": "HMG-CoA reductase (cholesterol)",
    "SREBF1": "SREBP1 — lipogenesis TF",
    "SREBF2": "SREBP2 — cholesterol TF",
}

# AKT inhibitors for drug correlation
AKT_DRUGS: list[str] = ["GSK2110183", "GDC-0068", "AZD5363"]


def _load_all() -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame,
    list[str], list[str],
]:
    """Load data, return expression, sensitivity, compound_info, target_ids, bg_ids."""
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
    border = "=" * 90
    logger.info(border)
    logger.info(title)
    logger.info(border)


# =========================================================================
# ANALYSIS 1: SLC2A8 co-expression with metabolic genes
# =========================================================================

def analysis_1_metabolic_coexpression(
    expression: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Which metabolic genes co-express with SLC2A8 and are elevated in target?"""
    _header("ANALYSIS 1: SLC2A8 METABOLIC CO-EXPRESSION FINGERPRINT")

    if "SLC2A8" not in expression.columns:
        logger.info("  SLC2A8 not found")
        return

    slc2a8 = expression["SLC2A8"]

    gene_groups = {
        "GLYCOLYSIS": GLYCOLYSIS_GENES,
        "GLUCOSE TRANSPORTERS": GLUCOSE_TRANSPORTERS,
        "LYSOSOME/ENDOSOME TRAFFICKING": LYSO_TRAFFICKING,
        "NUTRIENT SENSING / mTOR": NUTRIENT_SENSING,
        "LIPID SYNTHESIS": LIPID_GENES,
    }

    for group_name, genes in gene_groups.items():
        logger.info(f"\n  --- {group_name} ---")
        logger.info(
            f"  {'Gene':<12} {'SLC2A8 rho':>12} {'rho p':>12} "
            f"{'Target d':>10} {'diff p':>12}"
        )
        logger.info("  " + "-" * 70)

        for gene, desc in genes.items():
            if gene not in expression.columns:
                continue

            # Correlation with SLC2A8
            shared = sorted(
                set(slc2a8.dropna().index) & set(expression[gene].dropna().index)
            )
            if len(shared) < 50:
                continue
            rho, p_rho = stats.spearmanr(slc2a8[shared], expression.loc[shared, gene])

            # Differential expression
            t_vals = expression.loc[
                expression.index.isin(target_ids), gene
            ].dropna()
            b_vals = expression.loc[
                expression.index.isin(background_ids), gene
            ].dropna()

            if len(t_vals) < 5 or len(b_vals) < 5:
                continue

            _, p_diff = stats.ttest_ind(t_vals, b_vals)
            d = (t_vals.mean() - b_vals.mean()) / np.sqrt(
                (t_vals.std() ** 2 + b_vals.std() ** 2) / 2
            )

            sig_rho = _sig(p_rho)
            sig_diff = _sig(p_diff)
            logger.info(
                f"  {gene:<12} {rho:>+12.4f} {p_rho:>12.2e} "
                f"{d:>+10.3f} {p_diff:>12.2e} {sig_rho} {sig_diff}"
            )


# =========================================================================
# ANALYSIS 2: Metabolic genes → AKT drug sensitivity
# =========================================================================

def analysis_2_metabolic_drug_link(
    expression: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
) -> None:
    """Do metabolic genes predict AKT drug sensitivity better than PI3K genes?"""
    _header("ANALYSIS 2: METABOLIC GENES → AKT DRUG SENSITIVITY")

    # Candidate metabolic predictors
    candidates: dict[str, str] = {
        "SLC2A8": "GLUT8 (our biomarker)",
        "HK2": "hexokinase 2 (AKT-regulated)",
        "SLC2A4": "GLUT4 (AKT-regulated)",
        "SLC2A6": "GLUT6 (lysosomal — same compartment)",
        "SLC2A1": "GLUT1 (main glucose transporter)",
        "SLC38A9": "lysosomal arginine → mTOR",
        "LAMTOR1": "Ragulator (lysosomal mTOR)",
        "SCD": "SCD1 (stearoyl-CoA desaturase)",
        "FASN": "fatty acid synthase",
        "TFEB": "lysosome biogenesis TF",
        "MCOLN1": "mucolipin-1 (lysosomal Ca2+)",
        "LAMP1": "lysosomal membrane",
        "AKT1": "AKT1 (control — PI3K pathway)",
        "MTOR": "mTOR (control — PI3K pathway)",
    }

    logger.info("  Spearman correlation: gene expression vs AKT drug sensitivity")
    logger.info("  (negative rho = higher expression → more killing)")
    logger.info("")

    for drug_name in AKT_DRUGS:
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            logger.info(f"  {drug_name}: NOT FOUND")
            continue

        logger.info(f"  {drug_name}:")
        logger.info(f"  {'Gene':<12} {'Description':<35} {'rho':>8} {'p':>12}")
        logger.info("  " + "-" * 75)

        results = []
        for gene, desc in candidates.items():
            if gene not in expression.columns:
                continue
            shared = sorted(
                set(expression[gene].dropna().index) & set(sens.index)
            )
            if len(shared) < 30:
                continue
            rho, p = stats.spearmanr(expression.loc[shared, gene], sens[shared])
            results.append((gene, desc, rho, p))

        # Sort by absolute rho
        results.sort(key=lambda x: abs(x[2]), reverse=True)
        for gene, desc, rho, p in results:
            logger.info(
                f"  {gene:<12} {desc:<35} {rho:>+8.4f} {p:>12.2e} {_sig(p)}"
            )
        logger.info("")


# =========================================================================
# ANALYSIS 3: Lysosomal metabolic z-score vs AKT sensitivity
# =========================================================================

def analysis_3_lysosomal_metabolic_score(
    expression: pd.DataFrame,
    sensitivity: pd.DataFrame,
    compound_info: pd.DataFrame,
    name_col: str,
    id_col: str,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Compute a 'lysosomal metabolic' z-score and test vs AKT drugs."""
    _header("ANALYSIS 3: LYSOSOMAL METABOLIC SCORE vs AKT DRUGS")

    # Lysosomal metabolic signature: SLC2A8 + trafficking + nutrient sensing
    lyso_metab_genes = [
        "SLC2A8", "SLC2A6",  # lysosomal glucose transporters
        "LAMP1", "LAMP2",  # lysosomal membrane
        "MCOLN1",  # lysosomal Ca2+
        "ATP6V0A1", "ATP6V1A",  # V-ATPase (acidification)
        "SLC38A9",  # lysosomal arginine → mTOR
        "LAMTOR1", "LAMTOR2",  # ragulator (lysosomal mTOR platform)
        "RAB7A",  # late endosome
        "CTSD",  # cathepsin D
    ]

    available = [g for g in lyso_metab_genes if g in expression.columns]
    logger.info(f"  Lysosomal metabolic genes: {available}")

    subset = expression[available]
    zscored = subset.apply(stats.zscore, nan_policy="omit")
    lyso_score = zscored.mean(axis=1)

    # Differential in target vs background
    t_score = lyso_score[lyso_score.index.isin(target_ids)].dropna()
    b_score = lyso_score[lyso_score.index.isin(background_ids)].dropna()
    _, p_diff = stats.ttest_ind(t_score, b_score)
    d = (t_score.mean() - b_score.mean()) / np.sqrt(
        (t_score.std() ** 2 + b_score.std() ** 2) / 2
    )
    logger.info(f"  Target mean: {t_score.mean():.4f}, BG mean: {b_score.mean():.4f}")
    logger.info(f"  Cohen's d: {d:+.4f}, p: {p_diff:.4e} {_sig(p_diff)}")
    logger.info("")

    # Correlation with AKT drugs
    logger.info("  Lysosomal metabolic score vs AKT drugs:")
    for drug_name in AKT_DRUGS:
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            continue
        shared = sorted(set(lyso_score.dropna().index) & set(sens.index))
        if len(shared) < 30:
            continue
        rho, p = stats.spearmanr(lyso_score[shared], sens[shared])
        logger.info(f"    {drug_name:<20} rho={rho:+.4f}, p={p:.4e} {_sig(p)}")

    # Compare: SLC2A8 alone vs lysosomal score vs PI3K score
    logger.info("")
    logger.info("  HEAD-TO-HEAD: SLC2A8 alone vs lyso score vs PI3K score → AKT drugs")

    pi3k_genes = [g for g in [
        "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R2",
        "AKT1", "AKT2", "MTOR", "RPTOR", "RPS6KB1",
    ] if g in expression.columns]
    pi3k_z = expression[pi3k_genes].apply(stats.zscore, nan_policy="omit").mean(axis=1)

    slc2a8 = expression["SLC2A8"] if "SLC2A8" in expression.columns else None

    for drug_name in AKT_DRUGS:
        sens = _resolve_drug(drug_name, sensitivity, compound_info, name_col, id_col)
        if sens is None:
            continue

        logger.info(f"  {drug_name}:")
        for name, score in [
            ("SLC2A8 alone", slc2a8),
            ("Lyso metabolic", lyso_score),
            ("PI3K pathway", pi3k_z),
        ]:
            if score is None:
                continue
            shared = sorted(set(score.dropna().index) & set(sens.index))
            if len(shared) < 30:
                continue
            rho, p = stats.spearmanr(score[shared], sens[shared])
            logger.info(f"    {name:<20} rho={rho:+.4f}, p={p:.4e} {_sig(p)}")


# =========================================================================
# ANALYSIS 4: CRISPR dependency for lysosomal/metabolic genes
# =========================================================================

def analysis_4_metabolic_dependencies(
    crispr: pd.DataFrame,
    target_ids: list[str],
    background_ids: list[str],
) -> None:
    """Are target cells more dependent on lysosomal/metabolic genes?"""
    _header("ANALYSIS 4: METABOLIC GENE DEPENDENCIES (CRISPR)")

    test_genes: dict[str, str] = {
        "HK2": "hexokinase 2 (AKT-phosphorylated)",
        "SLC2A1": "GLUT1",
        "LAMTOR1": "Ragulator 1 (lysosomal mTOR)",
        "LAMTOR2": "Ragulator 2",
        "RRAGA": "Rag GTPase A (nutrient sensing)",
        "RRAGC": "Rag GTPase C",
        "SLC38A9": "lysosomal arginine transporter",
        "ATP6V1A": "V-ATPase (lysosomal acidification)",
        "RAB7A": "late endosome/lysosome Rab",
        "VPS35": "retromer (endosomal recycling)",
        "MCOLN1": "mucolipin-1",
        "TFEB": "lysosome biogenesis TF",
        "SCD": "stearoyl-CoA desaturase",
        "FASN": "fatty acid synthase",
        "ACLY": "ATP citrate lyase",
        "PIK3C3": "VPS34 (control — known hit)",
        "AKT1": "AKT1 (control)",
    }

    logger.info(f"  {'Gene':<12} {'Description':<40} {'Target dep':>12} {'BG dep':>10} {'d':>8} {'p':>12}")
    logger.info("  " + "-" * 100)

    for gene, desc in test_genes.items():
        if gene not in crispr.columns:
            continue

        t_vals = crispr.loc[crispr.index.isin(target_ids), gene].dropna()
        b_vals = crispr.loc[crispr.index.isin(background_ids), gene].dropna()

        if len(t_vals) < 5 or len(b_vals) < 5:
            continue

        _, p = stats.ttest_ind(t_vals, b_vals)
        d = (t_vals.mean() - b_vals.mean()) / np.sqrt(
            (t_vals.std() ** 2 + b_vals.std() ** 2) / 2
        )
        # Negative dependency = essential
        logger.info(
            f"  {gene:<12} {desc:<40} {t_vals.mean():>12.4f} {b_vals.mean():>10.4f} "
            f"{d:>+8.3f} {p:>12.2e} {_sig(p)}"
        )


# =========================================================================
# Main
# =========================================================================

def main() -> None:
    """Run metabolic exploration."""
    _header("H13f: METABOLIC MECHANISM EXPLORATION")

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
    logger.info("Data loaded")

    analysis_1_metabolic_coexpression(expression, target_ids, background_ids)
    analysis_2_metabolic_drug_link(
        expression, sensitivity, compound_info, name_col, id_col,
    )
    analysis_3_lysosomal_metabolic_score(
        expression, sensitivity, compound_info, name_col, id_col,
        target_ids, background_ids,
    )
    analysis_4_metabolic_dependencies(crispr, target_ids, background_ids)

    _header("H13f METABOLIC EXPLORATION COMPLETE")


if __name__ == "__main__":
    main()
