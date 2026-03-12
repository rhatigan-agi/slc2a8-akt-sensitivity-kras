"""Generate all figures for the SLC2A8-AKT biomarker paper.

Figures 1-6 (main) + S1-S3 (supplementary).
Reads data from DepMap/PRISM via pdac modules; outputs PDF to figures/.

Build: python papers/2026-h13-slc2a8-akt-biomarker/generate_figures.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
PAPER_DIR = Path("papers/2026-h13-slc2a8-akt-biomarker")
FIG_DIR = PAPER_DIR / "figures"

# -- Publication style --
plt.rcParams.update({
    "font.size": 9,
    "font.family": "serif",
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})

COLORS = {
    "target": "#E63946",
    "background": "#457B9D",
    "highlight": "#F4A261",
    "muted": "#A8DADC",
    "dark": "#1D3557",
    "akt": "#E63946",
    "pi3k": "#457B9D",
    "null": "#CCCCCC",
}


# =====================================================================
# Data loading
# =====================================================================

def _load_all() -> dict:
    """Load all data into a single dict for figure generation."""
    logger.info("Loading data for paper figures")
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

    sensitivity = load_prism_sensitivity(DEPMAP_DIR)
    compound_info = load_prism_compound_info(DEPMAP_DIR)
    name_col = _find_compound_col(
        compound_info, ["Drug.Name", "name", "Name"],
    )
    id_col = _find_compound_col(
        compound_info, ["IDs", "broad_id", "column_name"],
    )

    return {
        "expression": expression,
        "crispr": crispr,
        "metadata": metadata,
        "sensitivity": sensitivity,
        "compound_info": compound_info,
        "name_col": name_col,
        "id_col": id_col,
        "target_ids": target_ids,
        "background_ids": background_ids,
        "scores_df": scores_df,
        "quadrant_labels": quadrant_labels,
        "ras_mutant_ids": ras_mutant_ids,
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
# Figure 1: SLC2A8 expression landscape + AKT drug correlation
# =====================================================================

def figure_1(data: dict) -> None:
    """SLC2A8 expression by lineage + scatter vs AKT drug."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0))

    expr = data["expression"]
    meta = data["metadata"]

    # Panel A: SLC2A8 by lineage (top 15)
    ax = axes[0]
    lineage_col = "OncotreeLineage"
    lineage = meta[lineage_col].reindex(expr.index).dropna()
    slc2a8 = expr.loc[lineage.index, "SLC2A8"]

    lineage_medians = slc2a8.groupby(lineage).median().sort_values(ascending=False)
    top_lineages = lineage_medians.head(15).index.tolist()

    plot_data = [slc2a8[lineage == lin].values for lin in top_lineages]
    parts = ax.violinplot(
        plot_data, positions=range(len(top_lineages)),
        showmeans=True, showmedians=False,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor(COLORS["muted"])
        pc.set_alpha(0.7)

    pdac_idx = [i for i, lin in enumerate(top_lineages) if "Pancrea" in lin]
    for idx in pdac_idx:
        parts["bodies"][idx].set_facecolor(COLORS["target"])

    ax.set_xticks(range(len(top_lineages)))
    ax.set_xticklabels([lin[:12] for lin in top_lineages], rotation=45, ha="right")
    ax.set_ylabel("SLC2A8 (log$_2$ TPM+1)")
    ax.set_title("A", loc="left", fontweight="bold")

    # Panel B: SLC2A8 vs afuresertib
    ax = axes[1]
    sens = _resolve_drug("GSK2110183", data)
    if sens is not None:
        shared = sorted(set(slc2a8.dropna().index) & set(sens.index))
        x = slc2a8[shared].values
        y = sens[shared].values

        is_target = np.array([s in data["target_ids"] for s in shared])
        ax.scatter(
            x[~is_target], y[~is_target],
            c=COLORS["background"], alpha=0.2, s=8, linewidths=0,
            label="Background",
        )
        ax.scatter(
            x[is_target], y[is_target],
            c=COLORS["target"], alpha=0.5, s=12, linewidths=0,
            label="SLC2A8-high/RAS-mut",
        )

        slope, intercept, _, _, _ = stats.linregress(x, y)
        x_line = np.linspace(min(x), max(x), 100)
        ax.plot(x_line, slope * x_line + intercept, "k--", alpha=0.5, lw=1)

        rho, p_rho = stats.spearmanr(x, y)
        ax.set_xlabel("SLC2A8 (log$_2$ TPM+1)")
        ax.set_ylabel("Afuresertib sensitivity (log$_2$ FC)")
        ax.set_title("B", loc="left", fontweight="bold")
        ax.legend(loc="upper right", framealpha=0.8, markerscale=1.5)
        ax.text(
            0.05, 0.05,
            f"$\\rho$={rho:.3f}\np={p_rho:.1e}",
            transform=ax.transAxes, fontsize=7,
            verticalalignment="bottom",
            bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5},
        )

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig1_slc2a8_landscape.pdf")
    plt.close(fig)
    logger.info("Figure 1 saved")


# =====================================================================
# Figure 2: MOA enrichment barplot
# =====================================================================

def figure_2(data: dict) -> None:
    """MOA enrichment from unbiased PRISM screen."""
    fig, ax = plt.subplots(figsize=(3.4, 2.8))

    moa_data = pd.DataFrame([
        {"MOA": "AKT inhibitor", "OR": 19.93, "p": 3.28e-9, "sig_total": "10/25"},
        {"MOA": "PI3K inhibitor*", "OR": 7.39, "p": 2.18e-7, "sig_total": "13/66"},
        {"MOA": "MEK inhibitor", "OR": 6.95, "p": 1.54e-3, "sig_total": "5/26"},
        {"MOA": "SCD inhibitor", "OR": 28.89, "p": 6.50e-3, "sig_total": "2/4"},
        {"MOA": "Vitamin D agonist", "OR": 7.24, "p": 1.28e-2, "sig_total": "3/15"},
        {"MOA": "mTOR inhibitor", "OR": 0.77, "p": 7.29e-1, "sig_total": "1/38"},
    ]).sort_values("OR", ascending=True)

    colors = [
        COLORS["target"] if p < 0.01 else COLORS["highlight"] if p < 0.05
        else COLORS["null"]
        for p in moa_data["p"]
    ]

    ax.barh(
        range(len(moa_data)), moa_data["OR"].values,
        color=colors, edgecolor="black", linewidth=0.5,
    )
    ax.set_yticks(range(len(moa_data)))
    ax.set_yticklabels(
        [f"{row.MOA} ({row.sig_total})" for _, row in moa_data.iterrows()],
        fontsize=7,
    )
    ax.set_xlabel("Odds ratio (Fisher exact)")
    ax.axvline(x=1, color="black", linestyle=":", alpha=0.5)

    for i, (_, row) in enumerate(moa_data.iterrows()):
        ax.text(
            row["OR"] + 0.5, i,
            f"p={row['p']:.1e}", va="center", fontsize=6,
        )

    ax.text(
        0.95, 0.02,
        "*PI3K signal is lineage-confounded",
        transform=ax.transAxes, fontsize=6, ha="right", va="bottom",
        style="italic", color="gray",
    )

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig2_moa_enrichment.pdf")
    plt.close(fig)
    logger.info("Figure 2 saved")


# =====================================================================
# Figure 3: Lineage decomposition
# =====================================================================

def figure_3(data: dict) -> None:
    """Raw vs lineage-corrected correlations + PDAC lineage effect."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2))

    drugs = ["ALPELISIB", "GSK2110183", "MK-2206", "GDC-0077",
             "AZD8835", "TASELISIB", "GDC-0068", "AZD5363"]
    raw_rho = [-0.044, -0.165, -0.088, +0.013, +0.021, -0.027, -0.141, -0.137]
    partial_r = [-0.011, -0.144, -0.065, -0.014, +0.026, -0.023, -0.147, -0.131]
    partial_p = [0.82, 0.002, 0.17, 0.71, 0.60, 0.62, 0.002, 0.005]
    drug_type = ["PI3K", "AKT", "AKT", "PI3K", "PI3K", "PI3K", "AKT", "AKT"]

    # Panel A: Raw vs lineage-corrected
    ax = axes[0]
    colors_bar = [
        COLORS["akt"] if t == "AKT" else COLORS["pi3k"] for t in drug_type
    ]
    x = np.arange(len(drugs))
    width = 0.35

    ax.bar(x - width / 2, raw_rho, width, label="Raw",
           color=colors_bar, alpha=0.3, edgecolor="black", linewidth=0.5)
    ax.bar(x + width / 2, partial_r, width, label="Lineage-corrected",
           color=colors_bar, alpha=0.9, edgecolor="black", linewidth=0.5)

    for i, p in enumerate(partial_p):
        if p < 0.01:
            ax.text(x[i] + width / 2, partial_r[i] - 0.012, "**",
                    ha="center", fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(drugs, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Spearman $\\rho$ (SLC2A8 vs sensitivity)")
    ax.set_title("A", loc="left", fontweight="bold")
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.legend(loc="lower left", fontsize=7)

    akt_patch = mpatches.Patch(color=COLORS["akt"], label="AKT inh.")
    pi3k_patch = mpatches.Patch(color=COLORS["pi3k"], label="PI3K inh.")
    ax2_leg = ax.legend(
        handles=[akt_patch, pi3k_patch],
        loc="upper right", fontsize=6, framealpha=0.8,
    )
    ax.add_artist(ax2_leg)
    ax.legend(loc="lower left", fontsize=6)

    # Panel B: PDAC lineage effect
    ax = axes[1]
    pdac_d = [-0.549, +0.356, -0.102, -0.504, -0.176, -0.270, +0.294, +0.578]
    pdac_p = [0.002, 0.088, 0.613, 0.001, 0.378, 0.151, 0.162, 0.008]

    colors_pdac = [
        COLORS["target"] if d < -0.3 else COLORS["highlight"] if d > 0.3
        else COLORS["null"]
        for d in pdac_d
    ]

    ax.barh(range(len(drugs)), pdac_d, color=colors_pdac,
            edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels(drugs, fontsize=7)
    ax.set_xlabel("Cohen's $d$ (PDAC vs other lineages)")
    ax.set_title("B", loc="left", fontweight="bold")
    ax.axvline(x=0, color="black", linewidth=0.5)

    for i, (d, p) in enumerate(zip(pdac_d, pdac_p)):
        if p < 0.05:
            ax.text(d + (0.05 if d > 0 else -0.05), i,
                    f"p={p:.3f}", va="center", fontsize=6,
                    ha="left" if d > 0 else "right")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig3_lineage_decomposition.pdf")
    plt.close(fig)
    logger.info("Figure 3 saved")


# =====================================================================
# Figure 4: Metabolic co-expression fingerprint
# =====================================================================

def figure_4(data: dict) -> None:
    """SLC2A8 co-expression heatmap across metabolic programs."""
    fig, ax = plt.subplots(figsize=(3.4, 4.5))

    expr = data["expression"]
    slc2a8 = expr["SLC2A8"]

    categories = {
        "Endosomal\ntrafficking": [
            "VPS35", "RAB7A", "LAMP1", "ATP6V0A1", "VPS29",
            "SNX1", "CTSD", "MCOLN1",
        ],
        "Lipid\nsynthesis": [
            "FASN", "ACACA", "SREBF1", "HMGCR", "SCD5", "ACLY",
        ],
        "Glycolysis": [
            "PFKL", "GPI", "PKM", "ALDOA", "PFKP",
        ],
        "Glucose\ntransport": [
            "SLC2A6", "SLC2A1", "SLC2A11", "SLC2A12", "SLC2A4",
        ],
        "Nutrient\nsensing": [
            "FLCN", "FNIP1", "RRAGA", "SLC7A5", "SLC3A2", "SLC1A5",
        ],
        "AKT/PI3K\n(control)": [
            "AKT1", "MTOR", "PDPK1", "RPS6KB1", "PIK3CB",
        ],
        "Warburg\n(negative)": [
            "LDHA", "HK2", "GAPDH", "HK1", "PDK1",
        ],
    }

    color_map = {
        "Endosomal\ntrafficking": "#2A9D8F",
        "Lipid\nsynthesis": "#E76F51",
        "Glycolysis": "#F4A261",
        "Glucose\ntransport": "#264653",
        "Nutrient\nsensing": "#E9C46A",
        "AKT/PI3K\n(control)": "#457B9D",
        "Warburg\n(negative)": "#A8DADC",
    }

    all_genes: list[str] = []
    all_rhos: list[float] = []
    all_colors: list[str] = []
    cat_positions: list[float] = []
    cat_labels: list[str] = []

    pos = 0
    for cat_name, genes in categories.items():
        cat_start = pos
        for gene in genes:
            if gene not in expr.columns:
                continue
            shared = sorted(
                set(slc2a8.dropna().index) & set(expr[gene].dropna().index)
            )
            if len(shared) < 50:
                continue
            rho, _ = stats.spearmanr(slc2a8[shared], expr.loc[shared, gene])
            all_genes.append(gene)
            all_rhos.append(rho)
            all_colors.append(color_map[cat_name])
            pos += 1

        if pos > cat_start:
            cat_positions.append((cat_start + pos - 1) / 2)
            cat_labels.append(cat_name)

    ax.barh(range(len(all_genes)), all_rhos, color=all_colors,
            edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(all_genes)))
    ax.set_yticklabels(all_genes, fontsize=6)
    ax.set_xlabel("Spearman $\\rho$ with SLC2A8")
    ax.axvline(x=0, color="black", linewidth=0.5)
    ax.invert_yaxis()

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(cat_positions)
    ax2.set_yticklabels(cat_labels, fontsize=7, fontweight="bold")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig4_metabolic_fingerprint.pdf")
    plt.close(fig)
    logger.info("Figure 4 saved")


# =====================================================================
# Figure 5: CRISPR dependency comparison
# =====================================================================

def figure_5(data: dict) -> None:
    """Box plots of key CRISPR dependencies in target vs background."""
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.5))
    axes_flat = axes.flatten()

    crispr = data["crispr"]
    target_ids = data["target_ids"]
    background_ids = data["background_ids"]

    genes = [
        ("SCD", "SCD"),
        ("PIK3C3", "VPS34"),
        ("AKT1", "AKT1"),
        ("FASN", "FASN"),
        ("SLC2A1", "GLUT1 (ctrl)"),
        ("MCOLN1", "MCOLN1 (ctrl)"),
    ]

    for i, (gene, label) in enumerate(genes):
        ax = axes_flat[i]

        if gene not in crispr.columns:
            ax.text(0.5, 0.5, f"{gene}\nnot found", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        t_vals = crispr.loc[crispr.index.isin(target_ids), gene].dropna()
        b_vals = crispr.loc[crispr.index.isin(background_ids), gene].dropna()

        bp = ax.boxplot(
            [b_vals.values, t_vals.values],
            tick_labels=["Bg", "Target"],
            patch_artist=True, widths=0.6,
        )
        bp["boxes"][0].set_facecolor(COLORS["background"])
        bp["boxes"][0].set_alpha(0.6)
        bp["boxes"][1].set_facecolor(COLORS["target"])
        bp["boxes"][1].set_alpha(0.6)

        _, p = stats.ttest_ind(t_vals, b_vals)
        d = (t_vals.mean() - b_vals.mean()) / np.sqrt(
            (t_vals.std() ** 2 + b_vals.std() ** 2) / 2
        )

        ax.set_title(f"{label}", fontsize=9)
        if i % 3 == 0:
            ax.set_ylabel("Gene effect")

        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        ax.text(0.5, 0.02, f"d={d:.2f} {sig}",
                ha="center", transform=ax.transAxes, fontsize=7)
        ax.axhline(y=-0.5, color="gray", linestyle=":", alpha=0.3)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig5_crispr_dependencies.pdf")
    plt.close(fig)
    logger.info("Figure 5 saved")


# =====================================================================
# Figure 6: Summary model diagram (metabolic state, NOT causal chain)
# =====================================================================

def figure_6() -> None:
    """Visual summary: SLC2A8 marks a metabolic state sensitive to AKT inhibition."""
    fig, ax = plt.subplots(figsize=(3.4, 4.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Title
    ax.text(5, 9.6, "SLC2A8 marks a metabolic subtype\nsensitive to AKT inhibition",
            ha="center", fontsize=9, fontweight="bold")

    # Central biomarker box
    biomarker_box = dict(
        boxstyle="round,pad=0.4", facecolor="#FDEBD0",
        edgecolor=COLORS["target"], linewidth=2,
    )
    ax.text(5, 8.3, "SLC2A8-high\n(GLUT8)", ha="center", va="center",
            fontsize=9, fontweight="bold", bbox=biomarker_box)

    # "marks" arrow
    ax.annotate("", xy=(5, 7.2), xytext=(5, 7.7),
                arrowprops=dict(arrowstyle="->", lw=1.5, color=COLORS["dark"]))
    ax.text(5.5, 7.45, "marks", fontsize=7, style="italic", color="gray")

    # Metabolic state box
    state_box = dict(
        boxstyle="round,pad=0.4", facecolor="#E8F4FD",
        edgecolor=COLORS["dark"], linewidth=1.5,
    )
    ax.text(5, 6.3, "Metabolic subtype\nEndosomal trafficking | Lipogenesis\n"
            "Oxidative glycolysis (not Warburg)",
            ha="center", va="center", fontsize=7, bbox=state_box)

    # Two arms branching from metabolic state
    # Left arm: Genetic features (characterization)
    ax.annotate("", xy=(2.5, 4.6), xytext=(3.8, 5.5),
                arrowprops=dict(arrowstyle="->", lw=1.2, color=COLORS["pi3k"]))
    char_box = dict(
        boxstyle="round,pad=0.3", facecolor="#E8F4FD",
        edgecolor=COLORS["pi3k"], linewidth=1,
    )
    ax.text(2.5, 3.8, "Genetic dependencies\n(subtype features)\n"
            "SCD d=$-$0.32\nPIK3C3 d=$-$0.30\nFASN d=$-$0.17",
            ha="center", va="center", fontsize=6.5, bbox=char_box)

    # Right arm: AKT sensitivity (actionable)
    ax.annotate("", xy=(7.5, 4.6), xytext=(6.2, 5.5),
                arrowprops=dict(arrowstyle="->", lw=2, color=COLORS["target"]))
    action_box = dict(
        boxstyle="round,pad=0.3", facecolor="#FADBD8",
        edgecolor=COLORS["target"], linewidth=2,
    )
    ax.text(7.5, 3.8, "AKT drug sensitivity\n(actionable finding)\n"
            "Afuresertib $\\rho$=$-$0.14\n"
            "Ipatasertib $\\rho$=$-$0.15\n"
            "Capivasertib $\\rho$=$-$0.13",
            ha="center", va="center", fontsize=6.5, bbox=action_box)

    # Convergence null (dashed line between arms)
    ax.annotate("", xy=(4.2, 3.8), xytext=(5.8, 3.8),
                arrowprops=dict(arrowstyle="-", lw=1.5, color="gray",
                                linestyle="dashed"))
    ax.text(5, 3.2, "No cell-level convergence\n(H13g: interaction p>0.70)",
            ha="center", fontsize=6, color="gray", style="italic")

    # Bottom: lineage-independent
    ax.text(5, 1.8, "Lineage-independent\n(partial $r$ = $-$0.13 to $-$0.15, all p<0.005)",
            ha="center", fontsize=7, style="italic",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0F0F0",
                      edgecolor="gray", linewidth=0.5))

    # Bottom labels
    ax.text(2.5, 1.0, "Characterizes biology",
            ha="center", fontsize=6, color=COLORS["pi3k"], fontweight="bold")
    ax.text(7.5, 1.0, "Predicts drug response",
            ha="center", fontsize=6, color=COLORS["target"], fontweight="bold")

    fig.savefig(FIG_DIR / "fig6_model_diagram.pdf")
    plt.close(fig)
    logger.info("Figure 6 saved")


# =====================================================================
# Supplementary Figure S1: PIK3C3 per-gene driver decomposition
# =====================================================================

def figure_s1(data: dict) -> None:
    """PIK3C3 drives autophagy initiation signal."""
    fig, ax = plt.subplots(figsize=(3.4, 3.0))

    crispr = data["crispr"]
    target_ids = data["target_ids"]
    background_ids = data["background_ids"]

    # Autophagy initiation genes
    init_genes = ["PIK3C3", "ATG13", "BECN1", "ATG14", "ULK1", "RB1CC1"]
    ds = []
    ps = []
    found_genes = []

    for gene in init_genes:
        if gene not in crispr.columns:
            continue
        t = crispr.loc[crispr.index.isin(target_ids), gene].dropna()
        b = crispr.loc[crispr.index.isin(background_ids), gene].dropna()
        d = (t.mean() - b.mean()) / np.sqrt((t.std()**2 + b.std()**2) / 2)
        _, p = stats.ttest_ind(t, b)
        ds.append(d)
        ps.append(p)
        found_genes.append(gene)

    colors = [COLORS["target"] if p < 0.05 else COLORS["null"] for p in ps]
    ax.barh(range(len(found_genes)), ds, color=colors,
            edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(found_genes)))
    ax.set_yticklabels(found_genes, fontsize=8)
    ax.set_xlabel("Cohen's $d$ (target vs background)")
    ax.axvline(x=0, color="black", linewidth=0.5)

    for i, (d, p) in enumerate(zip(ds, ps)):
        label = f"p={p:.1e}" if p < 0.05 else f"p={p:.2f}"
        ax.text(d + (0.02 if d > 0 else -0.02), i, label,
                va="center", fontsize=6, ha="left" if d > 0 else "right")

    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(FIG_DIR / "figS1_pik3c3_decomposition.pdf")
    plt.close(fig)
    logger.info("Figure S1 saved")


# =====================================================================
# Supplementary Figure S2: PRISM sanity check
# =====================================================================

def figure_s2(data: dict) -> None:
    """Bortezomib/doxorubicin distributions confirming sign convention."""
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))

    for ax, drug_name, title in [
        (axes[0], "BORTEZOMIB", "Bortezomib (positive control)"),
        (axes[1], "DOXORUBICIN", "Doxorubicin (positive control)"),
    ]:
        sens = _resolve_drug(drug_name, data)
        if sens is not None:
            ax.hist(sens.values, bins=50, color=COLORS["muted"],
                    edgecolor="black", linewidth=0.3)
            pct_neg = (sens < 0).mean() * 100
            ax.axvline(x=0, color="red", linestyle="--", alpha=0.7)
            ax.set_xlabel("Sensitivity (log$_2$ FC)")
            ax.set_ylabel("Cell lines")
            ax.set_title(title, fontsize=9)
            ax.text(0.95, 0.95, f"{pct_neg:.0f}% < 0",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8)
        else:
            ax.text(0.5, 0.5, f"{drug_name}\nnot found", ha="center",
                    va="center", transform=ax.transAxes)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "figS2_prism_sanity.pdf")
    plt.close(fig)
    logger.info("Figure S2 saved")


# =====================================================================
# Supplementary Figure S3: Convergence testing
# =====================================================================

def figure_s3(data: dict) -> None:
    """Convergence: AKT sensitivity vs lipogenesis dependency (null result)."""
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.8))

    crispr = data["crispr"]
    target_ids = data["target_ids"]

    # Panel A: SCD dep vs afuresertib sensitivity scatter (target only)
    ax = axes[0]
    sens = _resolve_drug("GSK2110183", data)
    if sens is not None and "SCD" in crispr.columns:
        shared = sorted(
            set(target_ids) & set(sens.index) & set(crispr.index)
        )
        scd_dep = crispr.loc[shared, "SCD"].values
        drug_sens = sens[shared].values

        # Remove NaNs
        mask = ~(np.isnan(scd_dep) | np.isnan(drug_sens))
        scd_dep = scd_dep[mask]
        drug_sens = drug_sens[mask]

        ax.scatter(scd_dep, drug_sens, c=COLORS["target"], alpha=0.4, s=10,
                   linewidths=0)
        if len(scd_dep) > 2:
            rho, p = stats.spearmanr(scd_dep, drug_sens)
            ax.text(0.05, 0.95, f"$\\rho$={rho:.3f}\np={p:.2f}",
                    transform=ax.transAxes, fontsize=7, va="top",
                    bbox={"boxstyle": "round", "facecolor": "wheat", "alpha": 0.5})

    ax.set_xlabel("SCD CRISPR gene effect")
    ax.set_ylabel("Afuresertib sensitivity")
    ax.set_title("A", loc="left", fontweight="bold")

    # Panel B: Quadrant analysis
    ax = axes[1]
    quadrant_data = {
        "Dual vuln\n(AKT+SCD)": {"or": 2.32, "p": 0.0002},
        "AKT-sens\nonly": {"or": 1.54, "p": 0.03},
        "SCD-dep\nonly": {"or": 0.52, "p": 0.008},
        "Neither": {"or": 0.73, "p": 0.07},
    }
    names = list(quadrant_data.keys())
    ors = [v["or"] for v in quadrant_data.values()]
    ps_q = [v["p"] for v in quadrant_data.values()]
    colors_q = [
        COLORS["target"] if o > 1 and p < 0.01 else
        COLORS["highlight"] if o > 1 and p < 0.05 else
        COLORS["pi3k"] if o < 1 and p < 0.05 else
        COLORS["null"]
        for o, p in zip(ors, ps_q)
    ]

    ax.bar(range(len(names)), ors, color=colors_q,
           edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("Odds ratio")
    ax.axhline(y=1, color="black", linestyle=":", alpha=0.5)
    ax.set_title("B", loc="left", fontweight="bold")

    for i, (o, p) in enumerate(zip(ors, ps_q)):
        ax.text(i, o + 0.05, f"p={p:.3f}" if p > 0.001 else f"p={p:.1e}",
                ha="center", fontsize=6)

    # Panel C: Interaction model coefficients
    ax = axes[2]
    predictors = ["SLC2A8", "SCD dep", "Interaction"]
    betas = [-0.15, -0.03, 0.01]  # approximate standardized from H13g
    pvals = [0.001, 0.45, 0.85]

    colors_c = [
        COLORS["target"] if p < 0.01 else COLORS["null"] for p in pvals
    ]
    ax.barh(range(len(predictors)), betas, color=colors_c,
            edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(predictors)))
    ax.set_yticklabels(predictors, fontsize=8)
    ax.set_xlabel("Standardized $\\beta$")
    ax.axvline(x=0, color="black", linewidth=0.5)
    ax.set_title("C", loc="left", fontweight="bold")

    for i, (b, p) in enumerate(zip(betas, pvals)):
        label = f"p={p:.3f}" if p > 0.001 else f"p={p:.1e}"
        ax.text(b + (0.01 if b >= 0 else -0.01), i, label,
                va="center", fontsize=6, ha="left" if b >= 0 else "right")

    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(FIG_DIR / "figS3_convergence.pdf")
    plt.close(fig)
    logger.info("Figure S3 saved")


# =====================================================================
# Supplementary Figure S4: ROC curves + TCGA descriptive
# =====================================================================

def figure_s4(data: dict) -> None:
    """ROC curves for classification performance + TCGA-PAAD validation."""
    import json

    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.8))

    # Panel A: ROC curves
    ax = axes[0]
    roc_path = Path("results/h13h/roc_data.json")
    if roc_path.exists():
        with open(roc_path) as f:
            roc_data = json.load(f)
        drug_colors = {
            "Afuresertib": COLORS["target"],
            "Ipatasertib": COLORS["highlight"],
            "Capivasertib": COLORS["dark"],
        }
        for drug_label, roc in roc_data.items():
            ax.plot(
                roc["fpr"], roc["tpr"],
                color=drug_colors.get(drug_label, "gray"),
                lw=1.2,
                label=f"{drug_label} ({roc['auc']:.3f})",
            )
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, lw=0.8)
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_title("A", loc="left", fontweight="bold")
        ax.legend(loc="lower right", fontsize=6, framealpha=0.8)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
    else:
        ax.text(0.5, 0.5, "ROC data not found\nRun h13h first",
                ha="center", va="center", transform=ax.transAxes, fontsize=8)

    # Panel B: TCGA vs DepMap SLC2A8 distribution
    ax = axes[1]
    tcga_path = Path("results/h13h/tcga_descriptive.json")
    if tcga_path.exists():
        with open(tcga_path) as f:
            tcga_data = json.load(f)
        tcga_vals = np.array(tcga_data["slc2a8_values"])

        # DepMap PDAC
        expr = data["expression"]
        meta = data["metadata"]
        lineage_col = "OncotreeLineage"
        pdac_ids = meta.index[meta[lineage_col].str.contains("Pancrea", na=False)]
        depmap_vals = expr.loc[expr.index.isin(pdac_ids), "SLC2A8"].dropna().values

        ax.hist(tcga_vals, bins=20, alpha=0.6, color=COLORS["target"],
                label=f"TCGA (n={len(tcga_vals)})", density=True,
                edgecolor="black", linewidth=0.3)
        ax.hist(depmap_vals, bins=12, alpha=0.6, color=COLORS["background"],
                label=f"DepMap (n={len(depmap_vals)})", density=True,
                edgecolor="black", linewidth=0.3)
        ax.set_xlabel("SLC2A8 (log$_2$ TPM+1)")
        ax.set_ylabel("Density")
        ax.set_title("B", loc="left", fontweight="bold")
        ax.legend(fontsize=6, framealpha=0.8)
    else:
        ax.text(0.5, 0.5, "TCGA data not found\nRun h13h first",
                ha="center", va="center", transform=ax.transAxes, fontsize=8)

    # Panel C: TCGA co-expression replication
    ax = axes[2]
    tcga_coexpr = {
        "AKT1": {"rho": 0.212, "p": 0.005, "replicates": True},
        "FASN": {"rho": 0.180, "p": 0.016, "replicates": True},
        "MTOR": {"rho": -0.022, "p": 0.771, "replicates": False},
        "VPS35": {"rho": -0.171, "p": 0.023, "replicates": False},
        "SCD": {"rho": 0.031, "p": 0.684, "replicates": False},
    }
    genes = list(tcga_coexpr.keys())
    rhos = [tcga_coexpr[g]["rho"] for g in genes]
    repl = [tcga_coexpr[g]["replicates"] for g in genes]
    ps_t = [tcga_coexpr[g]["p"] for g in genes]
    colors_t = [COLORS["target"] if r else COLORS["null"] for r in repl]

    ax.barh(range(len(genes)), rhos, color=colors_t,
            edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlabel("Spearman $\\rho$ with SLC2A8 (TCGA)")
    ax.axvline(x=0, color="black", linewidth=0.5)
    ax.set_title("C", loc="left", fontweight="bold")

    for i, (r, p) in enumerate(zip(rhos, ps_t)):
        label = f"p={p:.3f}" if p > 0.001 else f"p={p:.1e}"
        ax.text(r + (0.02 if r >= 0 else -0.02), i, label,
                va="center", fontsize=6, ha="left" if r >= 0 else "right")

    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(FIG_DIR / "figS4_roc_tcga.pdf")
    plt.close(fig)
    logger.info("Figure S4 saved")


# =====================================================================
# Main
# =====================================================================

def main() -> None:
    """Generate all paper figures."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("Generating H13 paper figures to %s", FIG_DIR)

    data = _load_all()

    figure_1(data)
    figure_2(data)
    figure_3(data)
    figure_4(data)
    figure_5(data)
    figure_6()
    figure_s1(data)
    figure_s2(data)
    figure_s3(data)
    figure_s4(data)

    logger.info("All figures saved to %s", FIG_DIR)


if __name__ == "__main__":
    main()
