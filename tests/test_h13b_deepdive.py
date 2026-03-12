"""Unit tests for H13b autophagy initiation deep dive.

Tests cover refined gene set definitions, lysosomal correlation, focused
enrichment (initiation vs execution), tertile split, PRISM drug sensitivity,
lineage comparison, patient scoring, survival analysis, subtype correlation,
and a mini integration test.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pdac.h13_trehalose_vulnerability.enrichment import check_dependency_enrichment
from pdac.h13_trehalose_vulnerability.h13b_gene_sets import (
    AUTOPHAGY_DRUGS_TIER1,
    AUTOPHAGY_DRUGS_TIER2,
    AUTOPHAGY_EXECUTION,
    AUTOPHAGY_INITIATION,
    COMPARISON_LINEAGES,
    LYSOSOMAL_PROGRAM,
)
from pdac.h13_trehalose_vulnerability.h13b_lineage import (
    compute_slc2a8_lysosomal_correlation,
    get_lineage_model_ids,
    run_lineage_enrichment_comparison,
    run_tertile_enrichment,
)
from pdac.h13_trehalose_vulnerability.h13b_patient import (
    correlate_with_subtype,
    kaplan_meier_split,
    run_bootstrap_cox,
    run_survival_analysis,
    score_tcga_patients,
)
from pdac.h13_trehalose_vulnerability.h13b_prism import (
    identify_autophagy_compounds,
    compute_drug_sensitivity,
)
from pdac.h13_trehalose_vulnerability.scoring import compute_gene_set_zscore


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def rng() -> np.random.Generator:
    return np.random.default_rng(42)


@pytest.fixture()
def synthetic_expression(rng: np.random.Generator) -> pd.DataFrame:
    """100 cell lines with all initiation + execution + lysosomal genes.

    First 30 lines are SLC2A8-high, rest are SLC2A8-low.
    """
    n_lines = 100
    all_genes = (
        ["SLC2A8", "SLC2A1"]
        + AUTOPHAGY_INITIATION
        + AUTOPHAGY_EXECUTION
        + LYSOSOMAL_PROGRAM
        + [f"RANDOM_{i}" for i in range(20)]
    )
    # Deduplicate (LAMP1, LAMP2, TFEB appear in both execution and lysosomal)
    all_genes = list(dict.fromkeys(all_genes))

    ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    data = rng.normal(5.0, 2.0, size=(n_lines, len(all_genes)))

    # Make SLC2A8 bimodal
    data[:30, 0] += 3.0
    data[30:, 0] -= 1.0

    # Plant lysosomal correlation: SLC2A8-high lines have higher LAMP1/LAMP2
    lamp1_idx = all_genes.index("LAMP1")
    lamp2_idx = all_genes.index("LAMP2")
    data[:30, lamp1_idx] += 1.5
    data[:30, lamp2_idx] += 1.5

    return pd.DataFrame(data, index=ids, columns=all_genes)


@pytest.fixture()
def synthetic_crispr(rng: np.random.Generator) -> pd.DataFrame:
    """100 cell lines with planted initiation dependency signal.

    First 30 lines (SLC2A8-high) are MORE dependent on INITIATION genes
    but NOT on EXECUTION genes. This is the key H13b hypothesis.
    """
    n_lines = 100
    genes = list(dict.fromkeys(AUTOPHAGY_INITIATION + AUTOPHAGY_EXECUTION))
    ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    data = rng.normal(-0.2, 0.5, size=(n_lines, len(genes)))

    # Plant initiation dependency in first 30 lines
    for i, gene in enumerate(genes):
        if gene in AUTOPHAGY_INITIATION:
            data[:30, i] -= 0.5  # stronger negative = more dependent

    return pd.DataFrame(data, index=ids, columns=genes)


@pytest.fixture()
def ras_mutant_ids() -> list[str]:
    """First 60 cell lines are RAS-mutant."""
    return [f"ACH-{i:06d}" for i in range(60)]


@pytest.fixture()
def synthetic_prism(rng: np.random.Generator) -> pd.DataFrame:
    """100 cell lines x 10 compounds with planted sensitivity."""
    n_lines = 100
    compounds = [f"BRD-{i:04d}" for i in range(10)]
    ids = [f"ACH-{i:06d}" for i in range(n_lines)]
    data = rng.normal(-0.5, 1.0, size=(n_lines, len(compounds)))

    # Plant sensitivity for first 2 compounds (Tier 1) in SLC2A8-high lines
    data[:30, 0] -= 1.0
    data[:30, 1] -= 0.8

    return pd.DataFrame(data, index=ids, columns=compounds)


@pytest.fixture()
def synthetic_clinical(rng: np.random.Generator) -> pd.DataFrame:
    """80 patients with survival data."""
    n = 80
    df = pd.DataFrame({
        "patient_id": [f"TCGA-{i:04d}" for i in range(n)],
        "os_days": rng.exponential(500, n) + 30,
        "os_event": rng.binomial(1, 0.6, n).astype(float),
        "age": rng.normal(65, 10, n),
        "stage_numeric": rng.choice([1.0, 2.0, 3.0, 4.0], n),
    })
    return df


@pytest.fixture()
def synthetic_tcga_expression(rng: np.random.Generator) -> pd.DataFrame:
    """80 patients x genes covering all scoring sets."""
    from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection

    gsc = GeneSetCollection()
    all_genes = gsc.all_expression_genes()

    # Add Moffitt subtype genes for classify_moffitt_subtype
    from pdac.h03_latent_divergence.subtypes import BASAL_GENES, CLASSICAL_GENES
    all_genes = list(dict.fromkeys(all_genes + BASAL_GENES + CLASSICAL_GENES))

    n = 80
    ids = [f"TCGA-{i:04d}" for i in range(n)]
    data = rng.normal(5.0, 2.0, size=(n, len(all_genes)))

    # Plant vulnerability signal: first 40 patients have higher SLC2A8
    slc2a8_idx = all_genes.index("SLC2A8")
    data[:40, slc2a8_idx] += 2.0
    data[40:, slc2a8_idx] -= 1.0

    # Plant subtype signal: first 20 are basal-like
    for gene in BASAL_GENES:
        if gene in all_genes:
            idx = all_genes.index(gene)
            data[:20, idx] += 2.0
    for gene in CLASSICAL_GENES:
        if gene in all_genes:
            idx = all_genes.index(gene)
            data[20:, idx] += 2.0

    return pd.DataFrame(data, index=ids, columns=all_genes)


@pytest.fixture()
def synthetic_metadata() -> pd.DataFrame:
    """Metadata for 100 cell lines across 3 lineages."""
    ids = [f"ACH-{i:06d}" for i in range(100)]
    lineages = (
        ["Pancreas"] * 30
        + ["Bowel"] * 35
        + ["Esophagus/Stomach"] * 35
    )
    return pd.DataFrame(
        {"OncotreeLineage": lineages},
        index=ids,
    )


@pytest.fixture()
def synthetic_mutations() -> pd.DataFrame:
    """Mutations for first 60 cell lines in KRAS."""
    rows = []
    for i in range(60):
        rows.append({
            "ModelID": f"ACH-{i:06d}",
            "HugoSymbol": "KRAS",
            "VariantInfo": "missense",
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Test classes
# ---------------------------------------------------------------------------


class TestRefinedGeneSets:
    """Validate gene set definitions."""

    def test_no_overlap_initiation_execution(self) -> None:
        init_set = set(AUTOPHAGY_INITIATION)
        exec_set = set(AUTOPHAGY_EXECUTION)
        assert init_set.isdisjoint(exec_set), f"Overlap: {init_set & exec_set}"

    def test_initiation_size(self) -> None:
        assert len(AUTOPHAGY_INITIATION) == 6

    def test_execution_size(self) -> None:
        assert len(AUTOPHAGY_EXECUTION) == 8

    def test_no_duplicates_within_sets(self) -> None:
        assert len(AUTOPHAGY_INITIATION) == len(set(AUTOPHAGY_INITIATION))
        assert len(AUTOPHAGY_EXECUTION) == len(set(AUTOPHAGY_EXECUTION))
        assert len(LYSOSOMAL_PROGRAM) == len(set(LYSOSOMAL_PROGRAM))

    def test_lysosomal_program_size(self) -> None:
        assert len(LYSOSOMAL_PROGRAM) == 6

    def test_drug_tiers_defined(self) -> None:
        assert len(AUTOPHAGY_DRUGS_TIER1) >= 1
        assert len(AUTOPHAGY_DRUGS_TIER2) >= 1

    def test_lineages_defined(self) -> None:
        assert "pancreas" in COMPARISON_LINEAGES
        assert len(COMPARISON_LINEAGES) >= 3

    def test_all_genes_are_strings(self) -> None:
        for gene in AUTOPHAGY_INITIATION + AUTOPHAGY_EXECUTION + LYSOSOMAL_PROGRAM:
            assert isinstance(gene, str)
            assert len(gene) > 0


class TestLysosomalCorrelation:
    """Spearman correlation between SLC2A8 and lysosomal genes."""

    def test_computes_correlation(
        self, synthetic_expression: pd.DataFrame,
    ) -> None:
        result = compute_slc2a8_lysosomal_correlation(
            synthetic_expression, "SLC2A8", LYSOSOMAL_PROGRAM,
        )
        assert len(result) > 0
        assert "rho" in result.columns
        assert "p_value" in result.columns
        assert "fdr_rejected" in result.columns

    def test_detects_planted_correlation(
        self, synthetic_expression: pd.DataFrame,
    ) -> None:
        """LAMP1/LAMP2 were planted with positive correlation to SLC2A8."""
        result = compute_slc2a8_lysosomal_correlation(
            synthetic_expression, "SLC2A8", ["LAMP1", "LAMP2"],
        )
        # Planted correlation should be positive
        lamp1_rho = result.loc[result["gene"] == "LAMP1", "rho"].values[0]
        assert lamp1_rho > 0

    def test_missing_target_gene(self) -> None:
        expr = pd.DataFrame({"A": [1, 2, 3]}, index=["s1", "s2", "s3"])
        result = compute_slc2a8_lysosomal_correlation(expr, "MISSING", ["A"])
        assert result.empty

    def test_no_lysosomal_genes_found(
        self, synthetic_expression: pd.DataFrame,
    ) -> None:
        result = compute_slc2a8_lysosomal_correlation(
            synthetic_expression, "SLC2A8", ["NOT_A_GENE"],
        )
        assert result.empty


class TestFocusedEnrichment:
    """6-gene initiation set detects signal; 8-gene execution set does NOT."""

    def test_initiation_detects_planted_signal(
        self, synthetic_crispr: pd.DataFrame,
    ) -> None:
        target_ids = [f"ACH-{i:06d}" for i in range(30)]
        background_ids = [f"ACH-{i:06d}" for i in range(30, 100)]

        result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_INITIATION, "autophagy_initiation",
        )
        assert result.set_level_statistic < 0  # target more dependent
        assert result.set_level_p < 0.05

    def test_execution_does_not_detect_signal(
        self, synthetic_crispr: pd.DataFrame,
    ) -> None:
        target_ids = [f"ACH-{i:06d}" for i in range(30)]
        background_ids = [f"ACH-{i:06d}" for i in range(30, 100)]

        result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_EXECUTION, "autophagy_execution",
        )
        # Execution should NOT be significant (no planted signal)
        assert result.set_level_p > 0.05


class TestMeanDepShift:
    """Mean dependency shift: initiation should be more negative."""

    def test_initiation_shift_more_negative(
        self, synthetic_crispr: pd.DataFrame,
    ) -> None:
        target_ids = [f"ACH-{i:06d}" for i in range(30)]
        background_ids = [f"ACH-{i:06d}" for i in range(30, 100)]

        init_result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_INITIATION, "initiation",
        )
        exec_result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_EXECUTION, "execution",
        )

        assert init_result.set_level_statistic < exec_result.set_level_statistic


class TestTertileEnrichment:
    """Tertile split produces signal on planted data."""

    def test_tertile_runs_correctly(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        ras_mutant_ids: list[str],
    ) -> None:
        result = run_tertile_enrichment(
            synthetic_crispr, synthetic_expression, ras_mutant_ids,
            "SLC2A8", AUTOPHAGY_INITIATION, "initiation",
        )
        assert result.n_target > 0
        assert result.n_background > 0
        assert result.gene_set_name == "initiation_tertile"

    def test_tertile_detects_planted_signal(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        ras_mutant_ids: list[str],
    ) -> None:
        result = run_tertile_enrichment(
            synthetic_crispr, synthetic_expression, ras_mutant_ids,
            "SLC2A8", AUTOPHAGY_INITIATION, "initiation",
        )
        # Signal should be present (planted in first 30 which are SLC2A8-high)
        assert result.set_level_statistic < 0

    def test_tertile_missing_gene(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        ras_mutant_ids: list[str],
    ) -> None:
        result = run_tertile_enrichment(
            synthetic_crispr, synthetic_expression, ras_mutant_ids,
            "NOT_A_GENE", AUTOPHAGY_INITIATION, "initiation",
        )
        assert result.set_level_p == 1.0


class TestPrismSensitivity:
    """Drug sensitivity tests with synthetic PRISM data."""

    def test_detects_planted_sensitivity(
        self, synthetic_prism: pd.DataFrame,
    ) -> None:
        target_ids = [f"ACH-{i:06d}" for i in range(30)]
        background_ids = [f"ACH-{i:06d}" for i in range(30, 100)]
        compounds = synthetic_prism.columns.tolist()

        result = compute_drug_sensitivity(
            synthetic_prism, target_ids, background_ids, compounds,
        )
        assert len(result) > 0
        # First compound should be most significant (planted strongest)
        top = result.iloc[0]
        assert top["mean_diff"] < 0  # target more sensitive

    def test_insufficient_cell_lines(self) -> None:
        sensitivity = pd.DataFrame(
            {"BRD-0001": [1.0]}, index=["ACH-000000"],
        )
        result = compute_drug_sensitivity(
            sensitivity, ["ACH-000000"], [], ["BRD-0001"],
        )
        assert result.empty

    def test_identify_compounds(self) -> None:
        compound_info = pd.DataFrame({
            "name": ["SBI-0206965", "SAR405", "aspirin", "rapamycin"],
            "moa": ["ULK1 inhibitor", "VPS34 inhibitor", "COX inhibitor", "mTOR inhibitor"],
            "broad_id": ["BRD-1", "BRD-2", "BRD-3", "BRD-4"],
        })
        result = identify_autophagy_compounds(
            compound_info,
            {
                "ULK1_inhibitor": ["SBI-0206965"],
                "VPS34_inhibitor": ["SAR405"],
            },
        )
        assert len(result) >= 2
        assert "tier" in result.columns

    def test_tier_classification(self) -> None:
        compound_info = pd.DataFrame({
            "name": ["SBI-0206965", "chloroquine"],
            "moa": ["ULK1 inhibitor", "lysosomal inhibitor"],
            "broad_id": ["BRD-1", "BRD-2"],
        })
        result = identify_autophagy_compounds(
            compound_info,
            {
                "ULK1_inhibitor": ["SBI-0206965"],
                "lysosomal_inhibitor": ["chloroquine"],
            },
            moa_keywords=[],
        )
        assert len(result) == 2
        tiers = set(result["tier"])
        assert "ULK1_inhibitor" in tiers
        assert "lysosomal_inhibitor" in tiers


class TestLineageComparison:
    """Per-lineage enrichment comparison."""

    def test_returns_results_per_lineage(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        synthetic_metadata: pd.DataFrame,
        synthetic_mutations: pd.DataFrame,
    ) -> None:
        results = run_lineage_enrichment_comparison(
            synthetic_crispr, synthetic_expression,
            synthetic_mutations, synthetic_metadata,
            AUTOPHAGY_INITIATION, "initiation",
            COMPARISON_LINEAGES,
            ["KRAS"],
        )
        assert "pancreas" in results
        assert "bowel" in results
        assert "esophagus_stomach" in results

    def test_get_lineage_model_ids(
        self, synthetic_metadata: pd.DataFrame,
    ) -> None:
        ids = get_lineage_model_ids(synthetic_metadata, "Pancrea")
        assert len(ids) == 30

    def test_lineage_result_is_enrichment_result(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        synthetic_metadata: pd.DataFrame,
        synthetic_mutations: pd.DataFrame,
    ) -> None:
        from pdac.h13_trehalose_vulnerability.enrichment import EnrichmentResult

        results = run_lineage_enrichment_comparison(
            synthetic_crispr, synthetic_expression,
            synthetic_mutations, synthetic_metadata,
            AUTOPHAGY_INITIATION, "initiation",
            {"pancreas": "Pancrea"},
            ["KRAS"],
        )
        assert isinstance(results["pancreas"], EnrichmentResult)


class TestPatientScoring:
    """Vulnerability score computation on TCGA-like data."""

    def test_scores_computed(
        self, synthetic_tcga_expression: pd.DataFrame,
    ) -> None:
        from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection

        gsc = GeneSetCollection()
        scores = score_tcga_patients(
            synthetic_tcga_expression, gsc.scoring_sets(), gsc.target,
        )
        assert len(scores) == len(synthetic_tcga_expression)
        assert scores.name == "vulnerability_score"

    def test_high_slc2a8_higher_scores(
        self, synthetic_tcga_expression: pd.DataFrame,
    ) -> None:
        """First 40 patients have higher SLC2A8 -> should have higher scores."""
        from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection

        gsc = GeneSetCollection()
        scores = score_tcga_patients(
            synthetic_tcga_expression, gsc.scoring_sets(), gsc.target,
        )
        high_mean = scores.iloc[:40].mean()
        low_mean = scores.iloc[40:].mean()
        assert high_mean > low_mean


class TestSurvivalAnalysis:
    """Cox fit, KM split, bootstrap on synthetic clinical data."""

    def test_kaplan_meier_split(
        self, synthetic_clinical: pd.DataFrame,
    ) -> None:
        df = synthetic_clinical.set_index("patient_id")
        df["vulnerability_score"] = np.random.default_rng(42).normal(0, 1, len(df))

        result = kaplan_meier_split(df, "vulnerability_score")
        assert "log_rank_p" in result
        assert result["n_high"] > 0
        assert result["n_low"] > 0
        assert 0 <= result["log_rank_p"] <= 1

    def test_run_survival_analysis(
        self, synthetic_clinical: pd.DataFrame,
    ) -> None:
        df = synthetic_clinical.set_index("patient_id")
        df["vulnerability_score"] = np.random.default_rng(42).normal(0, 1, len(df))

        result = run_survival_analysis(df)
        assert "cox_result" in result
        assert "km_result" in result
        assert "c_index_improvement" in result

    def test_bootstrap_cox_runs(
        self, synthetic_clinical: pd.DataFrame,
    ) -> None:
        df = synthetic_clinical.set_index("patient_id")
        df["vulnerability_score"] = np.random.default_rng(42).normal(0, 1, len(df))

        result = run_bootstrap_cox(
            df, ["vulnerability_score", "age", "stage_numeric"],
            n_bootstrap=50,
        )
        assert "hr_median" in result
        assert "hr_ci_lower" in result
        assert "hr_ci_upper" in result
        assert "ci_excludes_1" in result

    def test_survival_missing_columns_raises(self) -> None:
        df = pd.DataFrame({"os_days": [100], "os_event": [1]})
        with pytest.raises(ValueError, match="Missing required columns"):
            run_survival_analysis(df)


class TestSubtypeCorrelation:
    """Moffitt classifier integration with vulnerability scores."""

    def test_correlation_computes(
        self,
        synthetic_tcga_expression: pd.DataFrame,
    ) -> None:
        from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection

        gsc = GeneSetCollection()
        scores = score_tcga_patients(
            synthetic_tcga_expression, gsc.scoring_sets(), gsc.target,
        )

        result = correlate_with_subtype(synthetic_tcga_expression, scores)
        assert "subtype_labels" in result
        assert "mann_whitney_p" in result
        assert "point_biserial_r" in result

    def test_subtype_labels_are_basal_classical(
        self,
        synthetic_tcga_expression: pd.DataFrame,
    ) -> None:
        from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection

        gsc = GeneSetCollection()
        scores = score_tcga_patients(
            synthetic_tcga_expression, gsc.scoring_sets(), gsc.target,
        )

        result = correlate_with_subtype(synthetic_tcga_expression, scores)
        labels = result["subtype_labels"]
        if not labels.empty:
            assert set(labels.unique()) <= {"basal", "classical"}

    def test_too_few_patients(self) -> None:
        expr = pd.DataFrame({"A": [1]}, index=["p1"])
        scores = pd.Series([1.0], index=["p1"], name="vulnerability_score")
        result = correlate_with_subtype(expr, scores)
        assert np.isnan(result["mann_whitney_p"])


class TestMiniPipeline:
    """Integration test for DepMap arm only."""

    def test_focused_enrichment_flow(
        self,
        synthetic_crispr: pd.DataFrame,
        synthetic_expression: pd.DataFrame,
        ras_mutant_ids: list[str],
    ) -> None:
        """End-to-end: scoring -> stratification -> focused enrichment."""
        from pdac.h13_trehalose_vulnerability.gene_sets import GeneSetCollection
        from pdac.h13_trehalose_vulnerability.scoring import (
            compute_vulnerability_score,
            stratify_by_target_and_ras,
        )

        gsc = GeneSetCollection()

        # Score
        scores = compute_vulnerability_score(
            synthetic_expression, gsc.target, gsc.scoring_sets(),
        )
        assert len(scores) == 100

        # Stratify
        labels = stratify_by_target_and_ras(scores, ras_mutant_ids)
        target_ids = labels[labels == "SLC2A8-high / RAS-mut"].index.tolist()
        background_ids = [i for i in scores.index if i not in target_ids]
        assert len(target_ids) > 0

        # Focused enrichment
        init_result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_INITIATION, "initiation",
        )
        exec_result = check_dependency_enrichment(
            synthetic_crispr, target_ids, background_ids,
            AUTOPHAGY_EXECUTION, "execution",
        )

        # Initiation should have stronger signal
        assert init_result.set_level_statistic < exec_result.set_level_statistic

        # Lysosomal correlation
        corr = compute_slc2a8_lysosomal_correlation(
            synthetic_expression, "SLC2A8", LYSOSOMAL_PROGRAM,
        )
        assert len(corr) > 0

    def test_falsification_checks(self) -> None:
        """Test falsification check functions with mock data."""
        from pdac.h13_trehalose_vulnerability.enrichment import EnrichmentResult
        from pdac.h13_trehalose_vulnerability.h13b_pipeline import (
            _check_f3b_focused_enrichment,
            _check_f5_prism,
            _check_f6_survival,
            _check_f7_lineage,
        )

        # F3b: pass when initiation p < 0.05 and execution p > 0.1
        init = EnrichmentResult(
            gene_set_name="init", n_target=20, n_background=80,
            genes_tested=AUTOPHAGY_INITIATION,
            gene_level_results=pd.DataFrame(),
            set_level_statistic=-0.3, set_level_p=0.01,
            fdr_significant_genes=["ATG13"],
        )
        exec_ctrl = EnrichmentResult(
            gene_set_name="exec", n_target=20, n_background=80,
            genes_tested=AUTOPHAGY_EXECUTION,
            gene_level_results=pd.DataFrame(),
            set_level_statistic=-0.05, set_level_p=0.5,
            fdr_significant_genes=[],
        )
        passed, _ = _check_f3b_focused_enrichment(init, exec_ctrl)
        assert passed

        # F3b: fail when execution is also significant
        exec_sig = EnrichmentResult(
            gene_set_name="exec", n_target=20, n_background=80,
            genes_tested=AUTOPHAGY_EXECUTION,
            gene_level_results=pd.DataFrame(),
            set_level_statistic=-0.2, set_level_p=0.03,
            fdr_significant_genes=[],
        )
        passed, _ = _check_f3b_focused_enrichment(init, exec_sig)
        assert not passed

        # F5: pass with tier 1 compound
        prism = pd.DataFrame({
            "compound": ["BRD-1", "BRD-2"],
            "p_value": [0.05, 0.3],
            "mean_diff": [-0.5, -0.1],
        })
        passed, _ = _check_f5_prism(prism, ["BRD-1"])
        assert passed

        # F5: fail with no sig
        prism_fail = pd.DataFrame({
            "compound": ["BRD-1"],
            "p_value": [0.5],
            "mean_diff": [-0.1],
        })
        passed, _ = _check_f5_prism(prism_fail, ["BRD-1"])
        assert not passed

        # F6: pass when CI excludes 1
        bootstrap = {
            "hr_median": 1.5, "hr_ci_lower": 1.1,
            "hr_ci_upper": 2.0, "ci_excludes_1": True,
        }
        passed, _ = _check_f6_survival(bootstrap)
        assert passed

        # F6: fail when CI includes 1
        bootstrap_fail = {
            "hr_median": 1.1, "hr_ci_lower": 0.8,
            "hr_ci_upper": 1.5, "ci_excludes_1": False,
        }
        passed, _ = _check_f6_survival(bootstrap_fail)
        assert not passed

        # F7: pass when PDAC is strongest
        pdac_result = EnrichmentResult(
            gene_set_name="pdac", n_target=10, n_background=20,
            genes_tested=AUTOPHAGY_INITIATION,
            gene_level_results=pd.DataFrame(),
            set_level_statistic=-0.3, set_level_p=0.01,
            fdr_significant_genes=[],
        )
        bowel_result = EnrichmentResult(
            gene_set_name="bowel", n_target=10, n_background=20,
            genes_tested=AUTOPHAGY_INITIATION,
            gene_level_results=pd.DataFrame(),
            set_level_statistic=-0.1, set_level_p=0.3,
            fdr_significant_genes=[],
        )
        passed, _ = _check_f7_lineage({"pancreas": pdac_result, "bowel": bowel_result})
        assert passed

        # F7: fail when bowel is stronger
        passed, _ = _check_f7_lineage({"pancreas": bowel_result, "bowel": pdac_result})
        assert not passed

    def test_h13b_results_verdict(self) -> None:
        """Test H13bResults verdict logic."""
        from pdac.h13_trehalose_vulnerability.h13b_pipeline import H13bResults

        results = H13bResults()
        results.falsification = {"F3b": True, "F5": True, "F6": True, "F7": True}
        assert results.verdict() == "PASS"

        results.falsification = {"F3b": True, "F5": True, "F6": False, "F7": False}
        assert results.verdict() == "PARTIAL"

        results.falsification = {"F3b": True, "F5": False, "F6": False, "F7": False}
        assert results.verdict() == "FAIL"
