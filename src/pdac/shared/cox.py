"""Shared Cox proportional-hazards model fitting for all hypothesis modules.

Provides a generic fit_cox function and CoxResult dataclass reusable across
H01, H04, H05, etc. without hypothesis-specific coupling.
"""

from dataclasses import dataclass

import pandas as pd
from lifelines import CoxPHFitter

from pdac.shared.logging import get_logger

logger = get_logger(__name__)


@dataclass
class CoxResult:
    """Summarized output from a Cox PH fit.

    Attributes:
        concordance_index: Harrell's C-index on the training data.
        hazard_ratios: Series of HR per covariate.
        p_values: Series of p-values per covariate.
        log_likelihood: Log partial likelihood at convergence.
    """

    concordance_index: float
    hazard_ratios: pd.Series
    p_values: pd.Series
    log_likelihood: float


def fit_cox(
    df: pd.DataFrame,
    covariates: list[str],
    duration_col: str = "os_days",
    event_col: str = "os_event",
    penalizer: float = 0.01,
) -> CoxResult:
    """Fit a multivariate Cox PH model with L2 ridge penalization.

    Uses a small L2 penalty to stabilize convergence when covariates have
    near-complete separation (e.g., KRAS mutated in ~95% of PDAC patients).

    Args:
        df: DataFrame containing survival data and covariates.
        covariates: List of covariate column names to include in the model.
        duration_col: Column name for survival duration (days).
        event_col: Column name for event indicator (1 = event, 0 = censored).
        penalizer: L2 ridge penalty strength. Small values (0.01) stabilize
            without materially affecting HR estimates.

    Returns:
        CoxResult with C-index, hazard ratios, and p-values.

    Raises:
        ValueError: If required columns are missing from df.
    """
    required = {duration_col, event_col} | set(covariates)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    logger.info(
        "Fitting Cox PH model",
        extra={"n_patients": len(df), "covariates": covariates},
    )

    fit_df = df[[duration_col, event_col, *covariates]].dropna()
    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=0.0)
    cph.fit(fit_df, duration_col=duration_col, event_col=event_col)

    result = CoxResult(
        concordance_index=cph.concordance_index_,
        hazard_ratios=cph.hazard_ratios_,
        p_values=cph.summary["p"],
        log_likelihood=cph.log_likelihood_,
    )
    logger.info(
        "Cox fit complete",
        extra={"c_index": round(result.concordance_index, 4)},
    )
    return result
