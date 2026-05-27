"""SHAP value computation for the L1 logistic fit (averaged across LOCO folds)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from opmd_msi.models.validation import LOCOResult


def _linear_shap_values(coefs: np.ndarray, X: np.ndarray, background: np.ndarray) -> np.ndarray:
    """SHAP for a linear model is (x - mean_background) * coef per feature.

    We avoid the SHAP runtime dependency for the demo path (and for CI without
    network access); the closed-form result is identical for linear models.
    """
    mean_bg = background.mean(axis=0)
    return (X - mean_bg) * coefs


def compute_mean_abs_shap(loco: LOCOResult) -> pd.DataFrame:
    """Average |SHAP| across LOCO folds for the L1 logistic models.

    Returns a DataFrame with columns: gene, mean_abs_shap.
    """
    if not loco.fold_fits_logistic:
        return pd.DataFrame(columns=["gene", "mean_abs_shap"])

    per_fold_abs_means: list[pd.Series] = []
    for cohort_id, fit in loco.fold_fits_logistic.items():
        if fit.train_X.empty:
            continue
        coefs = fit.coefs.values
        X = fit.train_X.values
        background = X  # use the training set as the background distribution
        shap = _linear_shap_values(coefs, X, background)
        mean_abs = pd.Series(np.abs(shap).mean(axis=0), index=fit.feature_columns)
        mean_abs.name = cohort_id
        per_fold_abs_means.append(mean_abs)

    if not per_fold_abs_means:
        return pd.DataFrame(columns=["gene", "mean_abs_shap"])
    matrix = pd.concat(per_fold_abs_means, axis=1)
    out = matrix.mean(axis=1).rename("mean_abs_shap").to_frame()
    out.index.name = "gene"
    return out.reset_index()


def compute_correlation_with_risk(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    risk_col: str = "risk_class",
    case_value: str = "high_risk_OPMD",
    control_value: str = "stable_OPMD",
) -> pd.DataFrame:
    """Pearson correlation between each gene's expression and binary risk label."""
    sub_meta = metadata[metadata[risk_col].isin([case_value, control_value])]
    if sub_meta.empty:
        return pd.DataFrame(columns=["gene", "corr_with_risk"])
    y = (sub_meta[risk_col] == case_value).astype(int).values
    X = expression.loc[sub_meta.index, genes].astype(float).values
    means_y = y.mean()
    std_y = y.std(ddof=0)
    means_x = X.mean(axis=0)
    std_x = X.std(axis=0, ddof=0)
    cov = ((X - means_x) * (y - means_y)[:, None]).mean(axis=0)
    denom = np.where(std_x * std_y == 0, 1.0, std_x * std_y)
    corr = cov / denom
    return pd.DataFrame({"gene": genes, "corr_with_risk": corr})
