"""Composite scoring → final ranked candidate list.

final_score(gene) =
        w_effect * sigmoid(z_effect)
      + w_shap   * sigmoid(z_shap)
      + w_corr   * sigmoid(z_corr)
      multiplied by the per-gene baseline-normal penalty multiplier.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-z))


def _zscore(s: pd.Series) -> pd.Series:
    if s.std(ddof=0) == 0:
        return pd.Series(np.zeros(len(s)), index=s.index)
    return (s - s.mean()) / s.std(ddof=0)


def compose_ranking(
    pooled_c3: pd.DataFrame,
    shap_importance: pd.DataFrame,
    corr_with_risk: pd.DataFrame,
    baseline_penalty: pd.DataFrame,
    surface_kept: list[str],
    w_effect: float = 0.40,
    w_shap: float = 0.40,
    w_corr: float = 0.20,
    top_n: int = 25,
) -> pd.DataFrame:
    """Return final ranking restricted to surface-localised, surviving genes."""
    df = pd.DataFrame({"gene": surface_kept})
    df = df.merge(
        pooled_c3[["gene", "pooled_estimate", "p_fdr", "I2",
                   "n_cohorts_same_direction", "n_cohorts_total"]],
        on="gene", how="left",
    )
    df = df.merge(shap_importance, on="gene", how="left")
    df = df.merge(corr_with_risk, on="gene", how="left")
    df = df.merge(
        baseline_penalty[["gene", "baseline_penalty_multiplier",
                          "excluded_for_high_baseline"]],
        on="gene", how="left",
    )
    # Drop anything excluded outright by the baseline filter (defensive — should
    # already have been removed at the filter step).
    df = df[df["excluded_for_high_baseline"] != True].copy()

    df["abs_pooled_c3"] = df["pooled_estimate"].abs().fillna(0)
    df["mean_abs_shap"] = df["mean_abs_shap"].fillna(0)
    df["abs_corr"] = df["corr_with_risk"].abs().fillna(0)

    df["z_effect"] = _zscore(df["abs_pooled_c3"])
    df["z_shap"] = _zscore(df["mean_abs_shap"])
    df["z_corr"] = _zscore(df["abs_corr"])

    df["score_unpenalized"] = (
        w_effect * _sigmoid(df["z_effect"].values)
        + w_shap   * _sigmoid(df["z_shap"].values)
        + w_corr   * _sigmoid(df["z_corr"].values)
    )
    df["composite_score"] = (
        df["score_unpenalized"] * df["baseline_penalty_multiplier"].fillna(1.0)
    )
    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df.head(top_n)
