"""Per-gene, per-cohort effect-size estimators."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class EffectEstimate:
    estimate: float    # Hedges' g, paired log2FC, or Cliff's delta
    variance: float    # variance of the estimate (for inverse-variance meta-analysis)
    n_case: int
    n_control: int
    p_value: float


def hedges_g(case: np.ndarray, control: np.ndarray) -> EffectEstimate:
    n1, n2 = len(case), len(control)
    if n1 < 2 or n2 < 2:
        return EffectEstimate(0.0, np.inf, n1, n2, 1.0)
    m1, m2 = float(np.mean(case)), float(np.mean(control))
    v1, v2 = float(np.var(case, ddof=1)), float(np.var(control, ddof=1))
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    if sp <= 0:
        return EffectEstimate(0.0, np.inf, n1, n2, 1.0)
    d = (m1 - m2) / sp
    # Small-sample correction
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * d
    # Variance of Hedges' g (Hedges & Olkin 1985)
    var_g = (n1 + n2) / (n1 * n2) + g**2 / (2 * (n1 + n2))
    # Two-sided Welch t-test for p
    t_stat, p = stats.ttest_ind(case, control, equal_var=False)
    return EffectEstimate(float(g), float(var_g), n1, n2, float(p))


def paired_log2_fold_change(case: np.ndarray, control: np.ndarray) -> EffectEstimate:
    """Within-subject paired log2 fold change. `case` and `control` are aligned."""
    n = len(case)
    if n != len(control) or n < 2:
        return EffectEstimate(0.0, np.inf, n, n, 1.0)
    diffs = case - control                 # already on log2 scale (post-harmonize)
    mean_d = float(np.mean(diffs))
    var_d = float(np.var(diffs, ddof=1) / n)
    t_stat, p = stats.ttest_rel(case, control)
    return EffectEstimate(mean_d, var_d, n, n, float(p))


def cliffs_delta(case: np.ndarray, control: np.ndarray) -> EffectEstimate:
    """Non-parametric effect size; robust to non-normal expression distributions."""
    n1, n2 = len(case), len(control)
    if n1 == 0 or n2 == 0:
        return EffectEstimate(0.0, np.inf, n1, n2, 1.0)
    # delta = (sum_{i,j} sign(x_i - y_j)) / (n1 * n2)
    diffs = case[:, None] - control[None, :]
    delta = float(np.sign(diffs).sum() / (n1 * n2))
    # Approximate variance via Cliff (1996); avoids zero for ties
    var = float((1 - delta**2) / (n1 + n2 - 1))
    # Mann-Whitney for p
    _, p = stats.mannwhitneyu(case, control, alternative="two-sided")
    return EffectEstimate(delta, var, n1, n2, float(p))


def compute_per_gene_effects(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    case_filter: str,
    control_filter: str,
    method: str = "hedges_g",
    paired_by: str | None = None,
) -> pd.DataFrame:
    """Compute per-gene effect sizes for one (cohort, contrast).

    `case_filter` and `control_filter` are pandas `.query()` expressions.
    Expression rows are samples; columns are genes. Metadata index must align
    with expression index.
    """
    case_ids = metadata.query(case_filter).index
    control_ids = metadata.query(control_filter).index
    if len(case_ids) == 0 or len(control_ids) == 0:
        return pd.DataFrame(columns=["gene", "estimate", "variance",
                                     "n_case", "n_control", "p_value"])

    estimator = {
        "hedges_g": hedges_g,
        "cliffs_delta": cliffs_delta,
        "paired_log2_fold_change": paired_log2_fold_change,
    }[method]

    rows = []
    if paired_by is not None and method == "paired_log2_fold_change":
        case_meta = metadata.loc[case_ids]
        ctrl_meta = metadata.loc[control_ids]
        common = sorted(set(case_meta[paired_by]) & set(ctrl_meta[paired_by]))
        if not common:
            return pd.DataFrame(columns=["gene", "estimate", "variance",
                                         "n_case", "n_control", "p_value"])
        case_aligned = case_meta.reset_index().set_index(paired_by).loc[common, "sample_id"]
        ctrl_aligned = ctrl_meta.reset_index().set_index(paired_by).loc[common, "sample_id"]
        case_expr = expression.loc[case_aligned.values]
        ctrl_expr = expression.loc[ctrl_aligned.values]
        for gene in expression.columns:
            est = estimator(case_expr[gene].values, ctrl_expr[gene].values)
            rows.append({"gene": gene, "estimate": est.estimate, "variance": est.variance,
                         "n_case": est.n_case, "n_control": est.n_control,
                         "p_value": est.p_value})
    else:
        case_expr = expression.loc[case_ids]
        ctrl_expr = expression.loc[control_ids]
        for gene in expression.columns:
            est = estimator(case_expr[gene].values, ctrl_expr[gene].values)
            rows.append({"gene": gene, "estimate": est.estimate, "variance": est.variance,
                         "n_case": est.n_case, "n_control": est.n_control,
                         "p_value": est.p_value})
    return pd.DataFrame(rows)
