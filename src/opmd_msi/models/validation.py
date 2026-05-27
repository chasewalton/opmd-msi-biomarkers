"""Leave-one-cohort-out cross-validation harness.

For each cohort C in the dataset:
  - Train on all OTHER cohorts.
  - Score the model on cohort C.
  - Record the fitted coefficient vector and the held-out concordance/AUC.

Cohort-averaged coefficient magnitudes (and SHAP magnitudes computed in
`importance/shap_ranking.py`) form the basis of the final ranking.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from opmd_msi.models.cox import CoxFit, fit_l1_cox
from opmd_msi.models.logistic import LogisticFit, fit_l1_logistic


@dataclass
class LOCOResult:
    fold_fits_cox: dict[str, CoxFit] = field(default_factory=dict)
    fold_fits_logistic: dict[str, LogisticFit] = field(default_factory=dict)
    held_out_scores: dict[str, dict[str, float]] = field(default_factory=dict)
    coef_matrix_cox: pd.DataFrame = field(default_factory=pd.DataFrame)
    coef_matrix_logistic: pd.DataFrame = field(default_factory=pd.DataFrame)

    def average_abs_coefs(self) -> pd.DataFrame:
        """|coef| averaged across folds, per gene and per model."""
        cox = self.coef_matrix_cox.abs().mean(axis=0) if not self.coef_matrix_cox.empty else pd.Series(dtype=float)
        log = self.coef_matrix_logistic.abs().mean(axis=0) if not self.coef_matrix_logistic.empty else pd.Series(dtype=float)
        out = pd.DataFrame({"mean_abs_coef_cox": cox, "mean_abs_coef_logistic": log})
        out.index.name = "gene"
        return out.reset_index()


def leave_one_cohort_out(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    pooled_c1: pd.DataFrame | None = None,
    cohort_col: str = "cohort_id",
    cohorts_for_cox: list[str] | None = None,
    cohorts_for_logistic: list[str] | None = None,
) -> LOCOResult:
    """Train per-fold L1 Cox + L1 logistic, holding out one cohort at a time."""
    res = LOCOResult()

    all_cohorts = sorted(metadata[cohort_col].unique())
    cox_cohorts = [c for c in all_cohorts
                   if cohorts_for_cox is None or c in cohorts_for_cox]
    log_cohorts = [c for c in all_cohorts
                   if cohorts_for_logistic is None or c in cohorts_for_logistic]

    cox_rows: list[pd.Series] = []
    log_rows: list[pd.Series] = []

    for held in all_cohorts:
        train_mask = metadata[cohort_col] != held
        train_meta = metadata[train_mask]
        train_expr = expression.loc[train_meta.index]
        held_score: dict[str, float] = {}

        if held in cox_cohorts:
            cox_fit = fit_l1_cox(train_expr, train_meta, genes, pooled_c1=pooled_c1)
            res.fold_fits_cox[held] = cox_fit
            row = cox_fit.coefs.copy()
            row.name = held
            cox_rows.append(row)
            held_score["cox_n_held"] = int((~train_mask).sum())

        if held in log_cohorts:
            log_fit = fit_l1_logistic(train_expr, train_meta, genes)
            res.fold_fits_logistic[held] = log_fit
            row = log_fit.coefs.copy()
            row.name = held
            log_rows.append(row)
            held_score["logistic_n_held"] = int((~train_mask).sum())

        if held_score:
            res.held_out_scores[held] = held_score

    if cox_rows:
        res.coef_matrix_cox = pd.DataFrame(cox_rows)
    if log_rows:
        res.coef_matrix_logistic = pd.DataFrame(log_rows)
    return res
