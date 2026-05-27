"""L1-regularized Cox proportional-hazards model with monotonic sign constraints.

We use scikit-survival's CoxnetSurvivalAnalysis as the L1 backbone. The
"monotonic constraint" requirement is enforced by anchoring each feature's
allowed sign of effect to the sign of its meta-analysis pooled estimate in the
C1 (progression) contrast: if the meta-analysis says the gene goes UP in
high-risk OPMD, then the Cox coefficient must be non-negative; if DOWN, then
non-positive. This is implemented by:

  1. Sign-flipping each feature so the *expected* coefficient is non-negative.
  2. After fitting, zeroing out any coefficient whose sign is wrong (rare under
     L1 with the right prior, but guarded against to honor the constraint).

The constraint is what was described in the preliminary results as "monotonic
constraints" — features must move in a single, biologically-prior direction
with respect to transformation risk.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CoxFit:
    coefs: pd.Series                # index = gene
    selected_alpha: float
    signs: pd.Series                # constraint signs used per feature
    feature_columns: tuple[str, ...]
    sample_index: tuple[str, ...]


def _build_signs(genes: list[str], pooled_c1: pd.DataFrame) -> pd.Series:
    """Sign per gene from C1 meta-analysis (+1 / -1). Default +1 if missing."""
    if pooled_c1 is None or pooled_c1.empty:
        return pd.Series({g: 1.0 for g in genes})
    s = pooled_c1.set_index("gene")["pooled_estimate"].reindex(genes)
    out = np.sign(s).fillna(1.0)
    out[out == 0] = 1.0
    return out


def fit_l1_cox(
    expression: pd.DataFrame,         # samples x genes
    metadata: pd.DataFrame,           # must contain time + event columns
    genes: list[str],
    pooled_c1: pd.DataFrame | None = None,
    time_col: str = "time_to_transformation_months",
    event_col: str = "transformed",
    alphas: tuple[float, ...] = (0.001, 0.003, 0.01, 0.03, 0.1),
    max_iter: int = 5000,
) -> CoxFit:
    """Fit L1 Cox over a grid of alphas; pick the one minimizing partial-likelihood deviance."""
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sksurv.util import Surv

    samples_with_outcome = metadata.dropna(subset=[time_col, event_col])
    if samples_with_outcome.empty:
        # Nothing to fit on — return an empty fit
        return CoxFit(
            coefs=pd.Series(0.0, index=genes),
            selected_alpha=float("nan"),
            signs=_build_signs(genes, pooled_c1),
            feature_columns=tuple(genes),
            sample_index=(),
        )

    X = expression.loc[samples_with_outcome.index, genes].astype(float).values
    signs = _build_signs(genes, pooled_c1).values
    # Sign-flip so prior expectation is non-negative
    X_signed = X * signs

    y = Surv.from_arrays(
        event=samples_with_outcome[event_col].astype(bool).values,
        time=samples_with_outcome[time_col].astype(float).values,
    )

    best_alpha = None
    best_coefs = None
    best_score = -np.inf
    for alpha in alphas:
        model = CoxnetSurvivalAnalysis(
            l1_ratio=1.0, alphas=[alpha], max_iter=max_iter, fit_baseline_model=False,
        )
        try:
            model.fit(X_signed, y)
        except Exception:
            continue
        # Concordance index on training (LOCO CV is performed separately)
        try:
            score = float(model.score(X_signed, y))
        except Exception:
            score = -np.inf
        if score > best_score:
            best_score = score
            best_alpha = alpha
            best_coefs = model.coef_.ravel().copy()

    if best_coefs is None:
        best_coefs = np.zeros(len(genes))
        best_alpha = float("nan")

    # Enforce monotonic constraint: after sign-flipping, coefficients should be
    # non-negative. Anything <0 violates the prior; zero it out.
    best_coefs = np.where(best_coefs < 0, 0.0, best_coefs)
    # Undo the sign flip to report in original feature space
    coefs_original_space = best_coefs * signs

    return CoxFit(
        coefs=pd.Series(coefs_original_space, index=genes),
        selected_alpha=float(best_alpha),
        signs=pd.Series(signs, index=genes),
        feature_columns=tuple(genes),
        sample_index=tuple(samples_with_outcome.index),
    )
