"""Random-effects meta-analysis across cohorts (DerSimonian-Laird)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def _dersimonian_laird(estimates: np.ndarray, variances: np.ndarray) -> dict[str, float]:
    """Random-effects pooled estimate using DerSimonian-Laird tau^2."""
    w_fixed = 1.0 / variances
    pooled_fixed = float((w_fixed * estimates).sum() / w_fixed.sum())
    Q = float(((estimates - pooled_fixed) ** 2 * w_fixed).sum())
    k = len(estimates)
    df = k - 1
    c = float(w_fixed.sum() - (w_fixed**2).sum() / w_fixed.sum())
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (variances + tau2)
    pooled = float((w_re * estimates).sum() / w_re.sum())
    se = float(np.sqrt(1.0 / w_re.sum()))
    z = pooled / se if se > 0 else 0.0
    p = float(2 * (1 - stats.norm.cdf(abs(z))))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    return {
        "pooled_estimate": pooled,
        "se": se,
        "z": z,
        "p": p,
        "tau2": tau2,
        "Q": Q,
        "I2": I2,
        "k_cohorts": k,
    }


def pool_across_cohorts(per_cohort: pd.DataFrame) -> pd.DataFrame:
    """Pool per-cohort effect sizes per gene.

    Input columns: gene, cohort_id, estimate, variance, n_case, n_control, p_value.
    Returns one row per gene with pooled stats and direction consistency.

    When only 1 cohort contributes (e.g. single-cohort C1 in the real-data
    run), the function still emits a row with the per-cohort estimate as the
    "pooled" value (k=1, tau²=0, I²=NaN) so downstream filters and rankings
    can use it. The n_cohorts_total field always reflects the real cohort count.
    """
    rows = []
    for gene, sub in per_cohort.groupby("gene"):
        est = sub["estimate"].values.astype(float)
        var = sub["variance"].values.astype(float)
        # Drop any rows with non-finite variance (e.g. zero-variance cohorts)
        mask = np.isfinite(est) & np.isfinite(var) & (var > 0)
        if mask.sum() == 0:
            continue
        if mask.sum() == 1:
            # single-cohort fallback — no meta-analysis, just propagate.
            e = float(est[mask][0])
            v = float(var[mask][0])
            se = float(np.sqrt(v))
            z = e / se if se > 0 else 0.0
            p = float(2 * (1 - stats.norm.cdf(abs(z))))
            meta = {"pooled_estimate": e, "se": se, "z": z, "p": p,
                    "tau2": 0.0, "Q": 0.0, "I2": float("nan"), "k_cohorts": 1}
            same_dir = 1
        else:
            meta = _dersimonian_laird(est[mask], var[mask])
            same_dir = int(np.sum(np.sign(est[mask]) == np.sign(meta["pooled_estimate"])))
        rows.append({
            "gene": gene,
            **meta,
            "n_cohorts_same_direction": same_dir,
            "n_cohorts_total": int(mask.sum()),
        })
    if not rows:
        return pd.DataFrame()
    pooled = pd.DataFrame(rows)
    # FDR adjust pooled p
    pooled["p_fdr"] = _benjamini_hochberg(pooled["p"].values)
    return pooled.sort_values("p_fdr").reset_index(drop=True)


def _benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = np.empty(n)
    for i, idx in enumerate(order):
        ranked[idx] = p[idx] * n / (i + 1)
    # Enforce monotonicity from largest to smallest
    sorted_p = ranked[order]
    for i in range(n - 2, -1, -1):
        sorted_p[i] = min(sorted_p[i], sorted_p[i + 1])
    out = np.empty(n)
    for i, idx in enumerate(order):
        out[idx] = min(sorted_p[i], 1.0)
    return out
