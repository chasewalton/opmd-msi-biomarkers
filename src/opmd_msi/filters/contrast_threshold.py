"""Filter on pooled effect-size strength and cross-cohort consistency."""
from __future__ import annotations

import numpy as np
import pandas as pd


def filter_by_contrast_strength(
    pooled_c3: pd.DataFrame,
    pooled_c1: pd.DataFrame | None = None,
    pooled_c2: pd.DataFrame | None = None,
    min_abs_pooled_estimate: float = 0.50,
    max_p_fdr: float = 0.05,
    min_cohorts_same_direction: int = 4,
    require_C1_consistent_direction: bool = True,
    require_C2_consistent_direction: bool = True,
) -> tuple[list[str], pd.DataFrame]:
    """Filter genes on pooled C3 strength + optional C1/C2 direction agreement."""
    if pooled_c3.empty:
        return [], pd.DataFrame()
    df = pooled_c3.copy()
    df["abs_estimate"] = df["pooled_estimate"].abs()
    df["passes_strength"] = df["abs_estimate"] >= min_abs_pooled_estimate
    df["passes_fdr"] = df["p_fdr"] <= max_p_fdr
    df["passes_consistency"] = df["n_cohorts_same_direction"] >= min_cohorts_same_direction

    # Direction agreement with C1 / C2 (if available).
    if require_C1_consistent_direction and pooled_c1 is not None and not pooled_c1.empty:
        c1_sign = pooled_c1.set_index("gene")["pooled_estimate"].apply(np.sign)
        df["c1_sign"] = df["gene"].map(c1_sign)
        df["passes_c1_dir"] = np.sign(df["pooled_estimate"]) == df["c1_sign"]
    else:
        df["passes_c1_dir"] = True
    if require_C2_consistent_direction and pooled_c2 is not None and not pooled_c2.empty:
        c2_sign = pooled_c2.set_index("gene")["pooled_estimate"].apply(np.sign)
        df["c2_sign"] = df["gene"].map(c2_sign)
        df["passes_c2_dir"] = np.sign(df["pooled_estimate"]) == df["c2_sign"]
    else:
        df["passes_c2_dir"] = True

    df["passes_all"] = (
        df["passes_strength"]
        & df["passes_fdr"]
        & df["passes_consistency"]
        & df["passes_c1_dir"]
        & df["passes_c2_dir"]
    )
    kept = df.loc[df["passes_all"], "gene"].tolist()
    return kept, df
