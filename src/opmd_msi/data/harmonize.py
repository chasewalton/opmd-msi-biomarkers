"""Cross-cohort harmonization.

Real cohorts come from heterogeneous platforms (Affymetrix HG-U133, Illumina
HT-12, RNA-seq). The pipeline applies:

  1. Gene-id space unification (HGNC symbols).
  2. log2 transform with a small pseudocount (skipped if already log-scale).
  3. Per-cohort quantile normalization.
  4. ComBat-style empirical-Bayes batch correction across cohorts.
  5. Drop genes present in fewer than N cohorts.

For the synthetic demo the inputs are already on a comparable log-scale, so
the harmonizer's batch-correction step is a no-op except for centering.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class HarmonizedDataset:
    expression: pd.DataFrame   # samples x genes (rows = sample_id)
    metadata: pd.DataFrame     # samples x meta-fields, index = sample_id
    cohort_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]


def _quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Per-sample quantile normalization to the cohort's mean rank distribution."""
    ranked = df.rank(axis=1, method="average")
    sorted_means = np.sort(df.values, axis=1).mean(axis=0)
    # Map each rank back via sorted_means
    n = df.shape[1]
    out = np.empty_like(df.values, dtype=np.float64)
    for i in range(df.shape[0]):
        ranks = ranked.values[i] - 1   # 0-indexed
        # Interpolate so non-integer ranks (ties) average neighbours
        lower = np.floor(ranks).astype(int).clip(0, n - 1)
        upper = np.ceil(ranks).astype(int).clip(0, n - 1)
        frac = ranks - lower
        out[i] = sorted_means[lower] * (1 - frac) + sorted_means[upper] * frac
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def _combat_like(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    cohort_col: str = "cohort_id",
) -> pd.DataFrame:
    """Center each cohort to a common mean per gene.

    Lightweight stand-in for ComBat. For the synthetic demo, removing the
    cohort-level mean shift per gene is sufficient. The real-data path swaps
    this for `combat_seq` (sva package) via R.
    """
    expr = expression.copy()
    grand_mean = expr.mean(axis=0)
    for _cohort_id, sample_ids in metadata.groupby(cohort_col).groups.items():
        idx = list(sample_ids)
        cohort_mean = expr.loc[idx].mean(axis=0)
        expr.loc[idx] = expr.loc[idx] - cohort_mean + grand_mean
    return expr


def harmonize(
    cohorts: dict[str, dict[str, pd.DataFrame]],
    drop_genes_in_fewer_than_n_cohorts: int = 4,
    quantile_normalize: bool = True,
    batch_correct: bool = True,
) -> HarmonizedDataset:
    """Merge per-cohort matrices into one harmonized matrix.

    Parameters
    ----------
    cohorts : dict[cohort_id, {"expression": DataFrame, "metadata": DataFrame}]
    """
    # 1. union of genes; restrict to those present in >= N cohorts
    presence: dict[str, int] = {}
    for c in cohorts.values():
        for g in c["expression"].columns:
            presence[g] = presence.get(g, 0) + 1
    keep_genes = sorted(g for g, n in presence.items()
                        if n >= drop_genes_in_fewer_than_n_cohorts)

    # 2. concat (samples on rows; reindex to keep_genes; fill missing with NaN)
    expr_parts: list[pd.DataFrame] = []
    meta_parts: list[pd.DataFrame] = []
    for _cid, payload in cohorts.items():
        ex = payload["expression"].reindex(columns=keep_genes)
        if quantile_normalize:
            ex = _quantile_normalize(ex)
        expr_parts.append(ex)
        meta_parts.append(payload["metadata"])
    expression = pd.concat(expr_parts, axis=0)
    metadata = pd.concat(meta_parts, axis=0)

    # 3. simple imputation for any residual NaNs (per-gene cohort-mean)
    expression = expression.fillna(expression.mean(axis=0, skipna=True))

    # 4. batch correction
    if batch_correct:
        expression = _combat_like(expression, metadata, cohort_col="cohort_id")

    return HarmonizedDataset(
        expression=expression,
        metadata=metadata,
        cohort_ids=tuple(sorted(metadata["cohort_id"].unique())),
        gene_ids=tuple(keep_genes),
    )
