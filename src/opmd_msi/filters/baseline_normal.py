"""Penalty filter for genes with elevated baseline expression in normal tissues.

Combines three reference panels (Human Protein Atlas, Bgee, Expression Atlas)
and requires agreement of >= 2/3 panels before a gene is considered "baseline-high."
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _max_non_oral_tpm_from_hpa(hpa: pd.DataFrame, oral_tissues: tuple[str, ...]) -> pd.Series:
    """Per-gene maximum nTPM across HPA tissues that are NOT in `oral_tissues`."""
    cols = [c for c in hpa.columns
            if c.startswith("Tissue RNA")
            and not any(o in c for o in oral_tissues)]
    if not cols:
        return pd.Series(dtype=float)
    return hpa.set_index("Gene")[cols].max(axis=1)


def _panel_high(series: pd.Series, threshold: float) -> pd.Series:
    return series.fillna(0) > threshold


def compute_baseline_penalty(
    genes: list[str],
    hpa_annotations: pd.DataFrame,
    bgee_annotations: pd.DataFrame,
    atlas_annotations: pd.DataFrame,
    high_baseline_tpm: float = 15.0,
    excluded_baseline_tpm: float = 50.0,
    min_panels_in_agreement: int = 2,
    oral_relevant_tissues: tuple[str, ...] = (
        "oral mucosa", "tongue", "skin", "esophagus",
    ),
) -> tuple[list[str], pd.DataFrame]:
    """Return (kept_genes, penalty_audit).

    `penalty_audit` includes a per-gene `baseline_penalty_multiplier` in [0, 1].
    A higher value means the gene is "cleaner" — i.e., low baseline in normals.
    Genes whose baseline exceeds `excluded_baseline_tpm` in >=2 panels are
    dropped entirely.
    """
    hpa_max = _max_non_oral_tpm_from_hpa(hpa_annotations, oral_relevant_tissues)
    bgee_max = bgee_annotations.set_index("Gene")["max_non_oral_normal_tpm"]
    atlas_max = atlas_annotations.set_index("Gene")["median_tpm_across_E-MTAB-2836"]

    rows = []
    kept: list[str] = []
    for g in genes:
        hpa_val = float(hpa_max.get(g, 0.0))
        bgee_val = float(bgee_max.get(g, 0.0))
        atlas_val = float(atlas_max.get(g, 0.0))
        # Number of panels that flag the gene as high-baseline / excluded
        panels_high = sum(v > high_baseline_tpm for v in (hpa_val, bgee_val, atlas_val))
        panels_excluded = sum(v > excluded_baseline_tpm for v in (hpa_val, bgee_val, atlas_val))
        max_val = max(hpa_val, bgee_val, atlas_val)
        penalty = 1.0 / (1.0 + np.log1p(max_val)) if panels_high >= min_panels_in_agreement else 1.0
        excluded = panels_excluded >= min_panels_in_agreement
        rows.append({
            "gene": g,
            "hpa_max_nontoral_tpm": hpa_val,
            "bgee_max_nontoral_tpm": bgee_val,
            "atlas_median_tpm": atlas_val,
            "panels_high": panels_high,
            "panels_excluded": panels_excluded,
            "baseline_penalty_multiplier": penalty,
            "excluded_for_high_baseline": excluded,
        })
        if not excluded:
            kept.append(g)
    return kept, pd.DataFrame(rows)
