"""Run the three harmonized contrasts (C1, C2, C3) across all cohorts."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from opmd_msi.contrasts.effect_size import compute_per_gene_effects
from opmd_msi.contrasts.meta_analysis import pool_across_cohorts
from opmd_msi.data.harmonize import HarmonizedDataset


@dataclass(frozen=True)
class ContrastResult:
    contrast_id: str
    per_cohort: pd.DataFrame      # long table: gene, cohort_id, estimate, variance, …
    pooled: pd.DataFrame          # one row per gene with meta-analysis stats


def load_contrasts_config(path: str | Path) -> dict:
    return yaml.safe_load(Path(path).read_text())["contrasts"]


def run_contrast(
    contrast_id: str,
    contrast_def: dict,
    harmonized: HarmonizedDataset,
    manifest_cohort_ids: list[str] | None = None,
) -> ContrastResult:
    """Compute per-cohort effects + pooled meta-analysis for one contrast."""
    expr = harmonized.expression
    meta = harmonized.metadata
    case_filter = contrast_def["case"]["filter"]
    control_filter = contrast_def["control"]["filter"]
    method = contrast_def["effect_size"]
    pairing = contrast_def.get("pairing")

    per_cohort_rows: list[pd.DataFrame] = []
    cohort_ids = manifest_cohort_ids or list(meta["cohort_id"].unique())
    for cid in cohort_ids:
        cohort_samples = meta[meta["cohort_id"] == cid].index
        if len(cohort_samples) == 0:
            continue
        sub_meta = meta.loc[cohort_samples]
        sub_expr = expr.loc[cohort_samples]
        df = compute_per_gene_effects(
            sub_expr, sub_meta,
            case_filter=case_filter,
            control_filter=control_filter,
            method=method,
            paired_by="patient_id" if pairing == "by_patient_id" else None,
        )
        if df.empty:
            continue
        df["cohort_id"] = cid
        per_cohort_rows.append(df)

    per_cohort = (pd.concat(per_cohort_rows, ignore_index=True)
                  if per_cohort_rows else pd.DataFrame())
    pooled = pool_across_cohorts(per_cohort) if not per_cohort.empty else pd.DataFrame()
    return ContrastResult(contrast_id=contrast_id, per_cohort=per_cohort, pooled=pooled)


def run_all_contrasts(
    harmonized: HarmonizedDataset,
    config_path: str | Path,
) -> dict[str, ContrastResult]:
    cfg = load_contrasts_config(config_path)
    return {cid: run_contrast(cid, cdef, harmonized) for cid, cdef in cfg.items()}
