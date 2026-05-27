"""End-to-end pipeline orchestrator.

Stages (mirrors the numbered scripts):
  01  Acquire data (real fetch OR synthetic generator).
  02  Harmonize across cohorts.
  03  Compute 3 contrasts (C1, C2, C3) + meta-analysis.
  04  Apply filters (surface localization, baseline-normal, contrast strength).
  05  Train models (L1 Cox, L1 logistic) under leave-one-cohort-out CV.
  06  Compute SHAP importance + composite score → ranked candidate list.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from opmd_msi.contrasts.contrasts import run_all_contrasts
from opmd_msi.data.fetchers.atlas import synthetic_atlas_annotations
from opmd_msi.data.fetchers.bgee import synthetic_bgee_annotations
from opmd_msi.data.fetchers.hpa import synthetic_hpa_annotations
from opmd_msi.data.harmonize import HarmonizedDataset, harmonize
from opmd_msi.data.synthetic import generate_all
from opmd_msi.filters.baseline_normal import compute_baseline_penalty
from opmd_msi.filters.contrast_threshold import filter_by_contrast_strength
from opmd_msi.filters.surface_localization import filter_surface_localized
from opmd_msi.importance.composite_score import compose_ranking
from opmd_msi.importance.shap_ranking import (
    compute_correlation_with_risk,
    compute_mean_abs_shap,
)
from opmd_msi.models.validation import leave_one_cohort_out


@dataclass
class PipelineResult:
    harmonized: HarmonizedDataset
    contrasts: dict
    surface_kept: list[str]
    surface_audit: pd.DataFrame
    baseline_kept: list[str]
    baseline_audit: pd.DataFrame
    contrast_filter_audit: pd.DataFrame
    candidate_genes: list[str]
    loco_result: object
    shap_importance: pd.DataFrame
    corr_with_risk: pd.DataFrame
    ranking: pd.DataFrame


def run_pipeline_synthetic(
    config_path: str | Path = "configs/default.yaml",
    contrasts_config_path: str | Path = "configs/contrasts.yaml",
    filters_config_path: str | Path = "configs/filters.yaml",
    seed: int = 42,
    out_dir: str | Path | None = None,
) -> PipelineResult:
    """Run the end-to-end pipeline on synthetic cohorts."""
    cfg = yaml.safe_load(Path(config_path).read_text())
    fcfg = yaml.safe_load(Path(filters_config_path).read_text())

    # 1. Generate synthetic data
    cohorts = generate_all(seed=seed)

    # 2. Harmonize
    harm_cfg = cfg.get("harmonization", {})
    harmonized = harmonize(
        cohorts,
        drop_genes_in_fewer_than_n_cohorts=harm_cfg.get(
            "drop_genes_in_fewer_than_n_cohorts", 4),
        quantile_normalize=harm_cfg.get("per_cohort_quantile_normalize", True),
        batch_correct=harm_cfg.get("remove_batch_effects", "combat") is not None,
    )

    # 3. Contrasts
    contrasts = run_all_contrasts(harmonized, contrasts_config_path)
    pooled_c1 = contrasts["C1"].pooled
    pooled_c2 = contrasts["C2"].pooled
    pooled_c3 = contrasts["C3"].pooled

    # 4. Filters
    hpa = synthetic_hpa_annotations()
    bgee = synthetic_bgee_annotations()
    atlas = synthetic_atlas_annotations()

    sfilt = fcfg["surface_localization"]
    surface_kept, surface_audit = filter_surface_localized(
        list(harmonized.gene_ids),
        hpa,
        required_locations=tuple(sfilt["required_subcellular_locations"]),
        required_secretome=tuple(sfilt.get("required_secretome_categories", ())),
        allow_secreted=sfilt.get("allow_secreted", False),
        allow_predicted_only=sfilt.get("allow_predicted_only", False),
    )

    bfilt = fcfg["baseline_normal_penalty"]
    baseline_kept, baseline_audit = compute_baseline_penalty(
        surface_kept, hpa, bgee, atlas,
        high_baseline_tpm=bfilt["high_baseline_tpm"],
        excluded_baseline_tpm=bfilt["excluded_baseline_tpm"],
        min_panels_in_agreement=bfilt["min_panels_in_agreement"],
        oral_relevant_tissues=tuple(bfilt.get("oral_relevant_tissues", ())),
    )

    cfilt = fcfg["lesion_to_normal_contrast"]
    contrast_kept, contrast_filter_audit = filter_by_contrast_strength(
        pooled_c3,
        pooled_c1=pooled_c1,
        pooled_c2=pooled_c2,
        min_abs_pooled_estimate=cfilt["min_hedges_g"],
        max_p_fdr=cfilt["max_meta_p_fdr"],
        min_cohorts_same_direction=cfilt["min_cohorts_with_same_direction"],
        require_C1_consistent_direction=cfilt["require_C1_consistent_direction"],
        require_C2_consistent_direction=cfilt["require_C2_consistent_direction"],
    )
    candidate_genes = sorted(set(baseline_kept) & set(contrast_kept))

    # 5. Models — LOCO CV
    loco = leave_one_cohort_out(
        harmonized.expression, harmonized.metadata, candidate_genes,
        pooled_c1=pooled_c1,
    )

    # 6. SHAP + composite scoring
    shap_imp = compute_mean_abs_shap(loco)
    corr = compute_correlation_with_risk(
        harmonized.expression, harmonized.metadata, candidate_genes,
    )
    rcfg = cfg.get("ranking", {}).get("composite_score", {})
    ranking = compose_ranking(
        pooled_c3=pooled_c3,
        shap_importance=shap_imp,
        corr_with_risk=corr,
        baseline_penalty=baseline_audit,
        surface_kept=candidate_genes,
        w_effect=rcfg.get("w_effect", 0.40),
        w_shap=rcfg.get("w_shap", 0.40),
        w_corr=rcfg.get("w_corr", 0.20),
        top_n=cfg.get("ranking", {}).get("top_n", 25),
    )

    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        ranking.to_csv(out / "ranking.csv", index=False)
        surface_audit.to_csv(out / "filter_surface_localization_audit.csv", index=False)
        baseline_audit.to_csv(out / "filter_baseline_normal_audit.csv", index=False)
        contrast_filter_audit.to_csv(out / "filter_contrast_strength_audit.csv", index=False)
        for cid, cres in contrasts.items():
            cres.per_cohort.to_csv(out / f"contrast_{cid}_per_cohort.csv", index=False)
            cres.pooled.to_csv(out / f"contrast_{cid}_pooled.csv", index=False)
        shap_imp.to_csv(out / "shap_importance.csv", index=False)
        corr.to_csv(out / "correlation_with_risk.csv", index=False)

    return PipelineResult(
        harmonized=harmonized,
        contrasts=contrasts,
        surface_kept=surface_kept,
        surface_audit=surface_audit,
        baseline_kept=baseline_kept,
        baseline_audit=baseline_audit,
        contrast_filter_audit=contrast_filter_audit,
        candidate_genes=candidate_genes,
        loco_result=loco,
        shap_importance=shap_imp,
        corr_with_risk=corr,
        ranking=ranking,
    )
