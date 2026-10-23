"""ENT candidate eligibility and evidence-based ranking."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from opmd_msi.contrasts.contrasts import load_contrasts_config
from opmd_msi.contrasts.meta_analysis import pool_across_cohorts


def _require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = columns - set(frame.columns)
    if missing:
        raise ValueError(f"{name} missing columns: {sorted(missing)}")


def load_ent_config(path: str | Path) -> dict:
    path = Path(path).resolve()
    config = yaml.safe_load(path.read_text())
    for name, value in config["references"].items():
        config["references"][name] = path.parent / value
    return config


def select_ent_candidates(
    pooled_c1: pd.DataFrame,
    pooled_c3: pd.DataFrame,
    evidence: pd.DataFrame,
    antibodies: pd.DataFrame,
    filters: dict,
) -> tuple[list[str], pd.DataFrame]:
    """Apply epithelial, protein-outcome, antibody and expression eligibility."""
    _require_columns(evidence, {"gene", "epithelial_localization",
                                "transformation_protein_evidence", "source_url"}, "Evidence")
    _require_columns(antibodies, {"gene", "antibody", "origin", "supplier_url"}, "Antibodies")
    if evidence.gene.duplicated().any():
        raise ValueError("Evidence must contain one row per gene")
    for frame, label in [(pooled_c1, "C1"), (pooled_c3, "C3")]:
        _require_columns(frame, {"gene", "pooled_estimate", "p", "p_fdr"}, label)
        if frame.gene.duplicated().any():
            raise ValueError(f"{label} must contain one row per gene")
    _require_columns(pooled_c3, {"n_cohorts_same_direction"}, "C3")

    antibody_mask = antibodies.origin.eq(filters["antibody_origin"])
    if filters.get("require_supplier_listing", True):
        antibody_mask &= antibodies.supplier_url.fillna("").str.startswith("https://")
    available = antibodies.loc[antibody_mask].groupby("gene").agg(
        antibodies=("antibody", lambda values: ";".join(sorted(set(values)))),
        supplier_urls=("supplier_url", lambda values: ";".join(sorted(set(values)))),
    )
    c1 = pooled_c1.set_index("gene")[["pooled_estimate", "p", "p_fdr"]].add_prefix("c1_")
    c3 = pooled_c3.set_index("gene")[["pooled_estimate", "p_fdr",
                                         "n_cohorts_same_direction"]].add_prefix("c3_")
    audit = c3.join(c1, how="outer").join(evidence.set_index("gene")).join(available)
    audit["passes_epithelial"] = audit.epithelial_localization.eq(filters["epithelial_localization"])
    audit["passes_protein_outcome"] = (
        audit.transformation_protein_evidence.eq(filters["transformation_protein_evidence"])
        & audit.source_url.fillna("").str.startswith("https://"))
    audit["passes_antibody"] = audit.antibodies.notna()
    audit["passes_progression"] = (
        (audit.c1_pooled_estimate.abs() >= filters["min_progression_abs_g"])
        & (audit.c1_p <= filters["max_progression_p"]))
    audit["passes_lesion"] = (
        (audit.c3_pooled_estimate >= filters["min_lesion_g"])
        & (audit.c3_p_fdr <= filters["max_lesion_fdr"])
        & (audit.c3_n_cohorts_same_direction >= filters["min_consistent_cohorts"]))
    audit["eligible"] = (audit.passes_epithelial & audit.passes_protein_outcome
                         & audit.passes_antibody & (audit.passes_progression | audit.passes_lesion))
    audit.index.name = "gene"
    audit = audit.reset_index().sort_values("gene").reset_index(drop=True)
    return audit.loc[audit.eligible, "gene"].tolist(), audit


def rank_ent_candidates(audit: pd.DataFrame, oral_background: pd.DataFrame,
                        ranking_config: dict) -> pd.DataFrame:
    """Rank eligible candidates by progression, lesion enrichment and oral background."""
    _require_columns(oral_background, {"gene", "level", "source_url"}, "Oral background")
    if oral_background.gene.duplicated().any():
        raise ValueError("Oral background must contain one row per gene")
    weights = np.array([ranking_config["progression_weight"], ranking_config["lesion_weight"],
                        ranking_config["oral_background_weight"]], dtype=float)
    if not np.isfinite(weights).all() or (weights < 0).any() or not np.isclose(weights.sum(), 1):
        raise ValueError("Ranking weights must be finite, nonnegative and sum to one")
    result = audit.loc[audit.eligible].copy()
    result = result.merge(oral_background.rename(columns={"level": "normal_oral_level",
                          "source_url": "oral_background_source"}), on="gene", how="left")
    progression = result.c1_pooled_estimate.abs()
    lesion = result.c3_pooled_estimate.clip(lower=0)
    result["progression_score"] = progression.rank(pct=True).where(progression > 0, 0).fillna(0)
    result["lesion_score"] = lesion.rank(pct=True).where(lesion > 0, 0).fillna(0)
    result["oral_background_score"] = result.normal_oral_level.map(
        {"Not detected": 1., "Low": .75, "Medium": .5, "High": .25})
    result["oral_background_observed"] = result.oral_background_score.notna()
    result["composite_score"] = (weights[0] * result.progression_score
                                  + weights[1] * result.lesion_score
                                  + weights[2] * result.oral_background_score.fillna(0))
    result = result.sort_values(["composite_score", "gene"], ascending=[False, True]).reset_index(drop=True)
    result["rank"] = np.arange(1, len(result) + 1)
    return result.head(int(ranking_config.get("top_n", 25)))


def write_ent_selection(pooled_c1: pd.DataFrame, pooled_c3: pd.DataFrame,
                        out_dir: str | Path, config_path: str | Path = "configs/ent.yaml",
                        expression: pd.DataFrame | None = None,
                        metadata: pd.DataFrame | None = None) -> pd.DataFrame:
    config = load_ent_config(config_path)
    refs = config["references"]
    _, audit = select_ent_candidates(pooled_c1, pooled_c3, pd.read_csv(refs["evidence"]),
                                     pd.read_csv(refs["antibodies"]), config["filters"])
    ranking = rank_ent_candidates(audit, pd.read_csv(refs["oral_background"]), config["ranking"])
    ranking["selection_type"] = config["selection_type"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (expression is None) != (metadata is None):
        raise ValueError("Expression and metadata must be supplied together")
    if expression is not None and metadata is not None and not ranking.empty:
        from opmd_msi.importance.shap_ranking import (
            compute_correlation_with_risk,
            compute_mean_abs_shap,
        )
        from opmd_msi.models.validation import leave_one_cohort_out

        genes = audit.loc[audit.eligible, "gene"].tolist()
        loco = leave_one_cohort_out(expression, metadata, genes, pooled_c1=pooled_c1)
        shap = compute_mean_abs_shap(loco)
        correlation = compute_correlation_with_risk(expression, metadata, genes)
        shap.to_csv(out / "shap_importance.csv", index=False)
        correlation.to_csv(out / "correlation_with_risk.csv", index=False)
        loco.coef_matrix_cox.to_parquet(out / "loco_coefs_cox.parquet")
        loco.coef_matrix_logistic.to_parquet(out / "loco_coefs_logistic.parquet")
        ranking = ranking.merge(shap, on="gene", how="left").merge(correlation, on="gene", how="left")
    audit.to_csv(out / "ent_filter_audit.csv", index=False)
    ranking.to_csv(out / "ranking.csv", index=False)
    ranking.head(config["ranking"]["promoted_top_n"]).to_csv(out / "ranking_top_promoted.csv", index=False)
    return ranking


def prioritize_saved_run(run_dir: str | Path, out_dir: str | Path,
                         config_path: str | Path = "configs/ent.yaml",
                         contrasts_config: str | Path = "configs/contrasts.yaml") -> pd.DataFrame:
    """Re-pool saved effects after excluding numerically constant cohort signals."""
    run = Path(run_dir)
    out = Path(out_dir)
    if run.resolve() == out.resolve():
        raise ValueError("Choose an output directory separate from the saved input run")
    expression = pd.read_parquet(run / "harmonized.expression.parquet")
    metadata = pd.read_parquet(run / "harmonized.metadata.parquet")
    contrasts = load_contrasts_config(contrasts_config)
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for name in ["C1", "C3"]:
        per = pd.read_parquet(run / f"contrast_{name}_per_cohort.parquet")
        valid = pd.Series(False, index=per.index)
        definition = contrasts[name]
        for cohort, rows in per.groupby("cohort_id"):
            meta = metadata[metadata.cohort_id.eq(cohort)]
            ids = meta.query(definition["case"]["filter"]).index.union(
                meta.query(definition["control"]["filter"]).index)
            if ids.empty:
                continue
            values = expression.loc[ids]
            spread = values.max() - values.min()
            tolerance = 64 * np.finfo(float).eps * values.abs().max().clip(lower=1)
            valid.loc[rows.index] = rows.gene.isin(spread[spread > tolerance].index)
        clean = per.loc[valid].copy()
        results[name] = pool_across_cohorts(clean)
        clean.to_parquet(out / f"contrast_{name}_per_cohort.parquet")
        results[name].to_parquet(out / f"contrast_{name}_pooled.parquet")
    return write_ent_selection(results["C1"], results["C3"], out, config_path,
                               expression=expression, metadata=metadata)
