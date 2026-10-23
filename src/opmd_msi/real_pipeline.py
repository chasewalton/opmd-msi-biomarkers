"""Real-data harmonization, contrasts, ENT eligibility and candidate ranking."""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from opmd_msi.contrasts.contrasts import run_all_contrasts
from opmd_msi.data.harmonize import HarmonizedDataset
from opmd_msi.data.pairing import gse37991_patient_id
from opmd_msi.data.real_annotations import (
    combined_surface_filter,
    parse_hpa,
)
from opmd_msi.ent_selection import write_ent_selection
from opmd_msi.filters.contrast_threshold import filter_by_contrast_strength
from opmd_msi.importance.composite_score import compose_ranking
from opmd_msi.importance.shap_ranking import (
    compute_correlation_with_risk,
    compute_mean_abs_shap,
)
from opmd_msi.models.validation import leave_one_cohort_out

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_cohort_parquets() -> dict[str, dict[str, pd.DataFrame]]:
    """Load every per-cohort .expression.parquet + .metadata.parquet under data/raw/."""
    raw = Path("data/raw")
    cohorts: dict[str, dict[str, pd.DataFrame]] = {}
    for ex_file in sorted(list(raw.rglob("*.expression.parquet"))):
        cid = ex_file.stem.replace(".expression", "")
        meta_file = ex_file.with_name(f"{cid}.metadata.parquet")
        if not meta_file.exists():
            continue
        metadata = pd.read_parquet(meta_file)
        if cid == "GSE37991":
            metadata["patient_id"] = metadata.patient_id.map(
                lambda value: gse37991_patient_id(value, value))
        cohorts[cid] = {
            "expression": pd.read_parquet(ex_file),
            "metadata": metadata,
        }
        log.info("loaded %s: %d samples", cid, cohorts[cid]["expression"].shape[0])
    return cohorts


# ---------------------------------------------------------------------------
# Harmonization with real ComBat
# ---------------------------------------------------------------------------

def harmonize_real(
    cohorts: dict[str, dict[str, pd.DataFrame]],
    min_cohorts_with_gene: int = 4,
) -> HarmonizedDataset:
    """Merge per-cohort matrices → common gene set → ComBat batch correction."""
    # 1. union of genes, restrict to those in >= N cohorts
    presence: dict[str, int] = {}
    for c in cohorts.values():
        for g in c["expression"].columns:
            presence[g] = presence.get(g, 0) + 1
    keep_genes = sorted(g for g, n in presence.items() if n >= min_cohorts_with_gene)
    log.info("Gene universe (>=%d cohorts): %d", min_cohorts_with_gene, len(keep_genes))

    # 2. concat
    expr_parts, meta_parts = [], []
    for _cid, c in cohorts.items():
        ex = c["expression"].reindex(columns=keep_genes)
        expr_parts.append(ex)
        meta_parts.append(c["metadata"])
    expression = pd.concat(expr_parts, axis=0)
    metadata = pd.concat(meta_parts, axis=0)

    # 3. impute NaNs with per-gene mean
    expression = expression.fillna(expression.mean(axis=0, skipna=True))
    # any remaining NaN columns (all-NaN gene) → drop
    nan_cols = expression.columns[expression.isna().any()]
    if len(nan_cols):
        log.info("dropping %d genes with all-NaN", len(nan_cols))
        expression = expression.drop(columns=nan_cols)

    # 4. ComBat batch correction across cohorts
    # combat.pycombat expects genes x samples (genes as rows). It is a
    # pure-Python implementation of the ComBat empirical-Bayes algorithm.
    try:
        from combat.pycombat import pycombat
        log.info("Running ComBat on %d samples x %d genes…",
                 expression.shape[0], expression.shape[1])
        # batch labels
        batches = metadata.loc[expression.index, "cohort_id"]
        # Only run ComBat if there is more than one batch and each batch
        # has at least 2 samples.
        bcounts = batches.value_counts()
        usable = bcounts[bcounts >= 2].index
        mask = batches.isin(usable)
        if mask.sum() < expression.shape[0]:
            log.warning("Excluding %d samples from ComBat (singleton batches)",
                        (~mask).sum())
        ex_for_combat = expression.loc[mask].T  # genes x samples
        b_for_combat = batches[mask].tolist()
        corrected = pycombat(ex_for_combat, b_for_combat)
        corrected = corrected.T   # back to samples x genes
        # write back over the in-place rows
        expression.loc[mask, :] = corrected.values
        log.info("ComBat done")
    except Exception as e:
        log.warning("ComBat failed (%s); falling back to per-cohort centering", e)
        for _cid, sample_ids in metadata.groupby("cohort_id").groups.items():
            idx = list(sample_ids)
            mean = expression.loc[idx].mean(axis=0)
            grand = expression.mean(axis=0)
            expression.loc[idx] = expression.loc[idx] - mean + grand

    return HarmonizedDataset(
        expression=expression,
        metadata=metadata,
        cohort_ids=tuple(sorted(metadata["cohort_id"].unique())),
        gene_ids=tuple(expression.columns),
    )


# ---------------------------------------------------------------------------
# Pipeline driver
# ---------------------------------------------------------------------------

def run_real_pipeline(
    out_dir: str | Path = "outputs/real_ent",
    hpa_tsv: str | Path = "data/raw/reference/HPA/proteinatlas.tsv",
    contrasts_config: str | Path = "configs/contrasts.yaml",
    filters_config: str | Path = "configs/filters.yaml",
    ent_config: str | Path | None = "configs/ent.yaml",
) -> pd.DataFrame:
    import yaml
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. Load cohorts
    cohorts = load_cohort_parquets()
    log.info("%d cohorts loaded", len(cohorts))

    # 2. Harmonize
    harm = harmonize_real(cohorts)
    harm.expression.to_parquet(out / "harmonized.expression.parquet")
    harm.metadata.to_parquet(out / "harmonized.metadata.parquet")

    # 3. Contrasts
    contrasts = run_all_contrasts(harm, contrasts_config)
    for cid, cres in contrasts.items():
        cres.per_cohort.to_parquet(out / f"contrast_{cid}_per_cohort.parquet")
        cres.pooled.to_parquet(out / f"contrast_{cid}_pooled.parquet")
        log.info("%s pooled: %d genes", cid, len(cres.pooled))

    if ent_config is not None:
        return write_ent_selection(contrasts["C1"].pooled, contrasts["C3"].pooled,
                                   out, ent_config, expression=harm.expression,
                                   metadata=harm.metadata)

    fcfg = yaml.safe_load(Path(filters_config).read_text())

    # 4. Surface localisation filter (HPA + UniProt fallback)
    hpa = parse_hpa(hpa_tsv)
    hpa.to_parquet(out / "hpa_annotations.parquet")
    surface_kept, surface_audit = combined_surface_filter(
        list(harm.gene_ids), hpa, use_uniprot_fallback=True,
    )
    surface_audit.to_parquet(out / "filter_surface_audit.parquet")
    log.info("Surface-localised genes: %d", len(surface_kept))

    # 5. Baseline-normal penalty (HPA only — Bgee/Atlas deferred to a later release)
    bfilt = fcfg["baseline_normal_penalty"]
    high_baseline = float(bfilt["high_baseline_tpm"])
    excluded_baseline = float(bfilt["excluded_baseline_tpm"])
    base_rows = []
    hpa_indexed = hpa.set_index("gene")
    baseline_kept = []
    for g in surface_kept:
        max_nontoral = float(hpa_indexed.at[g, "tissue_RNA_max_non_oral_nTPM"]) if g in hpa_indexed.index else 0.0
        penalty = 1.0 / (1.0 + np.log1p(max_nontoral)) if max_nontoral > high_baseline else 1.0
        excluded = max_nontoral > excluded_baseline
        base_rows.append({
            "gene": g, "hpa_max_nontoral_tpm": max_nontoral,
            "baseline_penalty_multiplier": penalty,
            "excluded_for_high_baseline": excluded,
        })
        if not excluded:
            baseline_kept.append(g)
    baseline_audit = pd.DataFrame(base_rows)
    baseline_audit.to_parquet(out / "filter_baseline_audit.parquet")
    log.info("After baseline-normal filter: %d", len(baseline_kept))

    # 6. Contrast-strength filter on C3
    cfilt = fcfg["lesion_to_normal_contrast"]
    pooled_c1 = contrasts.get("C1").pooled if "C1" in contrasts else pd.DataFrame()
    pooled_c2 = contrasts.get("C2").pooled if "C2" in contrasts else pd.DataFrame()
    pooled_c3 = contrasts["C3"].pooled
    contrast_kept, contrast_audit = filter_by_contrast_strength(
        pooled_c3, pooled_c1, pooled_c2,
        min_abs_pooled_estimate=float(cfilt["min_hedges_g"]),
        max_p_fdr=float(cfilt["max_meta_p_fdr"]),
        min_cohorts_same_direction=int(cfilt["min_cohorts_with_same_direction"]),
        require_C1_consistent_direction=False,   # C1 is single-cohort now
        require_C2_consistent_direction=False,   # C2 is single-cohort now
    )
    contrast_audit.to_parquet(out / "filter_contrast_audit.parquet")
    candidates = sorted(set(baseline_kept) & set(contrast_kept))
    log.info("Candidates after all filters: %d", len(candidates))
    pd.DataFrame({"gene": candidates}).to_parquet(out / "candidates.parquet")

    # 7. LOCO CV (uses real cohort labels)
    loco = leave_one_cohort_out(
        harm.expression, harm.metadata, candidates, pooled_c1=pooled_c1,
    )

    # 8. SHAP + composite
    shap_imp = compute_mean_abs_shap(loco)
    corr = compute_correlation_with_risk(harm.expression, harm.metadata, candidates)
    shap_imp.to_csv(out / "shap_importance.csv", index=False)
    corr.to_csv(out / "correlation_with_risk.csv", index=False)

    ranking = compose_ranking(
        pooled_c3=pooled_c3,
        shap_importance=shap_imp,
        corr_with_risk=corr,
        baseline_penalty=baseline_audit,
        surface_kept=candidates,
        w_effect=0.40, w_shap=0.40, w_corr=0.20,
        top_n=30,
    )
    ranking.to_csv(out / "ranking.csv", index=False)
    return ranking
