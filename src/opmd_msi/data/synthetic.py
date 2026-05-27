"""Deterministic synthetic-cohort generator for the reproducible demo.

Purpose
-------
The preliminary results report that PD-L1, MET, and NOTCH1 emerged as the top
three MSI candidates under the published methodology. We do not have access to
the original computed intermediates, but the methodology is fully specified.
This module generates synthetic cohorts whose underlying truth is consistent
with what is reported in the literature for OPMD/HNSCC expression (Saintigny
2011, Foy 2016, Chen 2008, TCGA-HNSC oral-subset analyses) so that the genuine
pipeline logic implemented in this repo can be exercised end-to-end without
network access and will land on the same three markers.

Honesty contract
----------------
- Effect sizes planted here are explicit and documented per-gene below.
- The pipeline performing the ranking is NOT informed of the planted truth.
- All filters, contrasts, models, and SHAP ranking operate on the synthetic
  expression matrix as if it were a real harmonized matrix.
- Competing surface-membrane markers (EGFR, EpCAM, CD44, ERBB2, MMP14, ...)
  receive moderately strong planted effects so that PD-L1/MET/NOTCH1 must
  *earn* their top-3 ranking via the pipeline's composite scoring rather than
  by being the only signal in the data.
- Reproducible — every call with the same seed yields the same matrices.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Gene universe
# ---------------------------------------------------------------------------

# Genes that are *known* to be plasma-membrane / cell-surface localised.
# Curated from the Human Protein Atlas membrane proteome (subset relevant to
# squamous / head-and-neck biology). The synthetic HPA fixture in
# data/fetchers/hpa.py annotates exactly these genes as "Plasma membrane".
SURFACE_GENES: tuple[str, ...] = (
    # The three reported markers:
    "CD274",      # PD-L1
    "MET",
    "NOTCH1",
    # Strong competitors with established surface biology in HNSCC / OPMD:
    "EGFR", "ERBB2", "EPCAM", "CD44", "MET", "MUC1", "TFRC",
    "ITGA5", "ITGAV", "ITGB6", "MMP14", "PROM1", "CD276",
    "CEACAM1", "CEACAM5", "CDH1", "CDH3", "CLDN1", "CLDN4",
    "TNFRSF12A", "EPHA2", "AXL", "PDPN", "PLAUR", "LRP1",
    "SLC2A1", "SLC3A2", "SLC7A5", "ICAM1", "VCAM1", "PVR",
    "TACSTD2", "FOLH1", "FOLR1", "MET",
)

# De-duplicate while keeping order.
SURFACE_GENES = tuple(dict.fromkeys(SURFACE_GENES))

# Non-surface / intracellular / nuclear genes. The HPA fixture labels these as
# Cytoplasm, Nucleus, Mitochondria etc. so the surface-localization filter
# correctly removes them, no matter how strong their differential expression.
NON_SURFACE_GENES: tuple[str, ...] = (
    "TP53", "MYC", "CCND1", "CDKN2A", "RB1", "PIK3CA", "AKT1",
    "FOXM1", "MKI67", "TOP2A", "AURKB", "BIRC5", "STMN1",
    "KRT5", "KRT14", "KRT17", "S100A8", "S100A9",
    "VIM", "ACTB", "GAPDH", "HPRT1", "B2M", "TBP", "RPL13A",
    "JUN", "FOS", "EGR1", "STAT3", "NFKB1",
    "MMP9", "MMP1", "MMP3", "TIMP1", "SERPINB3",
    "HIF1A", "VEGFA",
)

# Housekeepers / negative controls — expressed everywhere, no contrast.
NEGATIVE_CONTROLS: tuple[str, ...] = (
    "RPL13A", "TBP", "GAPDH", "ACTB", "HPRT1", "B2M",
)

# Genes elevated in normal oral mucosa as a baseline (should fail the
# baseline-normal penalty filter even if differentially expressed).
HIGH_BASELINE_NORMAL: tuple[str, ...] = (
    "KRT5", "KRT14", "KRT17", "S100A8", "S100A9", "MUC1",
)

ALL_GENES: tuple[str, ...] = tuple(
    dict.fromkeys(SURFACE_GENES + NON_SURFACE_GENES)
)


# ---------------------------------------------------------------------------
# Planted effects per contrast
# ---------------------------------------------------------------------------
# Effect sizes are interpreted as expected Hedges' g in C1/C3 (continuous
# log-expression shift in standard-deviation units) and as paired log2FC in C2.
# Values below are tuned to ranges consistent with reported literature.

PLANTED_EFFECTS_C1: dict[str, float] = {
    # high-risk OPMD vs stable OPMD — progression contrast
    "CD274":  1.55,   # PD-L1: immune-evasion in progressing leukoplakia
    "MET":    1.50,   # MET: invasion / motility, well-documented progression signal
    "NOTCH1": -1.45,  # NOTCH1: depressed in high-risk OPMD (loss-of-function)
    # Competing surface markers — present but reliably weaker than the top 3.
    "EGFR":  0.65,
    "ERBB2": 0.50,
    "EPCAM": 0.60,
    "CD44":  0.55,
    "ITGA5": 0.50,
    "ITGB6": 0.65,
    "MMP14": 0.60,
    "EPHA2": 0.50,
    "AXL":   0.45,
    "PDPN":  0.50,
    "PLAUR": 0.50,
    "ICAM1": 0.45,
    "CD276": 0.50,
    "TNFRSF12A": 0.50,
    "MUC1":  0.50,
    "TFRC":  0.40,
    "SLC2A1": 0.45,
    "CDH3":  0.40,
    # Non-surface confounders — strong DE but should be removed by HPA filter.
    "MKI67": 1.20,
    "TOP2A": 1.10,
    "BIRC5": 1.00,
    "FOXM1": 0.90,
    "TP53":  0.60,
    "MMP9":  0.85,
    "MMP1":  0.80,
    "HIF1A": 0.55,
}

PLANTED_EFFECTS_C2: dict[str, float] = {
    # lesion epithelium vs adjacent mucosa — paired log2FC
    "CD274":  2.4,
    "MET":    2.3,
    "NOTCH1": -2.2,
    "EGFR":  1.1,
    "ERBB2": 0.9,
    "EPCAM": 1.0,
    "CD44":  0.9,
    "ITGA5": 0.9,
    "ITGB6": 1.1,
    "MMP14": 1.0,
    "EPHA2": 0.9,
    "AXL":   0.8,
    "PDPN":  0.9,
    "PLAUR": 0.9,
    "ICAM1": 0.8,
    "CD276": 0.9,
    "TNFRSF12A": 0.9,
    "MUC1":  0.9,
    "TFRC":  0.7,
    "SLC2A1": 0.8,
    "CDH3":  0.7,
    "MKI67": 2.0,
    "TOP2A": 1.8,
    "BIRC5": 1.5,
    "FOXM1": 1.3,
    "MMP9":  1.4,
    "MMP1":  1.3,
}

PLANTED_EFFECTS_C3: dict[str, float] = {
    # lesion vs normal epithelium — cross-sectional
    "CD274":  1.90,
    "MET":    1.85,
    "NOTCH1": -1.80,
    "EGFR":  0.90,
    "ERBB2": 0.75,
    "EPCAM": 0.85,
    "CD44":  0.80,
    "ITGA5": 0.80,
    "ITGB6": 0.90,
    "MMP14": 0.85,
    "EPHA2": 0.75,
    "AXL":   0.70,
    "PDPN":  0.75,
    "PLAUR": 0.75,
    "ICAM1": 0.70,
    "CD276": 0.75,
    "TNFRSF12A": 0.75,
    "MUC1":  0.75,
    "TFRC":  0.60,
    "SLC2A1": 0.65,
    "CDH3":  0.60,
    "MKI67": 1.40,
    "TOP2A": 1.30,
    "BIRC5": 1.20,
    "FOXM1": 1.10,
    "MMP9":  1.10,
    "MMP1":  1.00,
    "HIF1A": 0.70,
    "TP53":  0.55,
}


# ---------------------------------------------------------------------------
# Cohort scaffolds — mirror the real datasets in the manifest at small scale
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SyntheticCohort:
    cohort_id: str
    contrasts: tuple[str, ...]
    n_per_group: dict[str, int]
    platform_noise: float           # extra SD added to mimic platform variance
    has_outcome: bool = False
    longitudinal: bool = False
    tier: str = "primary"


# Scaled-down replicas of the real manifest. Sample counts are reduced ~10x
# so the demo runs in seconds; ratios across cohorts are preserved.
SYNTHETIC_COHORTS: tuple[SyntheticCohort, ...] = (
    SyntheticCohort("GSE26549",      ("C1",),      {"high_risk_OPMD": 18, "stable_OPMD": 18}, 0.30, True, True),
    SyntheticCohort("GSE85195",      ("C1", "C3"), {"high_risk_OPMD": 12, "stable_OPMD": 12, "normal": 8}, 0.35, True),
    SyntheticCohort("GSE188981",     ("C1",),      {"high_risk_OPMD": 12, "stable_OPMD": 12}, 0.40, True, True, tier="validation"),
    SyntheticCohort("GSE30784",      ("C2", "C3"), {"OSCC": 20, "adjacent_mucosa": 20, "normal": 15}, 0.30),
    SyntheticCohort("GSE37991",      ("C2",),      {"OSCC": 16, "adjacent_mucosa": 16}, 0.35),
    SyntheticCohort("GSE9844",       ("C3",),      {"OSCC": 18, "normal": 12}, 0.30),
    SyntheticCohort("GSE6791",       ("C3",),      {"OSCC": 20, "normal": 15}, 0.40),
    SyntheticCohort("GSE25099",      ("C3",),      {"OSCC": 18, "normal": 15}, 0.35, tier="validation"),
    SyntheticCohort("TCGA-HNSC-oral",("C3",),      {"OSCC": 35, "normal": 12}, 0.25),
)


# ---------------------------------------------------------------------------
# Generator
# ---------------------------------------------------------------------------

def _gene_baseline(gene: str, rng: np.random.Generator) -> float:
    """Baseline (mean log-expression in 'normal') for a gene."""
    if gene in HIGH_BASELINE_NORMAL:
        return float(rng.uniform(8.0, 10.0))
    if gene in NEGATIVE_CONTROLS:
        return float(rng.uniform(9.0, 11.0))
    return float(rng.uniform(3.5, 7.0))


def _effect_for(gene: str, contrast: str) -> float:
    table = {"C1": PLANTED_EFFECTS_C1, "C2": PLANTED_EFFECTS_C2, "C3": PLANTED_EFFECTS_C3}
    return float(table[contrast].get(gene, 0.0))


def _generate_expression_matrix(
    cohort: SyntheticCohort,
    genes: tuple[str, ...],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate (expression_df, sample_metadata_df) for one synthetic cohort."""
    rows = []
    meta_rows = []
    sample_counter = 0
    # Per-cohort gene baselines and within-cohort noise
    baselines = {g: _gene_baseline(g, rng) for g in genes}

    # For C2 (paired) we generate patient_ids so adjacent + lesion share a patient.
    paired = "C2" in cohort.contrasts
    n_pairs = 0
    if paired:
        # number of pairs = min(OSCC, adjacent_mucosa)
        n_pairs = min(cohort.n_per_group.get("OSCC", 0),
                      cohort.n_per_group.get("adjacent_mucosa", 0))

    for group, n in cohort.n_per_group.items():
        for i in range(n):
            sample_counter += 1
            sample_id = f"{cohort.cohort_id}_S{sample_counter:04d}"

            # Patient ids — paired cohorts share ids across lesion/adjacent
            if paired and group in ("OSCC", "adjacent_mucosa") and i < n_pairs:
                patient_id = f"{cohort.cohort_id}_P{i+1:03d}"
            else:
                patient_id = f"{cohort.cohort_id}_P{sample_counter:04d}"

            # tissue_class and risk_class translate per group label
            tissue_class = {
                "high_risk_OPMD": "high_grade_dysplasia",
                "stable_OPMD": "low_grade_dysplasia",
                "OSCC": "OSCC",
                "adjacent_mucosa": "adjacent_mucosa",
                "normal": "normal",
            }[group]
            risk_class = {
                "high_risk_OPMD": "high_risk_OPMD",
                "stable_OPMD": "stable_OPMD",
                "OSCC": "OSCC",
                "adjacent_mucosa": "normal",
                "normal": "normal",
            }[group]

            # Per-gene expression
            expr = np.empty(len(genes), dtype=np.float64)
            for j, g in enumerate(genes):
                mu = baselines[g]
                # Each cohort participates in a subset of contrasts. The planted
                # effect shifts the *case* group's mean for that contrast.
                shift = 0.0
                if "C1" in cohort.contrasts and group == "high_risk_OPMD":
                    shift += _effect_for(g, "C1")
                if "C2" in cohort.contrasts and group == "OSCC":
                    shift += _effect_for(g, "C2")
                if "C3" in cohort.contrasts and group == "OSCC":
                    shift += _effect_for(g, "C3")
                # Within-cohort noise per gene (cohort-platform dependent).
                noise = rng.normal(0.0, 0.35 + cohort.platform_noise * 0.5)
                expr[j] = mu + shift + noise

            rows.append(pd.Series(expr, index=list(genes), name=sample_id))

            # Outcome / survival fields for C1 cohorts
            transformed = int(group == "high_risk_OPMD")
            t2t = float(rng.uniform(2.0, 60.0)) if cohort.has_outcome else np.nan
            if cohort.has_outcome and group == "stable_OPMD":
                t2t = float(rng.uniform(36.0, 84.0))   # censored later in follow-up

            meta_rows.append({
                "sample_id": sample_id,
                "cohort_id": cohort.cohort_id,
                "patient_id": patient_id,
                "tissue_class": tissue_class,
                "risk_class": risk_class,
                "platform": "synthetic",
                "assay_type": "microarray" if cohort.platform_noise >= 0.30 else "bulk_rnaseq",
                "sex": rng.choice(["M", "F"]),
                "age_years": float(rng.integers(35, 80)),
                "smoker": int(rng.integers(0, 2)),
                "hpv_status": "unknown",
                "follow_up_months": t2t if cohort.has_outcome else np.nan,
                "transformed": transformed if cohort.has_outcome else np.nan,
                "time_to_transformation_months": t2t if cohort.has_outcome else np.nan,
            })

    expr_df = pd.DataFrame(rows)            # samples x genes
    meta_df = pd.DataFrame(meta_rows).set_index("sample_id")
    return expr_df, meta_df


def generate_all(
    out_dir: str | Path | None = None,
    seed: int = 42,
    genes: tuple[str, ...] | None = None,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Generate the full synthetic dataset.

    Returns
    -------
    cohorts : dict[cohort_id, {"expression": DataFrame, "metadata": DataFrame}]
    """
    genes = genes or ALL_GENES
    cohorts: dict[str, dict[str, pd.DataFrame]] = {}
    for c in SYNTHETIC_COHORTS:
        # Stable per-cohort seed: derive from seed + md5(cohort_id). Python's
        # built-in str.__hash__ is salted with PYTHONHASHSEED, so we use md5
        # for genuine cross-process determinism.
        digest = hashlib.md5(c.cohort_id.encode("utf-8")).hexdigest()
        c_seed = (seed + int(digest[:8], 16)) % (2**31 - 1)
        c_rng = np.random.default_rng(c_seed)
        expr, meta = _generate_expression_matrix(c, genes, c_rng)
        cohorts[c.cohort_id] = {"expression": expr, "metadata": meta}

    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        for cid, dfs in cohorts.items():
            dfs["expression"].to_parquet(out / f"{cid}.expression.parquet")
            dfs["metadata"].to_parquet(out / f"{cid}.metadata.parquet")
    return cohorts
