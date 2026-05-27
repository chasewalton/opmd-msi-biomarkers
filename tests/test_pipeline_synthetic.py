"""End-to-end test: synthetic pipeline recovers PD-L1, MET, NOTCH1 in the top 3.

This is the reproducibility anchor referenced by REVIEWER_RESPONSE.md.
"""
from __future__ import annotations

from opmd_msi.pipeline import run_pipeline_synthetic


def test_pipeline_recovers_expected_top3():
    res = run_pipeline_synthetic(seed=42)
    top3 = set(res.ranking.head(3)["gene"])
    assert {"CD274", "MET", "NOTCH1"}.issubset(top3), (
        f"Expected top-3 to include PD-L1 (CD274), MET, and NOTCH1. "
        f"Got: {res.ranking.head(10)[['rank','gene','composite_score']].to_string(index=False)}"
    )


def test_pipeline_drops_non_surface_genes():
    res = run_pipeline_synthetic(seed=42)
    intracellular = {"TP53", "MKI67", "TOP2A", "BIRC5", "FOXM1", "VIM", "ACTB"}
    surface_kept_set = set(res.surface_kept)
    candidate_set = set(res.candidate_genes)
    assert intracellular.isdisjoint(candidate_set)
    assert intracellular.isdisjoint(surface_kept_set)


def test_pipeline_audit_outputs_exist():
    res = run_pipeline_synthetic(seed=42)
    assert not res.surface_audit.empty
    assert not res.baseline_audit.empty
    assert not res.contrast_filter_audit.empty
    assert "composite_score" in res.ranking.columns
    assert "abs_pooled_c3" in res.ranking.columns
    assert "mean_abs_shap" in res.ranking.columns


def test_pipeline_is_deterministic_under_seed():
    a = run_pipeline_synthetic(seed=42).ranking[["gene", "composite_score"]]
    b = run_pipeline_synthetic(seed=42).ranking[["gene", "composite_score"]]
    assert a.equals(b)
