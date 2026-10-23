"""Real-data cohort loading and paired contrast regression checks."""
from __future__ import annotations

import pandas as pd

from opmd_msi.contrasts.contrasts import load_contrasts_config, run_contrast
from opmd_msi.data.harmonize import HarmonizedDataset
from opmd_msi.data.pairing import gse37991_patient_id
from opmd_msi.real_pipeline import load_cohort_parquets


def test_gse37991_tumor_and_adjacent_titles_share_patient():
    assert gse37991_patient_id("T17 OSCC", "sample_T") == "GSE37991_17"
    assert gse37991_patient_id("N17 non-OSCC", "sample_N") == "GSE37991_17"
    assert gse37991_patient_id("T18", "sample_other") == "GSE37991_18"


def test_unknown_and_previously_normalized_ids_are_preserved():
    assert gse37991_patient_id("unrecognized", "sample_N") == "sample_N"
    assert gse37991_patient_id("GSE37991_17", "GSE37991_17") == "GSE37991_17"


def test_cached_tumor_adjacent_pairs_contribute_to_c2(tmp_path, monkeypatch):
    definition = load_contrasts_config("configs/contrasts.yaml")["C2"]
    ids = pd.Index(["tumor_1", "tumor_2", "normal_1", "normal_2"], name="sample_id")
    expression = pd.DataFrame({"TARGET": [3., 5., 2., 3.]}, index=ids)
    metadata = pd.DataFrame({
        "cohort_id": ["GSE37991"] * 4,
        "patient_id": ["T1", "T2", "N1", "N2"],
        "tissue_class": ["OSCC", "OSCC", "adjacent_mucosa", "adjacent_mucosa"],
    }, index=ids)
    raw = tmp_path / "data/raw/GEO/GSE37991"
    raw.mkdir(parents=True)
    expression.to_parquet(raw / "GSE37991.expression.parquet")
    metadata.to_parquet(raw / "GSE37991.metadata.parquet")
    monkeypatch.chdir(tmp_path)
    loaded = load_cohort_parquets()["GSE37991"]
    harmonized = HarmonizedDataset(
        expression=loaded["expression"], metadata=loaded["metadata"],
        cohort_ids=("GSE37991",), gene_ids=("TARGET",),
    )
    result = run_contrast("C2", definition, harmonized)
    assert result.per_cohort.n_case.tolist() == [2]
    assert result.pooled.pooled_estimate.tolist() == [1.5]
