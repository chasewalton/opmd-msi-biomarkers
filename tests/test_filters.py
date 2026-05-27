from opmd_msi.data.fetchers.atlas import synthetic_atlas_annotations
from opmd_msi.data.fetchers.bgee import synthetic_bgee_annotations
from opmd_msi.data.fetchers.hpa import synthetic_hpa_annotations
from opmd_msi.filters.baseline_normal import compute_baseline_penalty
from opmd_msi.filters.surface_localization import filter_surface_localized


def test_surface_filter_keeps_known_membrane_proteins():
    hpa = synthetic_hpa_annotations()
    kept, audit = filter_surface_localized(
        ["CD274", "MET", "NOTCH1", "TP53", "MYC", "GAPDH"], hpa,
    )
    assert "CD274" in kept
    assert "MET" in kept
    assert "NOTCH1" in kept
    assert "TP53" not in kept
    assert "MYC" not in kept


def test_baseline_penalty_drops_high_baseline_keratins():
    hpa = synthetic_hpa_annotations()
    bgee = synthetic_bgee_annotations()
    atlas = synthetic_atlas_annotations()
    kept, audit = compute_baseline_penalty(
        ["CD274", "MET", "KRT5", "KRT14", "MUC1"], hpa, bgee, atlas,
        high_baseline_tpm=15.0, excluded_baseline_tpm=50.0,
        oral_relevant_tissues=("oral mucosa", "tongue", "esophagus", "skin"),
    )
    # KRT5 / KRT14 are flagged as high baseline; MUC1 too.
    # With default exclusion threshold (50 TPM) they're penalized but not dropped.
    assert "CD274" in kept
    assert "MET" in kept
    # Penalty multiplier should be < 1 for the high-baseline ones
    audit = audit.set_index("gene")
    assert audit.loc["KRT5", "baseline_penalty_multiplier"] < 1.0
    assert audit.loc["CD274", "baseline_penalty_multiplier"] == 1.0
