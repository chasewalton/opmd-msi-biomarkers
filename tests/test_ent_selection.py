"""ENT eligibility, annotation handling and configurable ranking."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from opmd_msi.ent_selection import (
    load_ent_config,
    prioritize_saved_run,
    rank_ent_candidates,
    select_ent_candidates,
)


def inputs():
    genes = ["EPITHELIAL", "MEMBRANE", "IMMUNE", "UNREVIEWED", "NO_SUPPLY", "NO_SIGNAL"]
    c1 = pd.DataFrame({"gene": genes, "pooled_estimate": [.8, -.7, .8, .8, .8, .1],
                       "p": [.01, .02, .01, .01, .01, .8], "p_fdr": [.1] * 6})
    c3 = pd.DataFrame({"gene": genes, "pooled_estimate": [1.2, -.1, 1.2, 1.2, 1.2, .1],
                       "p": [.001] * 6, "p_fdr": [.01, .8, .01, .01, .01, .8],
                       "n_cohorts_same_direction": [5] * 6})
    evidence = pd.DataFrame({"gene": genes,
        "epithelial_localization": ["supported", "supported", "immune", "supported", "supported", "supported"],
        "transformation_protein_evidence": ["supported", "supported", "supported", "unreviewed", "supported", "supported"],
        "source_url": ["https://example.org/study"] * 6})
    antibodies = pd.DataFrame({"gene": genes, "antibody": ["mAb"] * 6,
        "origin": ["Humanised"] * 6,
        "supplier_url": ["https://example.org/product"] * 4 + [None, "https://example.org/product"]})
    config = load_ent_config("configs/ent.yaml")
    return c1, c3, evidence, antibodies, config


def test_eligibility_distinguishes_progression_from_lesion_expression():
    c1, c3, evidence, antibodies, config = inputs()
    genes, audit = select_ent_candidates(c1, c3, evidence, antibodies, config["filters"])
    assert set(genes) == {"EPITHELIAL", "MEMBRANE"}
    audit = audit.set_index("gene")
    assert audit.at["MEMBRANE", "passes_progression"]
    assert not audit.at["MEMBRANE", "passes_lesion"]
    assert not audit.at["IMMUNE", "passes_epithelial"]
    assert not audit.at["UNREVIEWED", "passes_protein_outcome"]
    assert not audit.at["NO_SUPPLY", "passes_antibody"]
    assert not audit.at["NO_SIGNAL", "eligible"]


def test_qualifying_new_candidate_is_not_restricted_to_a_fixed_gene_list():
    c1, c3, evidence, antibodies, config = inputs()
    evidence.loc[evidence.gene.eq("UNREVIEWED"), "transformation_protein_evidence"] = "supported"
    genes, _ = select_ent_candidates(c1, c3, evidence, antibodies, config["filters"])
    assert "UNREVIEWED" in genes


def test_unavailable_and_unreviewed_annotations_fail_closed():
    c1, c3, evidence, antibodies, config = inputs()
    evidence = evidence[evidence.gene.ne("EPITHELIAL")]
    antibodies.loc[antibodies.gene.eq("MEMBRANE"), "origin"] = "Genetically human"
    genes, _ = select_ent_candidates(c1, c3, evidence, antibodies, config["filters"])
    assert genes == []


def test_ranking_weights_and_missing_background_are_explicit():
    c1, c3, evidence, antibodies, config = inputs()
    _, audit = select_ent_candidates(c1, c3, evidence, antibodies, config["filters"])
    background = pd.DataFrame({"gene": ["MEMBRANE"], "level": ["Not detected"],
                               "source_url": ["https://example.org/background"]})
    settings = dict(config["ranking"], progression_weight=0., lesion_weight=0., oral_background_weight=1.)
    ranked = rank_ent_candidates(audit, background, settings)
    assert ranked.iloc[0].gene == "MEMBRANE"
    assert ranked.iloc[0].composite_score == 1.
    missing = ranked.set_index("gene").loc["EPITHELIAL"]
    assert not missing.oral_background_observed
    assert np.isnan(missing.oral_background_score)


def test_invalid_weights_and_duplicate_evidence_are_rejected():
    c1, c3, evidence, antibodies, config = inputs()
    with pytest.raises(ValueError, match="one row per gene"):
        select_ent_candidates(c1, c3, pd.concat([evidence, evidence]), antibodies, config["filters"])
    _, audit = select_ent_candidates(c1, c3, evidence, antibodies, config["filters"])
    background = pd.DataFrame(columns=["gene", "level", "source_url"])
    with pytest.raises(ValueError, match="sum to one"):
        rank_ent_candidates(audit, background, dict(config["ranking"], progression_weight=2.))


def test_saved_input_run_cannot_be_overwritten(tmp_path):
    with pytest.raises(ValueError, match="separate"):
        prioritize_saved_run(tmp_path, tmp_path)
