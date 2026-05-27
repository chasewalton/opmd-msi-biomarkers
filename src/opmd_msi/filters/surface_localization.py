"""Filter to retain only genes annotated as cell-surface localized."""
from __future__ import annotations

import pandas as pd


def filter_surface_localized(
    genes: list[str],
    hpa_annotations: pd.DataFrame,
    required_locations: tuple[str, ...] = ("Plasma membrane", "Cell surface"),
    required_secretome: tuple[str, ...] = ("Membrane proteins",),
    allow_secreted: bool = False,
    allow_predicted_only: bool = False,
) -> tuple[list[str], pd.DataFrame]:
    """Return (kept_genes, audit_df).

    The audit_df records the surface-status decision for every input gene so the
    reviewer can trace why each gene was kept or rejected.
    """
    # Normalize column names — accept several common HPA exports.
    loc_col = next(
        (c for c in ("Subcellular main location", "Subcellular location") if c in hpa_annotations.columns),
        None,
    )
    secr_col = next(
        (c for c in ("Secretome location", "Secretome category") if c in hpa_annotations.columns),
        None,
    )
    gene_col = next(
        (c for c in ("Gene", "gene", "Gene name") if c in hpa_annotations.columns),
        None,
    )
    if gene_col is None or loc_col is None:
        raise ValueError("HPA annotations missing required columns (Gene + Subcellular location)")

    indexed = hpa_annotations.set_index(gene_col)
    rows = []
    kept: list[str] = []
    for g in genes:
        if g not in indexed.index:
            rows.append({"gene": g, "surface_kept": False,
                         "reason": "missing from HPA"})
            continue
        loc = str(indexed.at[g, loc_col]) if loc_col else ""
        secr = str(indexed.at[g, secr_col]) if secr_col and secr_col in indexed.columns else ""

        loc_ok = any(rl in loc for rl in required_locations)
        secr_ok = (not required_secretome) or any(s in secr for s in required_secretome)

        if loc_ok and (secr_ok or allow_secreted):
            kept.append(g)
            rows.append({"gene": g, "surface_kept": True,
                         "subcellular": loc, "secretome": secr,
                         "reason": "annotated plasma membrane"})
        else:
            rows.append({"gene": g, "surface_kept": False,
                         "subcellular": loc, "secretome": secr,
                         "reason": "not annotated cell-surface"})
    audit = pd.DataFrame(rows)
    return kept, audit
