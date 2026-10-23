"""Real-data annotation parsers: HPA surface localization + normal baseline,
UniProt subcellular-location fallback for canonical transmembrane proteins.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)


def parse_hpa(proteinatlas_tsv: str | Path) -> pd.DataFrame:
    """Parse the HPA proteinatlas.tsv into a tidy annotation table.

    Returns columns: gene, subcellular_locations, secretome_location,
    tissue_RNA_max_non_oral_nTPM, hpa_passes_surface_filter.
    """
    df = pd.read_csv(proteinatlas_tsv, sep="\t", low_memory=False)
    # HPA can have duplicate rows per gene if multiple Ensembl IDs map to same symbol;
    # collapse by taking the first row per gene.
    df = df.drop_duplicates(subset=["Gene"], keep="first").reset_index(drop=True)
    out = pd.DataFrame({"gene": df["Gene"].astype(str)})
    out["subcellular_locations"] = df["Subcellular location"].fillna("").astype(str)
    out["secretome_location"] = df["Secretome location"].fillna("").astype(str)

    # Per-tissue RNA columns we care about
    oral_tissue_cols = [c for c in df.columns
                        if c.startswith("Tissue RNA")
                        and any(t in c.lower()
                                for t in ["oral mucosa", "tongue", "skin", "esophagus"])]
    non_oral_cols = [c for c in df.columns
                     if c.startswith("Tissue RNA")
                     and c not in oral_tissue_cols]
    # numeric coerce
    for c in oral_tissue_cols + non_oral_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    out["tissue_RNA_max_non_oral_nTPM"] = (
        df[non_oral_cols].max(axis=1) if non_oral_cols else 0.0
    )
    out["tissue_RNA_max_oral_nTPM"] = (
        df[oral_tissue_cols].max(axis=1) if oral_tissue_cols else 0.0
    )
    # HPA surface annotation: "Plasma membrane" or "Cell surface" in subcellular
    out["hpa_passes_surface_filter"] = out["subcellular_locations"].str.contains(
        r"Plasma membrane|Cell surface", regex=True, na=False
    )
    return out


# ---------------------------------------------------------------------------
# UniProt fallback for canonical transmembrane proteins
# ---------------------------------------------------------------------------

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# Hand-curated set of canonical transmembrane proteins relevant to OPMD/HNSCC.
# These are well-established type-I/II single-pass or multi-pass membrane
# receptors that HPA may annotate ambiguously (e.g. NOTCH1's intracellular
# domain localises to the nucleus and dominates HPA's antibody signal).
# Source: UniProt subcellular location annotations as of release 2024_01.
CANONICAL_TRANSMEMBRANE: dict[str, str] = {
    "NOTCH1": "Type I single-pass membrane protein (UniProt P46531)",
    "NOTCH2": "Type I single-pass membrane protein (UniProt Q04721)",
    "NOTCH3": "Type I single-pass membrane protein (UniProt Q9UM47)",
    "NOTCH4": "Type I single-pass membrane protein (UniProt Q99466)",
}


def uniprot_membrane_lookup(genes: list[str]) -> dict[str, bool]:
    """Quick UniProt query: returns {gene: True if subcellular location mentions membrane}.

    For the OPMD pipeline we keep this lightweight — we hard-code the known
    transmembrane receptors that HPA mislocalises (above), and only query
    UniProt for genes not in either HPA's plasma-membrane set or the curated
    canonical set. This keeps the network footprint small.
    """
    out: dict[str, bool] = {g: True for g in genes if g in CANONICAL_TRANSMEMBRANE}
    to_query = [g for g in genes if g not in out]
    if not to_query:
        return out
    # Batch query in chunks of 25 with reviewed=true filter
    for i in range(0, len(to_query), 25):
        batch = to_query[i:i+25]
        q = " OR ".join(f"gene:{g}" for g in batch)
        params = {
            "query": f"({q}) AND organism_id:9606 AND reviewed:true",
            "fields": "gene_names,cc_subcellular_location",
            "format": "tsv", "size": 100,
        }
        try:
            r = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
            r.raise_for_status()
            df = pd.read_csv(pd.io.common.StringIO(r.text), sep="\t")
        except Exception as e:
            log.warning("uniprot fetch failed for batch %d: %s", i, e)
            for g in batch:
                out[g] = False
            continue
        for g in batch:
            # find row whose gene_names contains g
            mask = df["Gene Names"].fillna("").str.contains(rf"\b{g}\b", regex=True, na=False)
            if not mask.any():
                out[g] = False
                continue
            loc = " ".join(df.loc[mask, "Subcellular location [CC]"].fillna("").astype(str)).lower()
            out[g] = ("cell membrane" in loc or "plasma membrane" in loc
                      or "single-pass" in loc or "multi-pass" in loc)
    return out


def combined_surface_filter(
    genes: list[str],
    hpa_annotations: pd.DataFrame,
    use_uniprot_fallback: bool = True,
) -> tuple[list[str], pd.DataFrame]:
    """Apply HPA → UniProt fallback surface filter.

    Returns (kept_genes, audit_df).
    """
    hpa = hpa_annotations.set_index("gene")
    rows = []
    pending_uniprot = []
    for g in genes:
        if g not in hpa.index:
            rows.append({"gene": g, "hpa_membrane": False,
                         "uniprot_membrane": None, "kept": False,
                         "reason": "missing from HPA"})
            continue
        hpa_pass = bool(hpa.at[g, "hpa_passes_surface_filter"])
        if hpa_pass:
            rows.append({"gene": g, "hpa_membrane": True,
                         "uniprot_membrane": None, "kept": True,
                         "reason": "HPA Plasma membrane / Cell surface"})
        else:
            pending_uniprot.append(g)

    if use_uniprot_fallback and pending_uniprot:
        uni = uniprot_membrane_lookup(pending_uniprot)
        for g in pending_uniprot:
            up = uni.get(g, False)
            rows.append({
                "gene": g, "hpa_membrane": False,
                "uniprot_membrane": up,
                "kept": up,
                "reason": ("UniProt membrane fallback"
                           if up else "Not membrane in HPA or UniProt"),
            })
    elif pending_uniprot:
        for g in pending_uniprot:
            rows.append({"gene": g, "hpa_membrane": False,
                         "uniprot_membrane": None, "kept": False,
                         "reason": "Not in HPA membrane (UniProt fallback disabled)"})

    audit = pd.DataFrame(rows).drop_duplicates(subset=["gene"], keep="first")
    kept = audit.loc[audit["kept"], "gene"].tolist()
    return kept, audit
