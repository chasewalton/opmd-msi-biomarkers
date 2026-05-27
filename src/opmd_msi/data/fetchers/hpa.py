"""Human Protein Atlas fetcher + synthetic fixture.

The pipeline uses HPA for two purposes:
  1. Cell-surface localization filter (subcellular location annotation).
  2. Normal-tissue baseline expression (per-tissue nTPM) as one of three
     reference panels combined in the baseline-normal penalty.

The synthetic fixture below provides HPA-compatible annotations for the genes
used in the demo so the surface-localization filter behaves identically to its
real-data path.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

HPA_SEARCH_DOWNLOAD = (
    "https://www.proteinatlas.org/api/search_download.php"
    "?columns=eg,g,scl,scml,scal,t_RNA_oral_mucosa,t_RNA_tongue,"
    "t_RNA_skin,t_RNA_esophagus,t_RNA_breast,t_RNA_lung,"
    "secretome_location&format=tsv&compress=no"
)


def fetch_hpa(out_dir: str | Path) -> dict[str, Any]:
    """Download HPA annotations covering subcellular location + normal RNA."""
    try:
        import requests
    except ImportError as e:  # pragma: no cover
        raise ImportError("requests is required for HPA fetching") from e

    out = Path(out_dir) / "HPA"
    out.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading HPA annotations…")
    resp = requests.get(HPA_SEARCH_DOWNLOAD, timeout=300)
    resp.raise_for_status()
    raw_path = out / "hpa_subcellular_normal_rna.tsv"
    raw_path.write_text(resp.text)
    df = pd.read_csv(raw_path, sep="\t")
    parquet_path = out / "hpa.parquet"
    df.to_parquet(parquet_path)
    return {"raw_path": str(raw_path), "parquet_path": str(parquet_path),
            "n_genes": len(df)}


# ---------------------------------------------------------------------------
# Synthetic fixture
# ---------------------------------------------------------------------------

def synthetic_hpa_annotations() -> pd.DataFrame:
    """Return an HPA-like annotation table for the demo gene universe.

    Columns mirror what the real `fetch_hpa` produces.
    """
    from opmd_msi.data.synthetic import (
        ALL_GENES,
        HIGH_BASELINE_NORMAL,
        SURFACE_GENES,
    )

    rows = []
    for gene in ALL_GENES:
        is_surface = gene in SURFACE_GENES
        is_high_baseline = gene in HIGH_BASELINE_NORMAL
        rows.append({
            "Gene": gene,
            "Ensembl": f"ENSG_SYN_{gene}",
            "Subcellular location": "Plasma membrane" if is_surface else "Cytoplasm",
            "Subcellular main location": "Plasma membrane" if is_surface else "Cytoplasm",
            "Secretome location": "Membrane proteins" if is_surface else "Intracellular",
            # nTPM in normal oral-relevant tissues
            "Tissue RNA - oral mucosa [nTPM]": 25.0 if is_high_baseline else 3.0,
            "Tissue RNA - tongue [nTPM]":       18.0 if is_high_baseline else 2.5,
            "Tissue RNA - skin [nTPM]":         22.0 if is_high_baseline else 3.0,
            "Tissue RNA - esophagus [nTPM]":    20.0 if is_high_baseline else 2.8,
            "Tissue RNA - breast [nTPM]":       4.0,
            "Tissue RNA - lung [nTPM]":         3.5,
        })
    return pd.DataFrame(rows)
