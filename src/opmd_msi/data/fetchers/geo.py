"""GEO fetcher.

Downloads SOFT/series-matrix files via the GEOparse library (opt-in extra).
Parses expression and sample metadata; writes parquet plus an unmodified copy
of the raw matrix for provenance.

Activated by `pip install -e .[fetch]`.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

# Mapping from GEO sample-level annotation columns to our harmonized schema.
# Each cohort has bespoke column names in the source SOFT file. The mappings
# here are conservative; tests should be added as new cohorts are wired up.
GEO_METADATA_MAPS: dict[str, dict[str, str]] = {
    "GSE26549": {
        "histology": "tissue_class",
        "transformation_status": "transformed",
        "time_to_transformation_(months)": "time_to_transformation_months",
    },
    "GSE85195": {
        "diagnosis": "tissue_class",
        "progression": "transformed",
    },
    "GSE30784": {
        "tissue": "tissue_class",
        "patient": "patient_id",
    },
    "GSE37991": {"tissue": "tissue_class", "patient": "patient_id"},
    "GSE9844":  {"tissue": "tissue_class"},
    "GSE6791":  {"tissue": "tissue_class", "hpv_status": "hpv_status"},
    "GSE25099": {"tissue": "tissue_class"},
    "GSE188981": {
        "histology": "tissue_class",
        "outcome": "transformed",
        "time_to_event_months": "time_to_transformation_months",
    },
}


def fetch_geo(accession: str, out_dir: str | Path) -> dict[str, Any]:
    """Download a GEO series and return paths to written artifacts.

    Raises ImportError if GEOparse is not installed (install with `[fetch]`).
    """
    try:
        import GEOparse  # type: ignore[import-not-found]
    except ImportError as e:  # pragma: no cover - opt-in dependency
        raise ImportError(
            "GEOparse is required for real GEO fetching. "
            "Install with: pip install -e '.[fetch]'"
        ) from e

    out = Path(out_dir) / accession
    out.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading GEO %s into %s", accession, out)

    gse = GEOparse.get_GEO(geo=accession, destdir=str(out), silent=True)
    # Combine sample expression matrices
    sample_exprs: list[pd.Series] = []
    sample_meta: list[dict[str, Any]] = []
    for gsm_name, gsm in gse.gsms.items():
        expr = gsm.table.set_index("ID_REF")["VALUE"].astype(float)
        expr.name = gsm_name
        sample_exprs.append(expr)
        m: dict[str, Any] = {"sample_id": gsm_name, "cohort_id": accession}
        for k, v in gsm.metadata.items():
            m[k] = v[0] if isinstance(v, list) and len(v) == 1 else v
        sample_meta.append(m)
    expression = pd.concat(sample_exprs, axis=1).T  # samples x probes
    metadata = pd.DataFrame(sample_meta).set_index("sample_id")

    # Apply per-cohort column rename
    rename_map = GEO_METADATA_MAPS.get(accession, {})
    metadata = metadata.rename(columns=rename_map)

    expression_path = out / "expression.raw.parquet"
    metadata_path = out / "metadata.raw.parquet"
    expression.to_parquet(expression_path)
    metadata.to_parquet(metadata_path)
    return {
        "accession": accession,
        "expression_path": str(expression_path),
        "metadata_path": str(metadata_path),
        "n_samples": expression.shape[0],
        "n_probes": expression.shape[1],
    }
