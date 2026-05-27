"""EBI Expression Atlas fetcher + synthetic fixture."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

ATLAS_BASELINE_API = "https://www.ebi.ac.uk/gxa/json/baseline"


def fetch_atlas_experiment(experiment_id: str, out_dir: str | Path) -> dict[str, Any]:
    """Download baseline expression for one Expression Atlas experiment."""
    try:
        import requests
    except ImportError as e:  # pragma: no cover
        raise ImportError("requests is required for Atlas fetching") from e
    out = Path(out_dir) / "Atlas"
    out.mkdir(parents=True, exist_ok=True)
    url = f"{ATLAS_BASELINE_API}/{experiment_id}"
    logger.info("Downloading Atlas experiment %s", experiment_id)
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    p = out / f"{experiment_id}.json"
    p.write_text(resp.text)
    return {"experiment_id": experiment_id, "path": str(p)}


def synthetic_atlas_annotations() -> pd.DataFrame:
    """Atlas-like baseline table for the demo gene universe."""
    from opmd_msi.data.synthetic import ALL_GENES, HIGH_BASELINE_NORMAL
    rows = []
    for gene in ALL_GENES:
        rows.append({
            "Gene": gene,
            "median_tpm_across_E-MTAB-2836": 16.5 if gene in HIGH_BASELINE_NORMAL else 2.7,
            "n_tissues_above_threshold": 5 if gene in HIGH_BASELINE_NORMAL else 0,
        })
    return pd.DataFrame(rows)
