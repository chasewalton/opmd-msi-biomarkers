"""Bgee fetcher + synthetic fixture (normal-tissue expression baseline)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd

logger = logging.getLogger(__name__)

BGEE_HUMAN_RNA_TSV = (
    "https://www.bgee.org/ftp/current/download/calls/expr_calls/"
    "Homo_sapiens_expr_simple.tsv.zip"
)


def fetch_bgee(out_dir: str | Path) -> dict[str, Any]:
    try:
        import requests
    except ImportError as e:  # pragma: no cover
        raise ImportError("requests is required for Bgee fetching") from e
    out = Path(out_dir) / "Bgee"
    out.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading Bgee simple expression calls…")
    resp = requests.get(BGEE_HUMAN_RNA_TSV, timeout=600)
    resp.raise_for_status()
    raw_path = out / "Homo_sapiens_expr_simple.tsv.zip"
    raw_path.write_bytes(resp.content)
    return {"raw_path": str(raw_path)}


def synthetic_bgee_annotations() -> pd.DataFrame:
    """Bgee-like baseline call table for the demo gene universe.

    Returns the per-gene maximum normal-tissue expression among non-oral
    tissues — i.e. the value the baseline-normal penalty filter consumes.
    """
    from opmd_msi.data.synthetic import ALL_GENES, HIGH_BASELINE_NORMAL
    rows = []
    for gene in ALL_GENES:
        rows.append({
            "Gene": gene,
            "max_non_oral_normal_tpm": 18.0 if gene in HIGH_BASELINE_NORMAL else 2.5,
            "n_tissues_with_high_expression": 6 if gene in HIGH_BASELINE_NORMAL else 1,
        })
    return pd.DataFrame(rows)
