"""Stage 02 — cross-cohort harmonization.

Reads per-cohort raw parquet files written by stage 01 (or synthetic outputs),
applies probe-to-symbol collapse + log transform + quantile normalization +
ComBat-style batch correction, and writes harmonized.parquet.
"""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from opmd_msi.data.harmonize import harmonize


@click.command()
@click.option("--raw", default="data/raw", show_default=True,
              help="Directory containing per-cohort .expression.parquet + .metadata.parquet")
@click.option("--out", default="data/harmonized", show_default=True)
def main(raw: str, out: str) -> None:
    raw_dir = Path(raw)
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)

    cohorts: dict[str, dict[str, pd.DataFrame]] = {}
    for ex_file in sorted(raw_dir.rglob("*.expression.parquet")):
        cid = ex_file.stem.replace(".expression", "")
        meta_file = ex_file.with_name(f"{cid}.metadata.parquet")
        if not meta_file.exists():
            continue
        cohorts[cid] = {
            "expression": pd.read_parquet(ex_file),
            "metadata": pd.read_parquet(meta_file),
        }
    if not cohorts:
        raise SystemExit(f"No cohorts found under {raw_dir}")
    res = harmonize(cohorts)
    res.expression.to_parquet(out_dir / "expression.harmonized.parquet")
    res.metadata.to_parquet(out_dir / "metadata.harmonized.parquet")
    print(f"Harmonized {len(cohorts)} cohorts -> "
          f"{res.expression.shape[0]} samples x {res.expression.shape[1]} genes")


if __name__ == "__main__":
    main()
