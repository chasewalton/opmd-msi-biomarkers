"""Stage 03 — compute the three harmonized contrasts and pooled meta-analysis."""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from opmd_msi.contrasts.contrasts import run_all_contrasts
from opmd_msi.data.harmonize import HarmonizedDataset


@click.command()
@click.option("--harmonized", default="data/harmonized", show_default=True)
@click.option("--contrasts-config", default="configs/contrasts.yaml", show_default=True)
@click.option("--out", default="outputs/contrasts", show_default=True)
def main(harmonized: str, contrasts_config: str, out: str) -> None:
    hp = Path(harmonized)
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    expr = pd.read_parquet(hp / "expression.harmonized.parquet")
    meta = pd.read_parquet(hp / "metadata.harmonized.parquet")
    h = HarmonizedDataset(
        expression=expr, metadata=meta,
        cohort_ids=tuple(sorted(meta["cohort_id"].unique())),
        gene_ids=tuple(expr.columns),
    )
    results = run_all_contrasts(h, contrasts_config)
    for cid, cres in results.items():
        cres.per_cohort.to_parquet(out_dir / f"{cid}_per_cohort.parquet")
        cres.pooled.to_parquet(out_dir / f"{cid}_pooled.parquet")
        print(f"{cid}: per-cohort rows={len(cres.per_cohort)}, "
              f"pooled genes={len(cres.pooled)}")


if __name__ == "__main__":
    main()
