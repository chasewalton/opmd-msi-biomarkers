"""Stage 05 — fit L1 Cox + L1 logistic under leave-one-cohort-out CV."""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd

from opmd_msi.models.validation import leave_one_cohort_out


@click.command()
@click.option("--filtered", default="outputs/filtered", show_default=True)
@click.option("--harmonized", default="data/harmonized", show_default=True)
@click.option("--contrasts", default="outputs/contrasts", show_default=True)
@click.option("--out", default="outputs/models", show_default=True)
def main(filtered: str, harmonized: str, contrasts: str, out: str) -> None:
    fdir = Path(filtered)
    hdir = Path(harmonized)
    cdir = Path(contrasts)
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_parquet(fdir / "candidates.parquet")["gene"].tolist()
    expr = pd.read_parquet(hdir / "expression.harmonized.parquet")
    meta = pd.read_parquet(hdir / "metadata.harmonized.parquet")
    pooled_c1 = pd.read_parquet(cdir / "C1_pooled.parquet")

    loco = leave_one_cohort_out(expr, meta, candidates, pooled_c1=pooled_c1)
    loco.coef_matrix_cox.to_parquet(out_dir / "loco_coefs_cox.parquet")
    loco.coef_matrix_logistic.to_parquet(out_dir / "loco_coefs_logistic.parquet")
    loco.average_abs_coefs().to_parquet(out_dir / "mean_abs_coefs.parquet")
    print(f"Trained {len(loco.fold_fits_cox)} Cox folds and "
          f"{len(loco.fold_fits_logistic)} logistic folds.")


if __name__ == "__main__":
    main()
