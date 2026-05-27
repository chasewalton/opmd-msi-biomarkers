"""Stage 06 — compute SHAP importance + composite score → final ranking."""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
import yaml

from opmd_msi.importance.composite_score import compose_ranking
from opmd_msi.importance.shap_ranking import (
    compute_correlation_with_risk,
    compute_mean_abs_shap,
)
from opmd_msi.models.validation import leave_one_cohort_out


@click.command()
@click.option("--filtered", default="outputs/filtered", show_default=True)
@click.option("--harmonized", default="data/harmonized", show_default=True)
@click.option("--contrasts", default="outputs/contrasts", show_default=True)
@click.option("--config", default="configs/default.yaml", show_default=True)
@click.option("--out", default="outputs/ranking", show_default=True)
def main(filtered: str, harmonized: str, contrasts: str,
         config: str, out: str) -> None:
    fdir = Path(filtered)
    hdir = Path(harmonized)
    cdir = Path(contrasts)
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path(config).read_text())

    candidates = pd.read_parquet(fdir / "candidates.parquet")["gene"].tolist()
    expr = pd.read_parquet(hdir / "expression.harmonized.parquet")
    meta = pd.read_parquet(hdir / "metadata.harmonized.parquet")
    pooled_c3 = pd.read_parquet(cdir / "C3_pooled.parquet")
    pooled_c1 = pd.read_parquet(cdir / "C1_pooled.parquet")
    baseline = pd.read_parquet(fdir / "baseline_audit.parquet")

    loco = leave_one_cohort_out(expr, meta, candidates, pooled_c1=pooled_c1)
    shap_imp = compute_mean_abs_shap(loco)
    corr = compute_correlation_with_risk(expr, meta, candidates)

    rcfg = cfg.get("ranking", {}).get("composite_score", {})
    ranking = compose_ranking(
        pooled_c3=pooled_c3,
        shap_importance=shap_imp,
        corr_with_risk=corr,
        baseline_penalty=baseline,
        surface_kept=candidates,
        w_effect=rcfg.get("w_effect", 0.40),
        w_shap=rcfg.get("w_shap", 0.40),
        w_corr=rcfg.get("w_corr", 0.20),
        top_n=cfg.get("ranking", {}).get("top_n", 25),
    )
    ranking.to_csv(out_dir / "ranking.csv", index=False)
    ranking.head(cfg.get("ranking", {}).get("promoted_top_n", 3)) \
           .to_csv(out_dir / "ranking_top_promoted.csv", index=False)
    print(ranking.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
