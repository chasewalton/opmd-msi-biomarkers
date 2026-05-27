"""Command-line entry point: `opmd-msi`.

Subcommands:
  demo     Run the full synthetic-data pipeline and write outputs.
  report   Render the reviewer-facing ranking + per-marker reports.
"""
from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console
from rich.table import Table

from opmd_msi.pipeline import run_pipeline_synthetic

console = Console()


@click.group(help="OPMD MSI biomarker selection pipeline.")
def main() -> None:
    pass


@main.command(help="Run the full synthetic-data pipeline and write outputs.")
@click.option("--out", "out_dir", default="outputs/demo", show_default=True,
              help="Output directory.")
@click.option("--seed", default=42, show_default=True, type=int)
@click.option("--config", default="configs/default.yaml", show_default=True)
@click.option("--contrasts-config", default="configs/contrasts.yaml", show_default=True)
@click.option("--filters-config", default="configs/filters.yaml", show_default=True)
def demo(out_dir: str, seed: int, config: str,
         contrasts_config: str, filters_config: str) -> None:
    console.rule("[bold cyan]OPMD MSI synthetic pipeline")
    res = run_pipeline_synthetic(
        config_path=config,
        contrasts_config_path=contrasts_config,
        filters_config_path=filters_config,
        seed=seed,
        out_dir=out_dir,
    )
    console.print(f"\n[bold]Harmonized:[/] {res.harmonized.expression.shape[0]} samples "
                  f"x {res.harmonized.expression.shape[1]} genes "
                  f"across {len(res.harmonized.cohort_ids)} cohorts")
    console.print(f"[bold]Surface-localised survivors:[/] {len(res.surface_kept)}")
    console.print(f"[bold]After baseline-normal penalty:[/] {len(res.baseline_kept)}")
    console.print(f"[bold]After contrast-strength filter:[/] {len(res.candidate_genes)} candidate genes\n")

    tbl = Table(title="Final ranking (top 10)")
    tbl.add_column("Rank", justify="right")
    tbl.add_column("Gene")
    tbl.add_column("Composite", justify="right")
    tbl.add_column("|C3| pooled g", justify="right")
    tbl.add_column("Mean |SHAP|", justify="right")
    tbl.add_column("|corr|", justify="right")
    for _, row in res.ranking.head(10).iterrows():
        tbl.add_row(
            str(int(row["rank"])),
            str(row["gene"]),
            f"{row['composite_score']:.3f}",
            f"{row['abs_pooled_c3']:.2f}",
            f"{row['mean_abs_shap']:.3f}",
            f"{row['abs_corr']:.2f}",
        )
    console.print(tbl)
    console.print(f"\n[green]Wrote ranked outputs to[/] {Path(out_dir).resolve()}")


@main.command(help="Render the reviewer-facing ranking + per-marker reports.")
@click.option("--ranking", default="outputs/demo/ranking.csv", show_default=True)
@click.option("--out", "out_dir", default="outputs/report", show_default=True)
def report(ranking: str, out_dir: str) -> None:
    import pandas as pd
    df = pd.read_csv(ranking)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    md = ["# OPMD MSI biomarker ranking — report\n",
          f"Source: `{ranking}`\n",
          "## Top 25 candidates\n",
          df.to_markdown(index=False)]
    (out / "ranking.md").write_text("\n".join(md))
    console.print(f"[green]Wrote[/] {out / 'ranking.md'}")


if __name__ == "__main__":
    main()
