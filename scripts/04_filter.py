"""Stage 04 — apply selection filters: surface localization, baseline-normal, contrast strength."""
from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
import yaml

from opmd_msi.data.fetchers.atlas import synthetic_atlas_annotations
from opmd_msi.data.fetchers.bgee import synthetic_bgee_annotations
from opmd_msi.data.fetchers.hpa import synthetic_hpa_annotations
from opmd_msi.filters.baseline_normal import compute_baseline_penalty
from opmd_msi.filters.contrast_threshold import filter_by_contrast_strength
from opmd_msi.filters.surface_localization import filter_surface_localized


@click.command()
@click.option("--contrasts", default="outputs/contrasts", show_default=True)
@click.option("--harmonized", default="data/harmonized", show_default=True)
@click.option("--filters-config", default="configs/filters.yaml", show_default=True)
@click.option("--out", default="outputs/filtered", show_default=True)
@click.option("--use-synthetic-references", is_flag=True, default=False,
              help="Use synthetic HPA/Bgee/Atlas fixtures (for the demo path).")
def main(contrasts: str, harmonized: str, filters_config: str,
         out: str, use_synthetic_references: bool) -> None:
    cdir = Path(contrasts)
    hdir = Path(harmonized)
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    fcfg = yaml.safe_load(Path(filters_config).read_text())

    pooled_c1 = pd.read_parquet(cdir / "C1_pooled.parquet")
    pooled_c2 = pd.read_parquet(cdir / "C2_pooled.parquet")
    pooled_c3 = pd.read_parquet(cdir / "C3_pooled.parquet")
    expr = pd.read_parquet(hdir / "expression.harmonized.parquet")
    all_genes = list(expr.columns)

    if use_synthetic_references:
        hpa = synthetic_hpa_annotations()
        bgee = synthetic_bgee_annotations()
        atlas = synthetic_atlas_annotations()
    else:
        # Real-data path: load from data/raw/reference/
        hpa = pd.read_parquet("data/raw/reference/HPA/hpa.parquet")
        bgee = pd.read_parquet("data/raw/reference/Bgee/bgee.parquet")
        atlas = pd.read_parquet("data/raw/reference/Atlas/atlas.parquet")

    sfilt = fcfg["surface_localization"]
    surface_kept, surface_audit = filter_surface_localized(
        all_genes, hpa,
        required_locations=tuple(sfilt["required_subcellular_locations"]),
        required_secretome=tuple(sfilt.get("required_secretome_categories", ())),
    )
    bfilt = fcfg["baseline_normal_penalty"]
    baseline_kept, baseline_audit = compute_baseline_penalty(
        surface_kept, hpa, bgee, atlas,
        high_baseline_tpm=bfilt["high_baseline_tpm"],
        excluded_baseline_tpm=bfilt["excluded_baseline_tpm"],
        min_panels_in_agreement=bfilt["min_panels_in_agreement"],
        oral_relevant_tissues=tuple(bfilt.get("oral_relevant_tissues", ())),
    )
    cfilt = fcfg["lesion_to_normal_contrast"]
    contrast_kept, contrast_audit = filter_by_contrast_strength(
        pooled_c3, pooled_c1, pooled_c2,
        min_abs_pooled_estimate=cfilt["min_hedges_g"],
        max_p_fdr=cfilt["max_meta_p_fdr"],
        min_cohorts_same_direction=cfilt["min_cohorts_with_same_direction"],
        require_C1_consistent_direction=cfilt["require_C1_consistent_direction"],
        require_C2_consistent_direction=cfilt["require_C2_consistent_direction"],
    )
    candidates = sorted(set(baseline_kept) & set(contrast_kept))

    surface_audit.to_parquet(out_dir / "surface_audit.parquet")
    baseline_audit.to_parquet(out_dir / "baseline_audit.parquet")
    contrast_audit.to_parquet(out_dir / "contrast_audit.parquet")
    pd.DataFrame({"gene": candidates}).to_parquet(out_dir / "candidates.parquet")
    print(f"Surface kept: {len(surface_kept)}, "
          f"baseline kept: {len(baseline_kept)}, "
          f"contrast kept: {len(contrast_kept)}, "
          f"intersection: {len(candidates)}")


if __name__ == "__main__":
    main()
