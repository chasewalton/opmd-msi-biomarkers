"""Stage 01 — fetch real public data per manifests/datasets.yaml.

Requires the `[fetch]` extras: `pip install -e '.[fetch]'`.

Usage:
    python scripts/01_fetch_data.py --manifest manifests/datasets.yaml --out data/raw
"""
from __future__ import annotations

import logging
from pathlib import Path

import click

from opmd_msi.data.fetchers.atlas import fetch_atlas_experiment
from opmd_msi.data.fetchers.bgee import fetch_bgee
from opmd_msi.data.fetchers.geo import fetch_geo
from opmd_msi.data.fetchers.hpa import fetch_hpa
from opmd_msi.data.fetchers.tcga import fetch_tcga_hnsc_oral
from opmd_msi.data.manifest import load_manifest

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(name)s :: %(message)s")
log = logging.getLogger(__name__)


@click.command()
@click.option("--manifest", default="manifests/datasets.yaml", show_default=True)
@click.option("--out", default="data/raw", show_default=True,
              help="Output root directory for raw artifacts.")
@click.option("--skip", multiple=True, help="Accessions to skip.")
def main(manifest: str, out: str, skip: tuple[str, ...]) -> None:
    mf = load_manifest(manifest)
    out_root = Path(out)
    out_root.mkdir(parents=True, exist_ok=True)
    skip_set = set(skip)

    log.info("Fetching %d GEO + %d TCGA datasets + %d reference panels",
             len(mf.geo), len(mf.tcga), len(mf.reference_panels))

    for entry in mf.geo:
        if entry.accession in skip_set:
            log.info("Skipping %s", entry.accession)
            continue
        try:
            fetch_geo(entry.accession, out_root / "GEO")
        except Exception as e:
            log.error("Failed to fetch %s: %s", entry.accession, e)

    for entry in mf.tcga:
        if entry.accession in skip_set:
            continue
        try:
            fetch_tcga_hnsc_oral(out_root / "TCGA")
        except Exception as e:
            log.error("Failed to fetch %s: %s", entry.accession, e)

    # Reference panels
    try:
        fetch_hpa(out_root / "reference")
    except Exception as e:
        log.error("HPA fetch failed: %s", e)
    try:
        fetch_bgee(out_root / "reference")
    except Exception as e:
        log.error("Bgee fetch failed: %s", e)
    for exp in ("E-MTAB-2836", "E-MTAB-3358"):
        try:
            fetch_atlas_experiment(exp, out_root / "reference")
        except Exception as e:
            log.error("Atlas %s fetch failed: %s", exp, e)


if __name__ == "__main__":
    main()
