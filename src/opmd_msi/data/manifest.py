"""Load and validate the dataset manifest."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class DatasetEntry:
    accession: str
    source: str
    platform: str | None
    tissue: str | None
    n_samples: int
    contrasts: tuple[str, ...]
    tier: str
    longitudinal: bool = False
    has_outcome: bool = False
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Manifest:
    geo: tuple[DatasetEntry, ...]
    tcga: tuple[DatasetEntry, ...]
    reference_panels: tuple[dict[str, Any], ...]
    counts: dict[str, Any]

    @property
    def all_cohorts(self) -> tuple[DatasetEntry, ...]:
        return self.geo + self.tcga

    def cohorts_for_contrast(self, contrast_id: str) -> tuple[DatasetEntry, ...]:
        return tuple(c for c in self.all_cohorts if contrast_id in c.contrasts)


def _to_entry(d: dict[str, Any]) -> DatasetEntry:
    known = {"accession", "source", "platform", "tissue", "n_samples", "contrasts",
             "tier", "longitudinal", "has_outcome"}
    extra = {k: v for k, v in d.items() if k not in known}
    return DatasetEntry(
        accession=d["accession"],
        source=d["source"],
        platform=d.get("platform"),
        tissue=d.get("tissue"),
        n_samples=int(d.get("n_samples") or d.get("n_used") or 0),
        contrasts=tuple(d.get("contrasts", ())),
        tier=d.get("tier", "primary"),
        longitudinal=bool(d.get("longitudinal", False)),
        has_outcome=bool(d.get("has_outcome", False)),
        extra=extra,
    )


def load_manifest(path: str | Path) -> Manifest:
    data = yaml.safe_load(Path(path).read_text())
    geo = tuple(_to_entry(d) for d in data.get("geo_datasets", []))
    tcga = tuple(_to_entry(d) for d in data.get("tcga_datasets", []))
    refs = tuple(data.get("reference_panels", []))
    counts = data.get("counts", {})
    return Manifest(geo=geo, tcga=tcga, reference_panels=refs, counts=counts)
