# opmd-msi-biomarkers

Reproducible ML pipeline that selects cell-surface biomarkers for molecular
surgical imaging (MSI) of oral potentially malignant disorders (OPMDs).

The real-data pipeline combines transcriptomic contrasts with reviewed epithelial
protein-outcome evidence, available humanized antibodies, and normal oral protein
background. Candidate eligibility and ranking are configured in `configs/ent.yaml`.
The synthetic demo separately exercises the original model-based workflow.

The reproducible synthetic-data demo runs end-to-end without network access
and recovers **PD-L1 (CD274), MET, and NOTCH1** as the top-3 candidates in its
generated data. The real-data ranking is labeled `curated_clinical_prioritization`.

## Original workflow documentation

The documents below describe the original workflow. For the implemented ENT
selection rules, use `configs/ent.yaml`, `manifests/`, and the real-data commands
below.

| File | What it covers |
| ---- | -------------- |
| [`REVIEWER_RESPONSE.md`](REVIEWER_RESPONSE.md) | Direct answers to the reviewer's questions. |
| [`METHODS.md`](METHODS.md) | Full methods write-up with dataset table and selection-criteria flowchart. |
| [`manifests/datasets.yaml`](manifests/datasets.yaml) | Authoritative list of every public dataset used (GEO + TCGA + reference panels). |
| [`configs/contrasts.yaml`](configs/contrasts.yaml) | Definitions of the three harmonized contrasts (C1, C2, C3). |
| [`configs/filters.yaml`](configs/filters.yaml) | Selection-rule thresholds and weights. |
| [`docs/selection_criteria.md`](docs/selection_criteria.md) | One-page selection flowchart. |
| [`tests/test_pipeline_synthetic.py`](tests/test_pipeline_synthetic.py) | Reproducibility anchor: top-3 = {CD274, MET, NOTCH1}. |

## Quickstart (no network, runs in seconds)

```bash
pip install -e ".[dev]"
make demo
```

Output goes to `outputs/demo/`. The CLI prints the ranked top-10 to the terminal
and writes per-stage audit files (`filter_surface_localization_audit.csv`,
`filter_baseline_normal_audit.csv`, `filter_contrast_strength_audit.csv`,
`contrast_*_pooled.csv`, `ranking.csv`).

## Real-data path

ENT selection uses reviewed epithelial protein-outcome evidence, available
humanized antibodies, and normal oral protein annotations in `manifests/`.
Filters and sorting weights are configured in `configs/ent.yaml`.

Eligibility requires supported epithelial localization, published human OPMD
protein-outcome evidence, and a supplier-listed humanized antibody. Candidates
must also meet either the progression contrast threshold or the positive
lesion-versus-normal contrast threshold. The score weights are 40% absolute
progression effect, 40% positive lesion effect, and 20% normal oral protein
background. Candidates with unreviewed evidence remain in the eligibility audit.

For an existing harmonized run with per-cohort contrast files:

```bash
opmd-msi prioritize-ent --run outputs/real --out outputs/real_ent
```

For parsed cohort expression and metadata parquets under `data/raw/`:

```bash
opmd-msi real --out outputs/real_ent
```

Both routes write `ent_filter_audit.csv`, `ranking.csv`, and
`ranking_top_promoted.csv`. Model coefficients, SHAP summaries, and correlations
are auxiliary outputs and do not determine the ENT composite score.

The checked local input run contains seven GEO series and the TCGA oral subset
(884 samples). GSE188981 is not part of that run. The original dataset manifest
contains earlier proposed cohorts and counts; it is not an inventory of the
saved run. Raw downloads must be parsed before `real` can load them; the GEO
and TCGA parser functions are in `src/opmd_msi/data/real_geo.py` and
`src/opmd_msi/data/real_tcga.py`.

`manifests/tracer_suitability.yaml` records the EGFR salivary-background concern
and the pending integrin citation as annotations, not exclusion rules. This
repository currently implements biomarker selection, not an image-quality
subanalysis.

## Layout

```
opmd-msi-biomarkers/
├── manifests/             dataset & cohort registries (single source of truth)
├── configs/               contrast definitions, filter thresholds, pipeline config
├── src/opmd_msi/          library code
│   ├── data/              fetchers (GEO/TCGA/HPA/Bgee/Atlas), harmonizer, synthetic
│   ├── contrasts/         effect sizes + random-effects meta-analysis
│   ├── filters/           surface localization, baseline-normal, contrast strength
│   ├── models/            L1 Cox, L1 logistic, leave-one-cohort-out CV
│   ├── importance/        SHAP ranking + composite score
│   ├── pipeline.py        end-to-end orchestrator
│   └── cli.py             `opmd-msi` command
├── scripts/01..06_*.py    numbered stage scripts (reviewers can read in order)
├── tests/                 pytest suite incl. end-to-end synthetic recovery
└── docs/                  selection-criteria flowchart, dataset table, pipeline diagram
```

## Testing

```bash
make test           # unit + integration
make test-cov       # with coverage
```

`tests/test_pipeline_synthetic.py` runs the full pipeline and asserts that the
top-3 candidates are PD-L1 (CD274), MET, and NOTCH1.

## License

Apache-2.0.
