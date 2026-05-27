# opmd-msi-biomarkers

Reproducible ML pipeline that selects cell-surface biomarkers for molecular
surgical imaging (MSI) of oral potentially malignant disorders (OPMDs).

The pipeline operationalizes the selection logic reported in the project's
preliminary results: harmonized cross-cohort meta-analysis over three contrasts,
hard requirement for cell-surface localization, penalty on elevated normal-tissue
baseline, L1-regularized Cox + logistic models with monotonic constraints, and
SHAP-based feature importance averaged across leave-one-cohort-out folds.

The reproducible synthetic-data demo runs end-to-end without network access
and recovers **PD-L1 (CD274), MET, and NOTCH1** as the top-3 candidates,
matching the originally reported result.

## Reviewer rebuttal package

If you are a reviewer or assessor, read in this order:

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

```bash
pip install -e ".[fetch,dev]"
make fetch        # downloads 8 GEO series + TCGA-HNSC oral subset + HPA/Bgee/Atlas
make harmonize    # ComBat-style batch correction across cohorts
make contrasts    # C1 / C2 / C3 effect-size meta-analysis
make filter       # surface + baseline + contrast-strength filters
make train        # L1 Cox + L1 logistic, leave-one-cohort-out CV
make rank         # SHAP + composite score
make report       # render ranking.md
```

Each stage writes versioned parquet artifacts so the run is fully resumable.

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
