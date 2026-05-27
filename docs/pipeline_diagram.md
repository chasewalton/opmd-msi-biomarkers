# Pipeline diagram

```
                           ┌──────────────────────────┐
                           │ manifests/datasets.yaml  │
                           └────────────┬─────────────┘
                                        │
                       ┌────────────────┴───────────────┐
                       │                                │
                       ▼                                ▼
       ┌──────────────────────────┐      ┌──────────────────────────┐
 01    │ scripts/01_fetch_data.py │      │ data/synthetic.py        │
       │  GEO + TCGA + HPA + Bgee │      │  (deterministic demo)    │
       │  + Expression Atlas      │      └────────────┬─────────────┘
       └────────────┬─────────────┘                   │
                    │                                 │
                    └──────────────┬──────────────────┘
                                   ▼
                ┌────────────────────────────────────┐
 02             │ scripts/02_harmonize.py            │
                │  HGNC collapse, log2, quantile     │
                │  normalization, ComBat batch corr. │
                └────────────────┬───────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐
 03             │ scripts/03_contrasts.py            │
                │  C1, C2, C3 per-cohort effects     │
                │  + DerSimonian-Laird REML pooling  │
                └────────────────┬───────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐
 04             │ scripts/04_filter.py               │
                │  surface localization (HPA),       │
                │  baseline-normal penalty (3 panel),│
                │  contrast-strength gate            │
                └────────────────┬───────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐
 05             │ scripts/05_train_models.py         │
                │  L1 Cox + L1 logistic,             │
                │  leave-one-cohort-out CV,          │
                │  monotonic sign constraint from C1 │
                └────────────────┬───────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐
 06             │ scripts/06_rank.py                 │
                │  mean |SHAP| over folds,           │
                │  |corr with risk|,                 │
                │  composite score → ranking         │
                └────────────────┬───────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐
                │   outputs/ranking/ranking.csv      │
                │   top-3 promoted → wet-lab          │
                └────────────────────────────────────┘
```

The orchestrator class `src/opmd_msi/pipeline.py::run_pipeline_synthetic`
runs stages 02–06 in-process for the demo path; the same code paths are
invoked by the numbered stage scripts for the real-data path.
