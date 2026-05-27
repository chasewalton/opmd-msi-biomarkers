# Selection-criteria flowchart

Single-page summary of how a candidate gene is promoted from "any expressed
gene in the harmonized matrix" to "top-3 MSI biomarker."

```
┌─────────────────────────────────────────────────────────────────────┐
│  Harmonized matrix: 1049 samples × ~17000 HGNC genes                │
│  (8 GEO + 1 TCGA-HNSC oral subset, ComBat batch-corrected)          │
└───────────────────────────────┬─────────────────────────────────────┘
                                │
                                ▼
                ┌───────────────────────────────┐
                │ Filter 1: cell-surface (HPA)  │  → drop intracellular,
                │  - Plasma membrane / surface  │     nuclear, secreted,
                │  - "Membrane proteins"        │     predicted-only.
                │  - Experimental evidence      │
                └───────────────┬───────────────┘
                                │
                                ▼
                ┌───────────────────────────────┐
                │ Filter 2: low normal baseline │  → drop if ≥2/3 panels
                │  - HPA, Bgee, Atlas (≥2/3)    │     exceed 50 nTPM;
                │  - Excl. oral-relevant        │     penalize if >15 nTPM.
                │    tissues                    │
                └───────────────┬───────────────┘
                                │
                                ▼
                ┌───────────────────────────────┐
                │ Filter 3: pooled C3 contrast  │  → |Hedges g| ≥ 0.50,
                │  with C1 + C2 direction       │     FDR ≤ 0.05,
                │  agreement                    │     ≥4 cohorts same dir.,
                │                               │     same sign in C1 & C2.
                └───────────────┬───────────────┘
                                │
                                ▼
                ┌───────────────────────────────┐
                │ Models: L1 Cox + L1 logistic  │  → monotonic sign
                │  with leave-one-cohort-out CV │     constraint from C1.
                └───────────────┬───────────────┘
                                │
                                ▼
                ┌───────────────────────────────┐
                │ Score: composite of           │
                │  • |pooled C3 effect|         │  → top-N by composite,
                │  • mean |SHAP| over folds     │     promoted top-3 to
                │  • |corr with risk label|     │     wet-lab validation.
                │  × baseline penalty mult.     │
                └───────────────────────────────┘
                                │
                                ▼
                       PD-L1, MET, NOTCH1
```

Numerical thresholds live in [`configs/filters.yaml`](../configs/filters.yaml);
scoring weights in [`configs/default.yaml`](../configs/default.yaml).
