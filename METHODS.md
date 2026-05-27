# Methods — ML-driven selection of MSI biomarkers for OPMD

## Overview

We assembled longitudinal and cross-sectional transcriptomic data on oral
potentially malignant disorders (OPMDs) and oral squamous cell carcinoma
(OSCC) from the Gene Expression Omnibus (GEO) and The Cancer Genome Atlas
(TCGA-HNSC). Three independent reference panels (Human Protein Atlas, Bgee,
EBI Expression Atlas) were used to constrain candidate biomarkers to cell-
surface proteins with low baseline expression in normal tissues — both
prerequisites for selective molecular surgical imaging (MSI). Effect sizes
from three harmonized contrasts (progression, paired lesion–adjacent, lesion–
normal) were combined via random-effects meta-analysis, gated by hard
filters, and integrated with L1-regularized survival and classification
models. SHAP-based feature importance averaged across leave-one-cohort-out
cross-validation folds, together with effect-size magnitude and direct
correlation with transformation, produced the composite ranking from which
the top three candidates — PD-L1 (CD274), MET, and NOTCH1 — were promoted.

## Table 1 — Datasets used

| ID | Source | Platform | Tissue | n | Contrast(s) | Tier | Reference |
| -- | ------ | -------- | ------ | -: | ----------- | ---- | --------- |
| GSE26549 | GEO | Affymetrix HG-U133 Plus 2.0 | Oral leukoplakia (longitudinal) | 86 | C1 | primary | Saintigny *et al.*, Cancer Prev Res 2011 |
| GSE85195 | GEO | Illumina HumanHT-12 v4 | Oral premalignant lesions | 62 | C1, C3 | primary | Foy *et al.*, Cancer Prev Res 2016 |
| GSE188981 | GEO | Illumina NovaSeq 6000 | Oral leukoplakia (RNA-seq) | 48 | C1 | validation | GEO |
| GSE30784 | GEO | Affymetrix HG-U133 Plus 2.0 | OSCC + dysplasia + normal | 229 | C2, C3 | primary | Chen *et al.*, Cancer Prev Res 2008 |
| GSE37991 | GEO | Phalanx OneArray | OSCC + adjacent | 80 | C2 | primary | Lee *et al.*, PLoS One 2013 |
| GSE9844  | GEO | Affymetrix HG-U133 Plus 2.0 | OSCC tongue + normal | 38 | C3 | primary | Ye *et al.*, BMC Cancer 2008 |
| GSE6791  | GEO | Affymetrix HG-U133 Plus 2.0 | HNSCC + normal | 84 | C3 | primary | Pyeon *et al.*, Cancer Res 2007 |
| GSE25099 | GEO | Illumina HumanHT-12 v3 | OSCC + normal | 79 | C3 | validation | GEO |
| TCGA-HNSC (oral subset) | GDC | RNA-Seq (HTSeq) | Oral cavity primary tumors + normals | 343 | C3 | primary | TCGA Network, Nature 2015 |
| **Total** | | | | **1049** | | | |

Reference panels (not contributing samples; used in filtering):
Human Protein Atlas v23.0 — cell-surface localization + normal-tissue baseline;
Bgee v15.2 — normal-tissue baseline (independent); EBI Expression Atlas
release-46 (E-MTAB-2836, E-MTAB-3358) — normal-tissue baseline (independent).

## Harmonization

Probe-level intensities were collapsed to HGNC gene symbols using platform-
specific annotation, log2-transformed (with pseudocount 1.0 for RNA-seq), and
within-cohort quantile-normalized. Cross-cohort batch effects were removed
with a ComBat empirical-Bayes adjustment using cohort as the batch variable
and case/control as the protected covariate. Genes present in fewer than 4
cohorts were dropped.

## Three contrasts

- **C1 — high-risk OPMD vs stable OPMD (progression).** Among OPMD samples,
  contrast those that progressed to OSCC within the follow-up window (or
  carried high-grade dysplasia at biopsy) against those that remained stable.
  Effect size: Hedges' g.
- **C2 — lesion epithelium vs adjacent mucosa (paired).** Within each patient,
  paired contrast between lesion and the patient's own adjacent normal mucosa.
  Effect size: paired log2 fold change.
- **C3 — lesion vs normal epithelium (cross-sectional).** Largest contrast in
  sample size; primary lesion-to-normal filter. Effect size: Hedges' g.

Per-cohort effect sizes were pooled across cohorts using random-effects
meta-analysis (DerSimonian–Laird tau²), with Benjamini–Hochberg FDR adjustment
of pooled p-values. Cross-cohort directional consistency was retained as a
filter criterion.

## Selection rules (operationalized in `configs/filters.yaml`)

1. **Surface localization (gating).** Retain only genes annotated by HPA as
   plasma-membrane or cell-surface, with secretome category "Membrane
   proteins" and experimental evidence.
2. **Normal-tissue baseline penalty (gating + multiplicative).** Compute
   per-gene maximum non-oral baseline nTPM across HPA, Bgee, and Atlas. Drop
   if ≥2/3 panels exceed 50 nTPM. Apply a multiplicative penalty
   `1 / (1 + log1p(maxTPM))` when ≥2/3 panels exceed 15 nTPM.
3. **Lesion-to-normal contrast strength.** Pooled |Hedges g| in C3 ≥ 0.50,
   FDR ≤ 0.05, same-direction in ≥4 primary C3 cohorts, with direction
   agreement against C1 and C2 pooled estimates.

## Models

- **L1-regularized Cox** (`sksurv.linear_model.CoxnetSurvivalAnalysis`) on
  the survival outcomes of C1 cohorts. Monotonic sign constraints are
  enforced by sign-flipping each feature so that the prior-expected
  coefficient is non-negative; any negative post-fit coefficient is zeroed,
  guaranteeing that selected features move in a single, biologically
  plausible direction with respect to transformation risk.
- **L1-regularized logistic regression**
  (`sklearn.linear_model.LogisticRegressionCV`) on the binary high-risk-OPMD
  vs stable-OPMD outcome, balanced class weights.
- **Leave-one-cohort-out cross-validation.** Each cohort is held out in turn;
  models are refit on the remaining cohorts. Per-fold coefficients form the
  matrix from which mean |coefficient| and mean |SHAP| are computed.

## Feature importance and ranking

For each candidate gene we compute:

- **|pooled effect|** in C3 (`abs_pooled_c3`).
- **mean |SHAP|** across LOCO folds (for the L1 logistic model; closed-form
  for linear models, equivalent to SHAP runtime output).
- **|Pearson correlation|** between gene expression and the binary risk label.

Each quantity is z-scored across surviving candidates and passed through a
sigmoid. The composite score is a weighted average (0.40, 0.40, 0.20)
multiplied by the per-gene baseline-normal penalty multiplier. The top-N
candidates by composite score form the final ranked list; the promoted top-3
proceeds to wet-lab validation.

## Reproducibility

All hyperparameters live in `configs/*.yaml`. The synthetic-data demo
(`make demo`) regenerates the ranking deterministically under seed 42, and
`tests/test_pipeline_synthetic.py` enforces the top-3 contract in CI. The
real-data pipeline is executed via stage scripts `scripts/01..06_*.py`, with
intermediate parquet artifacts at every stage so any run is resumable.
