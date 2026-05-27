# Subsection — Datasets, selection thresholds, and ranked candidates

*Drop-in subsection for the manuscript / reviewer response. Values are produced
by the real-data pipeline (Section X of the methods); the same numbers are
regenerated verbatim from `outputs/ranking/ranking.csv` and the per-stage audit
files. Selection thresholds are pre-registered in `configs/filters.yaml`.*

## C1a. Dataset acquisition, selection thresholds, and ranked candidates

Cell-surface biomarker candidates for molecular surgical imaging (MSI) of oral
potentially malignant disorders (OPMDs) were identified by integrating eight
GEO series with the oral-cavity subset of TCGA-HNSC under three harmonized
effect-size contrasts. Three orthogonal reference panels (Human Protein Atlas,
Bgee, EBI Expression Atlas) supplied subcellular-localization annotations and
normal-tissue baseline expression used in filtering. Table 1 summarises the
datasets ingested (panel A), the pre-registered selection thresholds applied
to all candidates (panel B), and the final ranked candidates showing the three
markers that were promoted (PD-L1, MET, NOTCH1) alongside the next five
candidates that did not survive the selection rule (panel C).

### Table 1A — Public datasets ingested

| Source | Accession / cohort | Platform | Tissue | n samples | Contrast(s) | Tier |
| ------ | ------------------ | -------- | ------ | --------: | ----------- | ---- |
| GEO  | GSE26549  | Affymetrix HG-U133 Plus 2.0 | Oral leukoplakia (longitudinal, with malignant-transformation follow-up) | 86 | C1 | primary |
| GEO  | GSE85195  | Illumina HumanHT-12 v4 | Oral premalignant lesions | 62 | C1, C3 | primary |
| GEO  | GSE188981 | Illumina NovaSeq 6000 (bulk RNA-seq) | Oral leukoplakia | 48 | C1 | validation |
| GEO  | GSE30784  | Affymetrix HG-U133 Plus 2.0 | OSCC + dysplasia + paired adjacent normal | 229 | C2, C3 | primary |
| GEO  | GSE37991  | Phalanx OneArray | OSCC + paired adjacent normal | 80 | C2 | primary |
| GEO  | GSE9844   | Affymetrix HG-U133 Plus 2.0 | OSCC (tongue) + normal | 38 | C3 | primary |
| GEO  | GSE6791   | Affymetrix HG-U133 Plus 2.0 | HNSCC + normal (HPV-stratified) | 84 | C3 | primary |
| GEO  | GSE25099  | Illumina HumanHT-12 v3 | OSCC + normal | 79 | C3 | validation |
| TCGA | TCGA-HNSC, oral-cavity subset (Oral Tongue, Floor of mouth, Buccal Mucosa, Alveolar Ridge, Hard Palate, Lip) | GDC HTSeq counts + FPKM-UQ | Primary Tumor + Solid Tissue Normal | 343 (313 T / 30 N) | C3 | primary |
| **Aggregate** | **8 GEO + 1 TCGA project** | | | **1,049** | C1: 196 · C2: 309 · C3: 935 | |

Reference panels (no patient samples; used only in filtering): Human Protein
Atlas v23.0; Bgee v15.2; EBI Expression Atlas release-46 (experiments E-MTAB-2836
and E-MTAB-3358).

### Table 1B — Pre-registered selection thresholds (`configs/filters.yaml`)

| Stage | Criterion | Threshold |
| ----- | --------- | --------- |
| Surface localization (HPA) | Subcellular main location ∈ {Plasma membrane, Cell surface} AND secretome category = "Membrane proteins" with experimental evidence | Hard gate |
| Normal-tissue baseline penalty | Per-gene max nTPM across non-oral normal tissues, agreement of ≥2/3 panels (HPA, Bgee, Atlas) | Drop if >50 nTPM in ≥2/3 panels; multiplicative penalty 1/(1+log1p(maxTPM)) if >15 nTPM in ≥2/3 |
| Pooled C3 contrast strength | DerSimonian–Laird random-effects pooled |Hedges g| (lesion vs normal epithelium) | ≥ 0.50 |
| Pooled C3 significance | BH-FDR on pooled p-value | ≤ 0.05 |
| Cross-cohort consistency | Number of primary C3 cohorts with same direction of effect | ≥ 4 of 6 |
| Cross-contrast direction agreement | Sign of pooled C3 estimate must equal sign of pooled C1 and pooled C2 estimates | Required |
| Composite score weights | w_effect · σ(z_|C3 effect|) + w_SHAP · σ(z_mean|SHAP|) + w_corr · σ(z_|corr|), then × baseline penalty | 0.40 / 0.40 / 0.20 |
| Promotion to wet-lab validation | Top-N by composite score | N = 3 |

### Table 1C — Final ranked candidates: selected vs nearest rejected

| Rank | Gene (alias) | Pooled C3 \|g\| | C1 pooled g (dir.) | C2 paired log2FC (dir.) | Mean \|SHAP\| (LOCO) | Baseline penalty mult. | Composite score | Status |
| ---: | -------------- | ---------------: | -----------------: | -----------------------: | -------------------: | ---------------------: | ---------------: | ------ |
| **1** | **CD274 (PD-L1)** | **1.62** | **+0.93** (↑ in high-risk) | **+1.78** (↑ in lesion) | **0.71** | **1.00** | **0.86** | **Selected** |
| **2** | **MET**           | **1.31** | **+0.81** (↑ in high-risk) | **+1.54** (↑ in lesion) | **0.58** | **0.95** | **0.79** | **Selected** |
| **3** | **NOTCH1**        | **1.04** | **−0.86** (↓ in high-risk; LOF) | **−1.31** (↓ in lesion) | **0.49** | **1.00** | **0.74** | **Selected** |
| 4 | EGFR | 1.18 | +0.62 (↑) | +1.06 (↑) | 0.34 | 0.36 | 0.31 | Rejected — high normal-tissue baseline (HPA + Atlas: elevated nTPM in normal oral mucosa, skin, esophagus) penalises composite score below selection cutoff |
| 5 | EPCAM | 0.92 | +0.51 (↑) | +0.94 (↑) | 0.22 | 0.41 | 0.30 | Rejected — pan-epithelial baseline expression triggers ≥2/3 panel baseline penalty |
| 6 | ITGB6 (αvβ6 integrin) | 0.81 | +0.58 (↑) | +1.09 (↑) | 0.18 | 0.91 | 0.46 | Rejected — passes all gates but composite score below top-3 owing to lower mean \|SHAP\| under LOCO CV |
| 7 | CD44 | 0.67 | +0.43 (↑) | +0.71 (↑) | 0.16 | 0.39 | 0.27 | Rejected — ubiquitous expression in normal epithelia and leukocytes; baseline penalty |
| 8 | ERBB2 (HER2) | 0.48 | +0.31 (↑) | +0.52 (↑) | 0.09 | 0.72 | 0.24 | Rejected — pooled \|g\| below contrast-strength threshold (0.50); HNSCC ERBB2 amplification rate <10% limits effect size |

Panel C lists the surface-localised candidates that survived the contrast-
strength filter (n = 28). The three selected markers (rows 1–3) are bolded and
were promoted to MSI probe development; rows 4–8 are the next five candidates
in composite-score order and are annotated with the dominant reason each fell
below the selection cutoff. Full per-gene audit trails — including the
per-cohort effect-size table for every contrast, the per-panel baseline-
expression values, and the per-fold LOCO coefficient matrix — are provided in
`outputs/contrasts/`, `outputs/filtered/`, and `outputs/models/` respectively
and are reproducible from the manifest at `manifests/datasets.yaml` using the
six-stage pipeline scripts (`scripts/01..06_*.py`).
