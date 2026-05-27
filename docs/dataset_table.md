# Dataset table (rendered from manifests/datasets.yaml)

This page is generated from `manifests/datasets.yaml` to give a publication-
ready summary of every public dataset the pipeline ingests.

## GEO series

| Accession | Platform | Tissue | n | Contrast(s) | Tier | Reference |
| --------- | -------- | ------ | -: | ----------- | ---- | --------- |
| GSE26549 | Affymetrix HG-U133 Plus 2.0 | Oral leukoplakia (longitudinal) | 86 | C1 | primary | Saintigny *et al.* Cancer Prev Res 2011 |
| GSE85195 | Illumina HumanHT-12 v4 | Oral premalignant lesions | 62 | C1, C3 | primary | Foy *et al.* Cancer Prev Res 2016 |
| GSE188981 | Illumina NovaSeq 6000 (bulk RNA-seq) | Oral leukoplakia | 48 | C1 | validation | GEO 2022 |
| GSE30784 | Affymetrix HG-U133 Plus 2.0 | OSCC + dysplasia + normal | 229 | C2, C3 | primary | Chen *et al.* Cancer Prev Res 2008 |
| GSE37991 | Phalanx OneArray | OSCC + adjacent | 80 | C2 | primary | Lee *et al.* PLoS One 2013 |
| GSE9844  | Affymetrix HG-U133 Plus 2.0 | OSCC tongue + normal | 38 | C3 | primary | Ye *et al.* BMC Cancer 2008 |
| GSE6791  | Affymetrix HG-U133 Plus 2.0 | HNSCC + normal | 84 | C3 | primary | Pyeon *et al.* Cancer Res 2007 |
| GSE25099 | Illumina HumanHT-12 v3 | OSCC + normal | 79 | C3 | validation | GEO 2011 |

**GEO total: 8 series, 706 samples.**

## TCGA

| Project | Data type | Sub-cohort | n used | Contrast | Tier | Reference |
| ------- | --------- | ---------- | -----: | -------- | ---- | --------- |
| TCGA-HNSC | RNA-Seq (HTSeq counts + FPKM-UQ) | Oral cavity (Oral Tongue, Floor of mouth, Buccal Mucosa, Alveolar Ridge, Hard Palate, Lip), Primary Tumor + Solid Tissue Normal | 343 (313 tumor / 30 normal) | C3 | primary | TCGA Network. Nature 2015;517:576-82. |

**TCGA total: 1 project, 343 samples (oral subset of 528 total TCGA-HNSC samples).**

## Reference panels (not contributing samples; used in filters)

| Panel | Release | Role | Fields consumed |
| ----- | ------- | ---- | --------------- |
| Human Protein Atlas | 23.0 | Cell-surface localization gate + normal-tissue baseline | Subcellular main location, Secretome location, Tissue RNA (nTPM) |
| Bgee | 15.2 | Normal-tissue baseline (independent of HPA) | Per-gene max nTPM in non-oral tissues |
| EBI Expression Atlas | release-46 | Normal-tissue baseline (independent of HPA + Bgee) | E-MTAB-2836 (32-tissue panel), E-MTAB-3358 (FANTOM5 subset) |

## Aggregate

- **Datasets:** 9 primary expression + 3 reference panels.
- **Patients / samples used:** 1049 human samples.
- **Contrasts:** 3 (C1: 196 samples; C2: 309 samples; C3: 935 samples).
- **Validation cohorts (LOCO held-out):** GSE188981 (C1), GSE25099 (C3).
