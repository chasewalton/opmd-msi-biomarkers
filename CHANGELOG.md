# Changelog

## 0.1.1 - 2026-10-23

- Add real-data and saved-run commands for ENT biomarker prioritization.
- Apply reviewed epithelial protein evidence, available humanized antibodies,
  and normal oral protein background to eligibility and ranking.
- Exclude numerically constant expression signals and retain valid
  single-cohort effects.
- Normalize paired GSE37991 tumor/adjacent labels to a shared patient ID,
  including previously parsed cohort files.
- Add eligibility tests and fix lint errors in the real-data modules.
