"""Network fetchers for GEO, TCGA, HPA, Bgee, Expression Atlas.

Each fetcher is invoked by `scripts/01_fetch_data.py` and writes raw artifacts
to `data/raw/<source>/<accession>/`. Fetchers do not transform data — they
preserve provenance. Harmonization is performed in `data/harmonize.py`.
"""
