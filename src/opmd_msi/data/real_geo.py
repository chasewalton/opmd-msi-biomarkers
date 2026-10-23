"""Real-data parser for GEO cohorts.

For each cohort:
  1. Load the cached SOFT family file.
  2. Parse per-sample expression from the platform-specific probe table.
  3. Collapse probe IDs to HGNC gene symbols using the GPL annotation that
     ships inside the same SOFT file.
  4. Harmonize sample metadata to the schema in manifests/cohorts.yaml.
  5. Write expression.parquet + metadata.parquet under data/raw/GEO/<ACC>/.

This file is intentionally cohort-specific where the metadata fields are
bespoke (the only way to map real GEO metadata is per-cohort).
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import GEOparse
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Probe → HGNC collapse
# ---------------------------------------------------------------------------

# Each GPL annotation table has a different column name for the gene symbol.
# We search columns in this priority order; first match wins.
_SYMBOL_COL_CANDIDATES: tuple[str, ...] = (
    "Gene Symbol", "Symbol", "gene_assignment", "Gene_Symbol",
    "ILMN_Gene", "Entrez_Gene_ID", "GENE_SYMBOL", "gene",
    "Gene symbol", "GeneSymbol", "GENE",
)


def _extract_symbol(value: str | float | None) -> str | None:
    """Pull the first HGNC-looking symbol from a possibly-decorated cell."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    s = str(value).strip()
    if not s or s.lower() in {"---", "na", "nan", "none"}:
        return None
    # Some platforms (Affy Human Gene 1.0 ST) pack multiple records as
    # "RefSeq // Symbol // Description // location" separated by " /// "
    if "//" in s:
        # take first sub-record's middle field if it looks like a symbol
        first = s.split("///")[0]
        parts = [p.strip() for p in first.split("//")]
        for p in parts:
            if re.fullmatch(r"[A-Z][A-Z0-9\-]{1,15}", p):
                return p
    # Some platforms use "GENE1 /// GENE2" between symbols directly
    if "///" in s:
        first = s.split("///")[0].strip()
        if re.fullmatch(r"[A-Z][A-Z0-9\-]{1,15}", first):
            return first
    if re.fullmatch(r"[A-Z][A-Z0-9\-]{1,15}", s):
        return s
    return None


def collapse_probes_to_symbols(
    expression_probes: pd.DataFrame,   # samples x probes
    gpl_table: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse probe-level expression to HGNC symbols by taking the per-symbol mean."""
    # find the symbol column
    sym_col = next((c for c in _SYMBOL_COL_CANDIDATES if c in gpl_table.columns), None)
    if sym_col is None:
        raise ValueError(f"No symbol column among candidates {_SYMBOL_COL_CANDIDATES}; "
                         f"available: {list(gpl_table.columns)[:10]}…")
    id_col = "ID" if "ID" in gpl_table.columns else gpl_table.columns[0]
    gpl = gpl_table[[id_col, sym_col]].copy()
    gpl[sym_col] = gpl[sym_col].apply(_extract_symbol)
    gpl = gpl.dropna(subset=[sym_col])

    # restrict probes to those with a mapped symbol
    keep = set(gpl[id_col].astype(str)) & set(expression_probes.columns.astype(str))
    if not keep:
        # try matching numeric vs string
        keep = set(gpl[id_col].astype(str)) & {str(c) for c in expression_probes.columns}
    if not keep:
        raise ValueError("No probes in expression matrix matched the GPL annotation IDs.")

    sub = expression_probes.loc[:, expression_probes.columns.astype(str).isin(keep)]
    # build probe → symbol mapping for the surviving probes
    mapping = dict(zip(gpl[id_col].astype(str), gpl[sym_col], strict=True))
    sub.columns = [mapping[str(c)] for c in sub.columns]
    # collapse multi-mapped symbols by mean
    collapsed = sub.T.groupby(level=0).mean().T
    return collapsed


# ---------------------------------------------------------------------------
# Per-cohort metadata harmonizers
# ---------------------------------------------------------------------------

def _kv_from_characteristics(g) -> dict[str, str]:
    """Parse 'key: value' entries from GEO characteristics_ch1 into a dict."""
    out: dict[str, str] = {}
    for raw in g.metadata.get("characteristics_ch1", []):
        if ":" in raw:
            k, v = raw.split(":", 1)
            out[k.strip().lower()] = v.strip()
    return out


def _meta_GSE26549(gsm_id: str, g) -> dict:
    c = _kv_from_characteristics(g)
    outcome = (c.get("outcome") or "").lower().strip()
    # Values observed in GSE26549: "oral cancer development" vs "no oral cancer development".
    transformed = 1 if outcome.startswith("oral cancer development") else 0
    survival_yrs = c.get("oral cancer-free survival time (years)") or c.get("oral cancer-free survival time")
    try:
        t2t_months = float(survival_yrs) * 12.0 if survival_yrs not in (None, "") else np.nan
    except ValueError:
        t2t_months = np.nan
    hist = (c.get("histology at baseline (hyperplasia versus dysplasia)") or "").lower()
    if "dysplasia" in hist:
        tissue_class = "high_grade_dysplasia" if "high" in hist else "low_grade_dysplasia"
    else:
        tissue_class = "low_grade_dysplasia"  # treat hyperplasia as low-grade per Saintigny grouping
    # Risk class: transformation status defines high-risk OPMD
    risk_class = "high_risk_OPMD" if transformed else "stable_OPMD"
    age = c.get("age")
    try:
        age_years = float(age) if age else np.nan
    except ValueError:
        age_years = np.nan
    smoking = (c.get("smoking habits") or "").lower()
    smoker = 1 if smoking in {"current", "former"} else (0 if smoking == "never" else np.nan)
    return {
        "sample_id": gsm_id,
        "cohort_id": "GSE26549",
        "patient_id": gsm_id,           # no patient ID exposed; treat samples as independent
        "tissue_class": tissue_class,
        "risk_class": risk_class,
        "platform": "GPL6244",          # Affymetrix Human Gene 1.0 ST
        "assay_type": "microarray",
        "sex": "F" if c.get("sex","").lower()=="female" else ("M" if c.get("sex","").lower()=="male" else "unknown"),
        "age_years": age_years,
        "smoker": smoker,
        "hpv_status": "unknown",
        "follow_up_months": t2t_months,
        "transformed": transformed,
        "time_to_transformation_months": t2t_months,
    }


def _meta_GSE85195(gsm_id: str, g) -> dict:
    """Foy 2016 — Gingivobuccal OPL + early-stage OSCC + 1 normal.
    Sample class is in source_name_ch1, NOT in characteristics_ch1.
    """
    sn = (g.metadata.get("source_name_ch1", [""])[0] or "").lower()
    if "tumor" in sn or "oscc" in sn:
        tissue_class, risk_class = "OSCC", "OSCC"
    elif "leukoplakia" in sn or "precancerous" in sn or "opl" in sn:
        tissue_class, risk_class = "high_grade_dysplasia", "high_risk_OPMD"
    elif "normal" in sn:
        tissue_class, risk_class = "normal", "normal"
    else:
        tissue_class, risk_class = "normal", "normal"
    c = _kv_from_characteristics(g)
    age_str = (c.get("age") or "").replace("y", "").strip()
    try:
        age_years = float(age_str)
    except ValueError:
        age_years = np.nan
    gender = (c.get("gender") or "").lower()
    sex = "F" if gender == "female" else ("M" if gender == "male" else "unknown")
    return {
        "sample_id": gsm_id, "cohort_id": "GSE85195", "patient_id": gsm_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL6480", "assay_type": "microarray",
        "sex": sex, "age_years": age_years, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


def _meta_GSE30784(gsm_id: str, g) -> dict:
    """Chen 2008 — OSCC + dysplasia + normal control. Status field gives class.
    No patient ID available for pairing; treated as independent C3 samples.
    """
    c = _kv_from_characteristics(g)
    status = (c.get("status") or "").lower()
    if status == "cancer":
        tissue_class, risk_class = "OSCC", "OSCC"
    elif status == "dysplasia":
        tissue_class, risk_class = "high_grade_dysplasia", "high_risk_OPMD"
    elif status == "control":
        tissue_class, risk_class = "normal", "normal"
    else:
        tissue_class, risk_class = "normal", "normal"
    age = (c.get("age") or "").strip()
    # Format e.g. "40-49"; take midpoint
    if "-" in age:
        try:
            lo, hi = age.split("-")
            age_years = (float(lo) + float(hi)) / 2
        except ValueError:
            age_years = np.nan
    else:
        try:
            age_years = float(age)
        except ValueError:
            age_years = np.nan
    sex_raw = (c.get("sex") or "").upper().strip()
    sex = "F" if sex_raw.startswith("F") else ("M" if sex_raw.startswith("M") else "unknown")
    return {
        "sample_id": gsm_id, "cohort_id": "GSE30784", "patient_id": gsm_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL570", "assay_type": "microarray",
        "sex": sex, "age_years": age_years, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


def _meta_GSE37991(gsm_id: str, g) -> dict:
    """Lee 2013 — 40 OSCC + 40 paired adjacent. Pair via title "T<i> OSCC" / "T<i> Adjacent".
    """
    sn = (g.metadata.get("source_name_ch1", [""])[0] or "").lower()
    title = (g.metadata.get("title", [""])[0] or "")
    # Title is like "T17 OSCC" or "T17 Adjacent"; patient_id = T17
    patient_id = title.split()[0] if title else gsm_id
    if "adjacent" in sn:
        tissue_class, risk_class = "adjacent_mucosa", "normal"
    elif "oscc" in sn or "carcinoma" in sn:
        tissue_class, risk_class = "OSCC", "OSCC"
    else:
        tissue_class, risk_class = "normal", "normal"
    return {
        "sample_id": gsm_id, "cohort_id": "GSE37991", "patient_id": patient_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL6883", "assay_type": "microarray",
        "sex": "unknown", "age_years": np.nan, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


def _meta_GSE9844(gsm_id: str, g) -> dict:
    """Ye 2008 — 26 oral tongue SCC + 12 normal. Class is in source_name_ch1.
    Characteristics are unstructured "Sex M, Age 62" — parse with regex.
    """
    sn = (g.metadata.get("source_name_ch1", [""])[0] or "").lower()
    if "normal" in sn or "control" in sn:
        tissue_class, risk_class = "normal", "normal"
    else:
        tissue_class, risk_class = "OSCC", "OSCC"
    raw_chars = " ".join(g.metadata.get("characteristics_ch1", []))
    sex_match = re.search(r"sex\s*[:,]?\s*([MF])", raw_chars, re.IGNORECASE)
    age_match = re.search(r"age\s*[:,]?\s*(\d+)", raw_chars, re.IGNORECASE)
    sex = sex_match.group(1).upper() if sex_match else "unknown"
    try:
        age_years = float(age_match.group(1)) if age_match else np.nan
    except ValueError:
        age_years = np.nan
    return {
        "sample_id": gsm_id, "cohort_id": "GSE9844", "patient_id": gsm_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL570", "assay_type": "microarray",
        "sex": sex, "age_years": age_years, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


def _meta_GSE6791(gsm_id: str, g) -> dict:
    """Pyeon 2007 — HNSCC + Cervical, with normals. We filter to H&N only
    via source_name_ch1; cervical samples are returned with risk_class='_drop'
    so they can be filtered out by the caller.
    """
    sn = (g.metadata.get("source_name_ch1", [""])[0] or "").lower()
    if "head and neck" not in sn:
        return {
            "sample_id": gsm_id, "cohort_id": "GSE6791", "patient_id": gsm_id,
            "tissue_class": "_drop", "risk_class": "_drop",
            "platform": "GPL570", "assay_type": "microarray",
            "sex": "unknown", "age_years": np.nan, "smoker": np.nan,
            "hpv_status": "unknown",
            "follow_up_months": np.nan, "transformed": np.nan,
            "time_to_transformation_months": np.nan,
        }
    if "normal" in sn:
        tissue_class, risk_class = "normal", "normal"
    else:
        tissue_class, risk_class = "OSCC", "OSCC"
    c = _kv_from_characteristics(g)
    try:
        age_years = float(c.get("age") or "nan")
    except ValueError:
        age_years = np.nan
    gender = (c.get("gender") or "").lower()
    sex = "F" if gender.startswith("f") else ("M" if gender.startswith("m") else "unknown")
    return {
        "sample_id": gsm_id, "cohort_id": "GSE6791", "patient_id": gsm_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL570", "assay_type": "microarray",
        "sex": sex, "age_years": age_years, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


def _meta_GSE25099(gsm_id: str, g) -> dict:
    """Reis 2011 — 57 OSCC + 22 normal. Class in source_name_ch1."""
    sn = (g.metadata.get("source_name_ch1", [""])[0] or "").lower()
    if "normal" in sn:
        tissue_class, risk_class = "normal", "normal"
    else:
        tissue_class, risk_class = "OSCC", "OSCC"
    return {
        "sample_id": gsm_id, "cohort_id": "GSE25099", "patient_id": gsm_id,
        "tissue_class": tissue_class, "risk_class": risk_class,
        "platform": "GPL5175", "assay_type": "microarray",
        "sex": "unknown", "age_years": np.nan, "smoker": np.nan,
        "hpv_status": "unknown",
        "follow_up_months": np.nan, "transformed": np.nan,
        "time_to_transformation_months": np.nan,
    }


COHORT_HARMONIZERS = {
    "GSE26549": _meta_GSE26549,
    "GSE85195": _meta_GSE85195,
    "GSE30784": _meta_GSE30784,
    "GSE37991": _meta_GSE37991,
    "GSE9844":  _meta_GSE9844,
    "GSE6791":  _meta_GSE6791,
    "GSE25099": _meta_GSE25099,
}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def process_cohort(accession: str, raw_dir: str | Path) -> dict[str, str]:
    """Parse one cohort end-to-end. Returns paths to written parquet files."""
    raw_dir = Path(raw_dir)
    soft = next(raw_dir.glob(f"{accession}_family.soft*"), None) or \
           next(raw_dir.glob(f"{accession}*.soft*"), None)
    if soft is None:
        raise FileNotFoundError(f"No SOFT file for {accession} under {raw_dir}")
    gse = GEOparse.get_GEO(filepath=str(soft), silent=True)

    # 1. assemble per-sample expression
    sample_series: list[pd.Series] = []
    for gsm_id, gsm in gse.gsms.items():
        if gsm.table is None or gsm.table.empty:
            log.warning("Empty table for %s", gsm_id)
            continue
        s = gsm.table.set_index("ID_REF")["VALUE"].astype(float)
        s.name = gsm_id
        sample_series.append(s)
    if not sample_series:
        raise RuntimeError(f"No sample tables in {accession}")
    expression_probes = pd.concat(sample_series, axis=1).T   # samples x probes

    # 2. probe → HGNC collapse using the GPL annotation in the SOFT file
    gpl_id = list(gse.gpls.keys())[0]
    gpl_table = gse.gpls[gpl_id].table
    expression_symbols = collapse_probes_to_symbols(expression_probes, gpl_table)

    # 3. harmonize metadata
    harmonizer = COHORT_HARMONIZERS.get(accession)
    if harmonizer is None:
        raise NotImplementedError(f"No metadata harmonizer registered for {accession}")
    meta_rows = [harmonizer(gsm_id, gsm) for gsm_id, gsm in gse.gsms.items()]
    metadata = pd.DataFrame(meta_rows).set_index("sample_id")
    # Drop samples flagged "_drop" by per-cohort harmonizer (e.g. GSE6791 cervical).
    metadata = metadata[metadata["risk_class"] != "_drop"].copy()

    # restrict expression to harmonized samples
    expression_symbols = expression_symbols.reindex(metadata.index)

    # 4. write parquet
    expr_path = raw_dir / f"{accession}.expression.parquet"
    meta_path = raw_dir / f"{accession}.metadata.parquet"
    expression_symbols.to_parquet(expr_path)
    metadata.to_parquet(meta_path)
    log.info("Wrote %s: %d samples x %d genes", accession,
             expression_symbols.shape[0], expression_symbols.shape[1])
    return {
        "accession": accession,
        "expression_path": str(expr_path),
        "metadata_path": str(meta_path),
        "n_samples": expression_symbols.shape[0],
        "n_genes": expression_symbols.shape[1],
        "platform": gpl_id,
    }
