"""Real-data parser for TCGA-HNSC oral subset (STAR counts -> harmonized parquet)."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def _parse_star_counts(path: Path) -> pd.Series:
    """Read one STAR-counts TSV from GDC. Returns Series indexed by Ensembl gene id."""
    # Format: skip the first 6 header lines (## comments + column names);
    # columns: gene_id, gene_name, gene_type, unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded
    # The first 4 rows after the header are N_unmapped, N_multimapping, etc.
    df = pd.read_csv(path, sep="\t", skiprows=1)
    # remove summary rows that start with N_
    df = df[~df["gene_id"].astype(str).str.startswith("N_")]
    # Use tpm_unstranded as expression measure; log2(tpm+1) downstream
    s = df.set_index("gene_name")["tpm_unstranded"].astype(float)
    # collapse duplicate gene_name (multiple Ensembl IDs per symbol) by mean
    s = s.groupby(level=0).mean()
    return s


def process_tcga_hnsc_oral(raw_dir: str | Path) -> dict:
    """Parse all downloaded TCGA-HNSC oral RNA-seq STAR counts files.

    Returns dict with paths to written expression.parquet + metadata.parquet.
    """
    raw_dir = Path(raw_dir)
    files_dir = raw_dir / "files"
    manifest = json.loads((raw_dir / "rna_seq_files.json").read_text())
    # exclude metastatic
    def sample_types(f): return [s.get("sample_type") for c in f.get("cases", []) for s in c.get("samples", [])]
    manifest = [f for f in manifest if "Metastatic" not in sample_types(f)]

    series: list[pd.Series] = []
    meta_rows: list[dict] = []
    n_ok = 0
    n_miss = 0
    for f in manifest:
        fid = f["file_id"]
        # match any file in files_dir starting with fid
        candidates = list(files_dir.glob(f"{fid}*"))
        if not candidates:
            n_miss += 1
            continue
        try:
            s = _parse_star_counts(candidates[0])
        except Exception as e:
            log.warning("parse fail %s: %s", fid, e)
            n_miss += 1
            continue
        sample_id = f"TCGA-{fid[:8]}"
        s.name = sample_id
        series.append(s)
        case = f["cases"][0]
        sample = case["samples"][0]
        st = sample.get("sample_type", "")
        too = case.get("diagnoses", [{}])[0].get("tissue_or_organ_of_origin", "")
        if st == "Primary Tumor":
            tissue_class, risk_class = "OSCC", "OSCC"
        elif st == "Solid Tissue Normal":
            tissue_class, risk_class = "normal", "normal"
        else:
            tissue_class, risk_class = "_drop", "_drop"
        meta_rows.append({
            "sample_id": sample_id,
            "cohort_id": "TCGA-HNSC-oral",
            "patient_id": case.get("submitter_id", sample_id),
            "tissue_class": tissue_class, "risk_class": risk_class,
            "platform": "RNA-seq (STAR counts; TPM)",
            "assay_type": "bulk_rnaseq",
            "sex": "unknown", "age_years": np.nan, "smoker": np.nan,
            "hpv_status": "unknown",
            "follow_up_months": np.nan, "transformed": np.nan,
            "time_to_transformation_months": np.nan,
            "tissue_or_organ_of_origin": too,
        })
        n_ok += 1

    log.info("TCGA-HNSC-oral: parsed=%d missing=%d", n_ok, n_miss)
    # Build expression matrix: samples x genes
    expression = pd.concat(series, axis=1).T
    # log2(tpm + 1)
    expression = np.log2(expression.fillna(0) + 1.0)
    metadata = pd.DataFrame(meta_rows).set_index("sample_id")
    metadata = metadata[metadata["risk_class"] != "_drop"]
    expression = expression.reindex(metadata.index)

    out_expr = raw_dir / "TCGA-HNSC-oral.expression.parquet"
    out_meta = raw_dir / "TCGA-HNSC-oral.metadata.parquet"
    expression.to_parquet(out_expr)
    metadata.to_parquet(out_meta)
    return {
        "accession": "TCGA-HNSC-oral",
        "expression_path": str(out_expr),
        "metadata_path": str(out_meta),
        "n_samples": expression.shape[0],
        "n_genes": expression.shape[1],
        "n_tumor": int((metadata["tissue_class"] == "OSCC").sum()),
        "n_normal": int((metadata["tissue_class"] == "normal").sum()),
    }
