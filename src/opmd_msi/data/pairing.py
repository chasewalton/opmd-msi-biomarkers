"""Patient labels for the paired tumor/adjacent samples in GSE37991."""
from __future__ import annotations

import re


def gse37991_patient_id(label: str, fallback: str) -> str:
    """Map T1 OSCC and N1 non-OSCC to one patient in this paired GEO series.

    The cached SOFT design specifies 40 pair-wise tumor/non-tumor biopsies;
    sample titles use T1..T40 and N1..N40. Preserve other IDs unchanged.
    """
    match = re.fullmatch(r"[TN](\d+)(?:\s+.*)?", str(label).strip(), re.IGNORECASE)
    return f"GSE37991_{int(match.group(1))}" if match else fallback
