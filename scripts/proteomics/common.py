"""Shared output contract of the proteomics cohort builders (scripts/proteomics/build_<study>.py).

Each builder reads data/proteomics/<study>/raw and calls ``write_cohort``, which writes

  data/proteomics/<study>/<study>_matrix.csv    proteins (UniProt key, rows) x samples (columns)
  data/proteomics/<study>/<study>_features.csv  the adapter's feature table (key, source ids)
  data/meta/<study>_meta.csv                    one row per matrix column, ClinicalSchema columns
                                                first, cohort specific columns after them

and returns the summary row used by docs/proteomics_cohorts.md. The metadata must pass
``chorus.schema.validate_metadata`` (strict) and have exactly the matrix columns as sample_id.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from chorus.adapters import Adapted
from chorus.schema import ClinicalSchema, validate_metadata

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "proteomics"
META = ROOT / "data" / "meta"


def swissprot_accessions(path: Path = ROOT / "data" / "ref" / "hgnc_complete_set.txt") -> set:
    """Reviewed (Swiss-Prot) UniProt accessions of human protein coding genes, from the HGNC complete
    set; passed to the ms_proteomics adapter as ``preferred_accessions``."""
    t = pd.read_csv(path, sep="\t", usecols=["uniprot_ids"], dtype=str).dropna()
    return {a.strip() for ids in t["uniprot_ids"] for a in ids.split("|") if a.strip()}


def raw_dir(study: str) -> Path:
    return DATA / study / "raw"


def write_cohort(study: str, adapted: Adapted, meta: pd.DataFrame, extra: Optional[Dict] = None) -> Dict:
    """Validate and write one cohort. ``adapted.X`` is samples x proteins; ``meta`` is indexed
    by sample_id in the same order. Returns the summary row."""
    schema = ClinicalSchema()
    X = adapted.X
    if not X.index.equals(pd.Index(meta.index)):
        raise ValueError(f"{study}: metadata rows do not match the matrix samples")
    meta = meta.copy()
    meta["sample_id"] = meta.index
    meta["study_id"] = study
    meta["omics_layer"] = adapted.layer
    cols = [c for c in schema.columns if c in meta.columns] + [c for c in meta.columns if c not in schema.columns]
    for c in schema.required:
        if c not in meta.columns:
            raise ValueError(f"{study}: metadata lacks required column {c}")
    meta = meta[cols]
    validate_metadata(meta, schema, strict=True)
    out = DATA / study
    out.mkdir(parents=True, exist_ok=True)
    M = X.T
    M.index.name = "uniprot"
    M.to_csv(out / f"{study}_matrix.csv")
    adapted.features.to_csv(out / f"{study}_features.csv", index=False)
    META.mkdir(parents=True, exist_ok=True)
    meta.to_csv(META / f"{study}_meta.csv", index=False)
    return summary_row(study, X, meta, adapted, extra)


def summary_row(study: str, X: pd.DataFrame, meta: pd.DataFrame, adapted: Adapted, extra: Optional[Dict] = None) -> Dict:
    phase = meta.get("sample_phase", pd.Series("unknown", index=meta.index)).astype(str)
    st = meta["lc_status"].astype(str)
    post = phase != "acute"
    row = {
        "study_id": study,
        "layer": adapted.layer,
        "n_samples": len(meta),
        "n_subjects": meta["subject_id"].nunique(),
        "n_LC": int(((st == "LC") & post).sum()),
        "n_COVID_recovered": int(((st == "COVID_recovered") & post).sum()),
        "n_never_infected": int((st == "never_infected").sum()),
        "n_acute": int((phase == "acute").sum()),
        "n_other": int((~st.isin(["LC", "COVID_recovered", "never_infected"]) & post).sum()),
        "n_LC_subjects": meta.loc[(st == "LC") & post, "subject_id"].nunique(),
        "n_COVID_recovered_subjects": meta.loc[(st == "COVID_recovered") & post, "subject_id"].nunique(),
        "n_proteins": X.shape[1],
        "sample_matrix": ",".join(sorted(meta["sample_matrix"].dropna().astype(str).unique())),
        "platform": ",".join(sorted(meta.get("platform", pd.Series(dtype=str)).dropna().astype(str).unique())),
    }
    row.update(extra or {})
    return row
