"""Cohort containers and loaders.

A ``Cohort`` is one study, one omics layer: a samples by features matrix plus the harmonised
metadata rows for those samples. A ``CohortSet`` holds several cohorts of the same layer and
provides the pooled view used by harmonisation and evaluation, with de-duplication of
subjects that appear in more than one study.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .schema import ClinicalSchema, validate_metadata


@dataclass
class Cohort:
    study_id: str
    layer: str
    X: pd.DataFrame  # samples x features, index = sample_id
    meta: pd.DataFrame  # index = sample_id, harmonised schema columns
    feature_space: str = "ensembl_gene"  # ensembl_gene | hgnc | uniprot | refmet | cpg
    is_counts: bool = False  # raw counts (bulk RNA) versus already continuous values
    source: Dict[str, str] = field(default_factory=dict)  # accession, processing provenance

    def __post_init__(self) -> None:
        if not self.X.index.equals(self.meta.index):
            common = self.X.index.intersection(self.meta.index)
            if len(common) == 0:
                raise ValueError(f"{self.study_id}: no overlap between matrix rows and metadata")
            self.X = self.X.loc[common]
            self.meta = self.meta.loc[common]
        self.X = self.X.apply(pd.to_numeric, errors="coerce")

    @property
    def n(self) -> int:
        return self.X.shape[0]

    @property
    def features(self) -> pd.Index:
        return self.X.columns

    def subset(self, mask: pd.Series) -> "Cohort":
        idx = self.X.index[mask.loc[self.X.index].astype(bool).values]
        return Cohort(self.study_id, self.layer, self.X.loc[idx], self.meta.loc[idx], self.feature_space, self.is_counts, dict(self.source))

    @classmethod
    def from_files(
        cls,
        study_id: str,
        layer: str,
        matrix_path: str | Path,
        meta_path: str | Path,
        features_in_rows: bool = True,
        is_counts: bool = False,
        feature_space: str = "ensembl_gene",
        strip_ensembl_version: bool = True,
        **source: str,
    ) -> "Cohort":
        """Load a cohort from a matrix file (csv/tsv/parquet) and a metadata file.

        Matrices downloaded from GEO are usually features in rows; set ``features_in_rows``
        accordingly. Ensembl version suffixes (ENSG00000123456.7) are removed by default so
        that studies quantified against different annotation releases share identifiers.
        """
        X = _read_table(matrix_path)
        if features_in_rows:
            X = X.T
        meta = _read_table(meta_path)
        if "sample_id" in meta.columns:
            meta = meta.set_index("sample_id", drop=False)
        if strip_ensembl_version and feature_space == "ensembl_gene":
            X.columns = [c.split(".")[0] if str(c).startswith("ENSG") else c for c in X.columns]
            X = X.T.groupby(level=0).sum().T if X.columns.duplicated().any() else X
        meta = meta.copy()
        meta["study_id"] = study_id
        meta["omics_layer"] = layer
        validate_metadata(meta, ClinicalSchema(), strict=False)
        return cls(study_id, layer, X, meta, feature_space, is_counts, dict(source))


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=0)


@dataclass
class CohortSet:
    cohorts: List[Cohort]

    def __post_init__(self) -> None:
        layers = {c.layer for c in self.cohorts}
        if len(layers) > 1:
            raise ValueError(f"a CohortSet holds one omics layer; got {layers}")
        ids = [c.study_id for c in self.cohorts]
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate study_id in CohortSet")

    @property
    def layer(self) -> str:
        return self.cohorts[0].layer

    @property
    def study_ids(self) -> List[str]:
        return [c.study_id for c in self.cohorts]

    def __getitem__(self, study_id: str) -> Cohort:
        for c in self.cohorts:
            if c.study_id == study_id:
                return c
        raise KeyError(study_id)

    def shared_features(self, min_fraction_studies: float = 1.0) -> pd.Index:
        """Features present in at least ``min_fraction_studies`` of the cohorts."""
        counts: Dict[str, int] = {}
        for c in self.cohorts:
            for f in c.features:
                counts[f] = counts.get(f, 0) + 1
        need = int(np.ceil(min_fraction_studies * len(self.cohorts)))
        return pd.Index([f for f, k in counts.items() if k >= need])

    def pooled(self, features: Optional[Iterable[str]] = None, fill: float = np.nan) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Stack cohorts into one matrix and one metadata table (study_id column kept)."""
        feats = pd.Index(features) if features is not None else self.shared_features()
        Xs, Ms = [], []
        for c in self.cohorts:
            Xc = c.X.reindex(columns=feats, fill_value=fill)
            Xc.index = [f"{c.study_id}::{s}" for s in Xc.index]
            Mc = c.meta.copy()
            Mc.index = Xc.index
            Xs.append(Xc)
            Ms.append(Mc)
        return pd.concat(Xs), pd.concat(Ms)

    def deduplicate_subjects(self, priority: Optional[List[str]] = None) -> "CohortSet":
        """Drop subjects that appear in more than one study, keeping the first study in ``priority``.

        Consortia such as INCOV/HAARVI, ISARIC/PHOSP and the Mount Sinai biobank feed several
        papers. Pooling the same person twice inflates performance in leave-one-study-out
        evaluation because the held-out subject has a twin in the training set.
        """
        priority = priority or self.study_ids
        seen: set = set()
        out: List[Cohort] = []
        for sid in priority:
            c = self[sid]
            keep = ~c.meta["subject_id"].astype(str).isin(seen)
            seen.update(c.meta["subject_id"].astype(str).tolist())
            out.append(c.subset(keep))
        return CohortSet(out)

    def summary(self) -> pd.DataFrame:
        rows = []
        for c in self.cohorts:
            st = c.meta["lc_status"].value_counts()
            rows.append(
                {
                    "study_id": c.study_id,
                    "layer": c.layer,
                    "n_samples": c.n,
                    "n_features": len(c.features),
                    "LC": int(st.get("LC", 0)),
                    "COVID_recovered": int(st.get("COVID_recovered", 0)),
                    "never_infected": int(st.get("never_infected", 0)),
                    "matrix": ",".join(sorted(set(c.meta.get("sample_matrix", pd.Series(dtype=str)).dropna().astype(str)))),
                    "accession": c.source.get("accession", ""),
                }
            )
        return pd.DataFrame(rows)
