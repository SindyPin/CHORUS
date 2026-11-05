"""Plasma and serum proteomics and metabolomics preprocessing before cross study harmonisation.

Affinity (Olink PEA). NPX values are relative log2 units whose absolute level differs
between panels and projects; ratios between samples within a project are preserved (Olink
FAQ). Without bridge samples the vendor states that datasets must not be placed on one NPX
scale (OlinkAnalyze bridging vignette), and consortia (SCALLOP) pool effect sizes, not NPX.
CHORUS therefore (i) restricts to assays shared across panel generations by UniProt ID,
(ii) drops assays with a high below-LOD fraction, and (iii) applies within study
standardisation or rank based operators, treating any pooled model as a rank level
signature rather than a concentration model.

Mass spectrometry. Prefer uniform re-search of raw files (quantms, DIA-NN, FragPipe) with
one FASTA and one FDR procedure. Correct at protein level with median centring; avoid
imputing then running ComBat (HarmonizR, proBatch). Expect protein overlap between studies in
the low hundreds and missingness that is not at random.
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..io import Cohort, CohortSet


def olink_overlap(cs: CohortSet, max_below_lod_fraction: float = 0.25, lod_flag_col: Optional[str] = None) -> pd.Index:
    """UniProt assays present in every Olink cohort, excluding assays mostly below LOD.

    If a cohort stores a below-LOD indicator matrix in ``cohort.source['lod_matrix_path']``
    this function does not read it; pass precomputed fractions through ``lod_flag_col`` in
    cohort.X.attrs instead. The simple default drops assays whose within study variance is
    near zero, which is what a high below-LOD fraction produces after LOD substitution.
    """
    keep: Optional[set] = None
    for c in cs.cohorts:
        sd = c.X.std(axis=0)
        frac_const = (c.X.round(3).nunique(axis=0) <= 2).astype(float)
        ok = set(sd.index[(sd > 1e-6) & (frac_const < 1)])
        keep = ok if keep is None else keep & ok
    return pd.Index(sorted(keep or []))


def within_study_zscore(cs: CohortSet, features: Optional[Iterable[str]] = None) -> CohortSet:
    """Per study z-score (label free). Equivalent to the center_scale operator but applied
    once at preprocessing so that downstream operators act on comparable scales."""
    feats = pd.Index(features) if features is not None else cs.shared_features()
    out = []
    for c in cs.cohorts:
        X = c.X.reindex(columns=feats)
        Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1).replace(0, np.nan)
        out.append(Cohort(c.study_id, c.layer, Z.fillna(0.0), c.meta, c.feature_space, False, {**c.source, "transform": "within_study_z"}))
    return CohortSet(out)


def ms_median_normalise(cs: CohortSet, min_detect_fraction: float = 0.7, log_transform: bool = True) -> CohortSet:
    """Sample wise median centring of log2 intensities and removal of proteins quantified in
    fewer than ``min_detect_fraction`` of samples within a study. Missing values are kept as
    NaN so that a downstream model can encode missingness explicitly."""
    out = []
    for c in cs.cohorts:
        X = c.X.astype(float)
        if log_transform and (X.stack().dropna() > 100).mean() > 0.5:
            X = np.log2(X.where(X > 0))
        detect = X.notna().mean(axis=0)
        X = X.loc[:, detect >= min_detect_fraction]
        med = X.median(axis=1)
        X = X.sub(med, axis=0) + med.median()
        out.append(Cohort(c.study_id, c.layer, X, c.meta, "uniprot", False, {**c.source, "transform": "log2_median_centred"}))
    return CohortSet(out)


def missingness_indicators(X: pd.DataFrame, min_missing_fraction: float = 0.05) -> pd.DataFrame:
    """Binary indicators for proteins with meaningful missingness, to append as features."""
    frac = X.isna().mean(axis=0)
    cols = frac.index[frac >= min_missing_fraction]
    return X[cols].isna().astype(float).add_suffix("__missing")


def map_identifiers(X: pd.DataFrame, mapping: pd.DataFrame, src_col: str, dst_col: str, agg: str = "mean") -> pd.DataFrame:
    """Rename and aggregate features through an identifier mapping table.

    Used for Ensembl gene to UniProt protein (diagonal integration of transcriptomic and
    proteomic cohorts at the gene product level), Olink ID to UniProt, metabolite names to
    RefMet. Many to one mappings are aggregated with ``agg``.
    """
    m = mapping[[src_col, dst_col]].dropna().drop_duplicates()
    m = m[m[src_col].isin(X.columns)]
    renamed = X[m[src_col].values]
    renamed.columns = m[dst_col].values
    return renamed.T.groupby(level=0).agg(agg).T
