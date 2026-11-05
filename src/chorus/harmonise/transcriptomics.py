"""Bulk transcriptomics preprocessing before cross study harmonisation.

The expected input is a gene level count matrix per study, ideally produced by one uniform
pipeline from FASTQ (nf-core/rnaseq with Salmon, one GENCODE release). Author supplied
matrices are accepted but flagged in provenance because pipeline differences are themselves
a study effect (Wang et al. 2018; Arora et al. 2020).
"""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..io import Cohort, CohortSet


def filter_low_counts(cs: CohortSet, min_cpm: float = 1.0, min_fraction_samples: float = 0.5) -> pd.Index:
    """Genes with CPM >= min_cpm in at least ``min_fraction_samples`` of samples in EVERY study.

    Filtering per study and intersecting avoids keeping genes that are only expressed in one
    platform (for example globin or neutrophil genes in whole blood versus PBMC).
    """
    keep: Optional[set] = None
    for c in cs.cohorts:
        if not c.is_counts:
            genes = set(c.features)
        else:
            lib = c.X.sum(axis=1).replace(0, np.nan)
            cpm = c.X.div(lib, axis=0) * 1e6
            frac = (cpm >= min_cpm).mean(axis=0)
            genes = set(frac.index[frac >= min_fraction_samples])
        keep = genes if keep is None else keep & genes
    return pd.Index(sorted(keep or []))


def log_cpm(X: pd.DataFrame, prior_count: float = 2.0) -> pd.DataFrame:
    """log2 counts per million with a prior count, library size from all genes in X."""
    lib = X.sum(axis=1).astype(float)
    cpm = (X.astype(float) + prior_count / 2).div(lib + prior_count, axis=0) * 1e6
    return np.log2(cpm)


def to_continuous(cs: CohortSet, features: Optional[Iterable[str]] = None) -> CohortSet:
    """Return a CohortSet where count cohorts are converted to log CPM (others unchanged)."""
    feats = pd.Index(features) if features is not None else None
    out = []
    for c in cs.cohorts:
        X = c.X if feats is None else c.X.reindex(columns=feats)
        if c.is_counts:
            X = log_cpm(c.X)  # library size on all genes, then subset
            X = X if feats is None else X.reindex(columns=feats)
        out.append(Cohort(c.study_id, c.layer, X, c.meta, c.feature_space, False, {**c.source, "transform": "log2cpm" if c.is_counts else c.source.get("transform", "as_provided")}))
    return CohortSet(out)


def select_tissue(cs: CohortSet, matrices: Iterable[str] = ("whole_blood",)) -> CohortSet:
    """Keep only samples from the given sample matrices (for example whole blood only).

    Whole blood and PBMC profiles separate cleanly (Gautam et al. 2019) and neutrophil to
    lymphocyte ratio dominates whole blood variance (O'Connell 2023). The default analysis
    pools within matrix; a sensitivity analysis pools across matrices with a matrix covariate.
    """
    allowed = set(matrices)
    out = []
    for c in cs.cohorts:
        mask = c.meta["sample_matrix"].astype(str).isin(allowed)
        if mask.sum() >= 6:
            out.append(c.subset(mask))
    return CohortSet(out)


def composition_covariates(X_logcpm: pd.DataFrame, markers: Optional[dict] = None) -> pd.DataFrame:
    """Crude marker based blood composition scores (mean log CPM of marker genes).

    These are used as covariates or as protected variables in ``AddonComBat``. For a
    production analysis replace with a deconvolution estimate (CIBERSORTx, MuSiC, BayesPrism)
    computed per study. Default markers use HGNC symbols, so the matrix must be in symbol
    space or a mapping must be applied first.
    """
    markers = markers or {
        "neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "FPR1", "MNDA"],
        "lymphocyte": ["CD3E", "CD3D", "CD79A", "MS4A1", "CD2"],
        "erythroid": ["HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1"],
        "platelet": ["PPBP", "PF4", "GP9", "ITGA2B", "TUBB1"],
    }
    out = {}
    for k, genes in markers.items():
        g = [x for x in genes if x in X_logcpm.columns]
        out[k] = X_logcpm[g].mean(axis=1) if g else pd.Series(np.nan, index=X_logcpm.index)
    return pd.DataFrame(out)
