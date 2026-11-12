"""Sex inferred from expression, used as a positive control label and as a biological
conservation check of harmonisation.

Within each cohort a two component Gaussian mixture is fitted to log2 CPM of XIST (female,
X inactivation) and RPS4Y1 (male, chrY). The component with the higher mean XIST minus
RPS4Y1 is called female. A sample is ``ambiguous`` when its posterior for the assigned
component is below ``min_posterior`` or when it lies outside that component (squared
Mahalanobis distance above the chi-square quantile ``1 - outlier_p``), for example a sample
with neither XIST nor RPS4Y1 or with both. When a matrix has no XIST row the split uses
RPS4Y1 alone (the component with the higher RPS4Y1 mean is male) and the method string says
so. Inference is per cohort because library chemistry shifts both genes between cohorts.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd
from scipy import stats

XIST = "ENSG00000229807"
RPS4Y1 = "ENSG00000129824"


def infer_sex(
    logcpm: pd.DataFrame,
    min_posterior: float = 0.95,
    outlier_p: float = 0.001,
    reg_covar: float = 0.25,
    seed: int = 0,
) -> pd.DataFrame:
    """Infer sex for the samples (rows) of one cohort from log2 CPM of XIST and RPS4Y1.

    Returns one row per sample: inferred_sex (female, male or ambiguous), posterior,
    mahalanobis2 (to the assigned component), reason (empty, low_posterior or outlier),
    xist_logcpm, rps4y1_logcpm and inferred_sex_method. ``reg_covar`` (log2 units squared)
    keeps a component whose gene sits at the detection floor from collapsing to zero variance.
    """
    from sklearn.mixture import GaussianMixture

    if RPS4Y1 not in logcpm.columns:
        raise ValueError("RPS4Y1 is not in the matrix; sex cannot be inferred")
    genes = [g for g in (XIST, RPS4Y1) if g in logcpm.columns]
    Z = logcpm[genes].astype(float)
    if Z.isna().any().any():
        raise ValueError("missing XIST or RPS4Y1 values")
    gm = GaussianMixture(n_components=2, covariance_type="full", reg_covar=reg_covar, n_init=5, random_state=seed).fit(Z.values)
    post = gm.predict_proba(Z.values)
    if XIST in genes:
        score = gm.means_[:, 0] - gm.means_[:, 1]  # XIST minus RPS4Y1
        method = "gmm2:XIST+RPS4Y1"
    else:
        score = -gm.means_[:, 0]  # RPS4Y1 only: lower mean is female
        method = "gmm2:RPS4Y1 (no XIST row in the matrix)"
    female = int(np.argmax(score))
    comp = post.argmax(axis=1)
    d2 = np.array([_mahalanobis2(Z.values[i], gm.means_[c], gm.covariances_[c]) for i, c in enumerate(comp)])
    cut = stats.chi2.ppf(1 - outlier_p, df=len(genes))
    pmax = post.max(axis=1)
    reason = np.where(pmax < min_posterior, "low_posterior", np.where(d2 > cut, "outlier", ""))
    sex = np.where(comp == female, "female", "male")
    sex = np.where(reason == "", sex, "ambiguous")
    return pd.DataFrame(
        {
            "inferred_sex": sex,
            "posterior": pmax,
            "mahalanobis2": d2,
            "reason": reason,
            "xist_logcpm": logcpm[XIST].astype(float) if XIST in genes else np.nan,
            "rps4y1_logcpm": logcpm[RPS4Y1].astype(float),
            "inferred_sex_method": method,
        },
        index=logcpm.index,
    )


def _mahalanobis2(x: np.ndarray, mu: np.ndarray, cov: np.ndarray) -> float:
    d = x - mu
    return float(d @ np.linalg.solve(cov, d))


def gene_chromosomes(table: str | Path) -> pd.Series:
    """Chromosome of each Ensembl gene (version removed) from a GENCODE GTF or from a tx2gene
    table with gene_id and chrom columns (scripts/salmon_quant_ena.py writes one). A gene id
    that appears on two chromosomes (the chrX and chrY pseudoautosomal copies) keeps both,
    joined by a comma."""
    table = Path(table)
    if table.name.endswith((".gtf", ".gtf.gz")):
        rows = []
        opener = gzip.open if table.suffix == ".gz" else open
        with opener(table, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if f[2] == "gene":
                    rows.append((re.search(r'gene_id "([^"]+)"', f[8]).group(1), f[0]))
        t = pd.DataFrame(rows, columns=["gene_id", "chrom"])
    else:
        t = pd.read_csv(table, sep="\t", usecols=["gene_id", "chrom"])
    t["gene_id"] = t["gene_id"].str.split(".").str[0]
    return t.drop_duplicates().groupby("gene_id")["chrom"].agg(lambda c: ",".join(sorted(set(c))))


def sex_chromosome_genes(chroms: pd.Series, features: Optional[Iterable[str]] = None) -> List[str]:
    """Genes on chrX or chrY (including the pseudoautosomal genes), optionally within ``features``."""
    on_xy = chroms.index[chroms.str.contains(r"\bchr[XY]\b", regex=True)]
    if features is None:
        return sorted(on_xy)
    keep = set(on_xy)
    return [f for f in features if f in keep]
