"""Programme scores: a bridge between layers measured in different cohorts.

A programme is a gene set (MSigDB Hallmark or Reactome). Its members are mapped to the feature
space of each layer: Ensembl gene ids for RNA, UniProt accessions for proteins, through the
Entrez ids of the gene sets and the HGNC complete set (``map_gene_sets``). Within each cohort
every feature is standardised over that cohort's samples (mean 0, SD 1, no labels used), and
the score of a programme in a sample is the mean of the standardised values of the members the
cohort measures. A programme is scored in a cohort only when at least ``min_members`` members
are measured there; otherwise its score is missing for the whole cohort. Because the scoring of
a cohort uses that cohort's samples alone and no labels, the scores of a held out cohort are
the same whether or not it was held out.

A member is measured in a cohort when the feature has at least one value there and a non-zero
standard deviation. A multi-protein key (Olink and SomaScan assays of complexes, "P1_P2")
belongs to a programme when any of its accessions maps to a member gene, and counts as one
member.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


def read_gmt(path: str | Path) -> Dict[str, List[str]]:
    """Gene set name -> member ids (third column onwards of a GMT file)."""
    out: Dict[str, List[str]] = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            out[parts[0]] = [g for g in parts[2:] if g]
    return out


def map_gene_sets(sets: Dict[str, List[str]], hgnc: pd.DataFrame, collection: str = "") -> pd.DataFrame:
    """Long table of programme members, one row per member gene and UniProt accession.

    ``sets`` holds Entrez ids; ``hgnc`` is the HGNC complete set (columns entrez_id, symbol,
    ensembl_gene_id, uniprot_ids). Columns: programme, collection, entrez_id, symbol,
    ensembl_gene_id, uniprot (missing where HGNC has no entry or no identifier)."""
    t = hgnc[["entrez_id", "symbol", "ensembl_gene_id", "uniprot_ids"]].copy()
    t["entrez_id"] = pd.to_numeric(t["entrez_id"], errors="coerce").astype("Int64").astype(str)
    t = t[t["entrez_id"] != "<NA>"].drop_duplicates("entrez_id")
    t["uniprot"] = t["uniprot_ids"].fillna("").astype(str).str.split("|")
    t = t.explode("uniprot")
    t["uniprot"] = t["uniprot"].str.strip().replace("", None)
    long = pd.DataFrame([(name, e) for name, genes in sets.items() for e in dict.fromkeys(genes)], columns=["programme", "entrez_id"])
    out = long.merge(t[["entrez_id", "symbol", "ensembl_gene_id", "uniprot"]], on="entrez_id", how="left")
    out.insert(1, "collection", collection)
    return out[["programme", "collection", "entrez_id", "symbol", "ensembl_gene_id", "uniprot"]].reset_index(drop=True)


def programme_members(table: pd.DataFrame, features: Iterable[str], feature_space: str = "ensembl_gene") -> Dict[str, List[str]]:
    """Programme -> the features of ``features`` that are members, in the feature space of a layer
    (``ensembl_gene`` or ``uniprot``; a multi-protein key is a member when any accession is)."""
    feats = list(dict.fromkeys(str(f) for f in features))
    col = "uniprot" if feature_space == "uniprot" else "ensembl_gene_id"
    by_id: Dict[str, set] = {}
    for p, i in zip(table["programme"], table[col]):
        if isinstance(i, str):
            by_id.setdefault(i, set()).add(p)
    out: Dict[str, List[str]] = {p: [] for p in pd.unique(table["programme"])}
    for f in feats:
        progs = set()
        for part in (f.split("_") if feature_space == "uniprot" else [f]):
            progs |= by_id.get(part, set())
        for p in progs:
            out[p].append(f)
    return out


def score_programmes(X: pd.DataFrame, studies: pd.Series, members: Dict[str, List[str]], min_members: int = 5) -> tuple:
    """Programme scores (samples x programmes) and coverage (one row per cohort and programme).

    Within each cohort of ``studies`` the features are standardised over the cohort's samples and
    a programme's score is the mean of its measured members (``min_members`` or more, else
    missing in that cohort). Coverage columns: study_id, programme, n_members (members in the
    layer's feature space), n_measured, fraction_measured, scored."""
    scores, cov = [], []
    for s in pd.unique(studies):
        idx = studies.index[studies == s]
        Xs = X.loc[idx].astype(float)
        sd = Xs.std(axis=0, ddof=1)
        measured = Xs.columns[Xs.notna().any(axis=0).values & (sd > 0).values]
        Z = (Xs[measured] - Xs[measured].mean(axis=0)) / sd[measured]
        mset = set(measured)
        S = {}
        for p, feats in members.items():
            m = [f for f in feats if f in mset]
            ok = len(m) >= min_members
            S[p] = Z[m].mean(axis=1) if ok else pd.Series(np.nan, index=idx)
            cov.append({"study_id": s, "programme": p, "n_members": len(feats), "n_measured": len(m), "fraction_measured": len(m) / len(feats) if feats else np.nan, "scored": ok})
        scores.append(pd.DataFrame(S, index=idx))
    return pd.concat(scores).reindex(X.index), pd.DataFrame(cov)


def programme_correlation(A: pd.DataFrame, B: pd.DataFrame, min_samples: int = 10) -> pd.DataFrame:
    """Per programme, Spearman and Pearson correlation across the samples shared by two score
    matrices of the same people (two layers), over the programmes scored in both."""
    from scipy import stats

    idx = A.index.intersection(B.index)
    rows = []
    for p in A.columns.intersection(B.columns):
        a, b = A.loc[idx, p], B.loc[idx, p]
        ok = a.notna() & b.notna()
        if ok.sum() < min_samples:
            continue
        rho, prho = stats.spearmanr(a[ok], b[ok])
        r, pr = stats.pearsonr(a[ok], b[ok])
        rows.append({"programme": p, "n_samples": int(ok.sum()), "spearman": rho, "spearman_p": prho, "pearson": r, "pearson_p": pr})
    return pd.DataFrame(rows)


def common_programmes(S: pd.DataFrame, studies: pd.Series, cohorts: Iterable[str]) -> List[str]:
    """Programmes with a score in every sample of every listed cohort."""
    m = studies.isin(set(cohorts))
    ok = S.loc[m].notna().all(axis=0)
    return list(S.columns[ok.values])

