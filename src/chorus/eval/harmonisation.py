"""Harmonisation diagnostics on the leave-one-study-out training pool.

For every fold of the benchmark (held-out study S) the training pool is the labelled samples
of the other studies, exactly as ``run_loso`` builds it. Each harmoniser is fitted as in the
benchmark (``fit_harmoniser``) and the training pool is compared before (log CPM as loaded)
and after harmonisation with two families of metrics.

Technical mixing, on the filtered genes of the config:
  silhouette_study        mean silhouette width with study as the cluster label, Euclidean,
                          on the top ``n_pcs`` principal components of the pool after each
                          gene is standardised across the pool
  knn_other_study         mean over samples of the fraction of the ``k`` nearest neighbours
                          (same PC space) that come from another study; ``expected_random``
                          is that fraction under random mixing, mean of (N - n_s) / (N - 1)

Biological conservation:
  lc_effect_spearman      per training study, Spearman correlation across the filtered genes
                          of the within-study LC effect (Hedges g, LC minus COVID_recovered)
                          before versus after harmonisation; also their mean over studies
  sex_auroc_XIST          AUROC of XIST for inferred female versus male over the pool
  sex_auroc_RPS4Y1        AUROC of RPS4Y1 for inferred male versus female over the pool
  silhouette_sex          mean silhouette width by inferred sex on XIST and RPS4Y1 (each
                          standardised across the pool), on those of the two genes available
  marker_score_spearman   neutrophil and lymphocyte scores (mean over the marker genes of the
                          values standardised across the pool), Spearman correlation of the
                          score before versus after, over the pool and within each study (mean)

XIST, RPS4Y1 and marker genes that the expression filter removes are harmonised together
with the filtered genes (``load_config(extra_features=...)``); a gene enters a fold only when
every sample the harmoniser is fitted on has a value (GSE270045 has no XIST row), so the
mixing and LC effect metrics are computed on the filtered genes only and the sex and marker
metrics report how many genes they used.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, silhouette_samples
from sklearn.neighbors import NearestNeighbors

from ..features.invariant import per_study_effects
from ..io import CohortSet
from ..schema import binary_label
from ..sex import RPS4Y1, XIST
from .loso import LosoConfig, _covariate_matrix, fit_harmoniser

# Ensembl ids of the marker genes in chorus.harmonise.transcriptomics.composition_covariates
MARKERS: Dict[str, Dict[str, str]] = {
    "neutrophil": {"FCGR3B": "ENSG00000162747", "CSF3R": "ENSG00000119535", "CXCR2": "ENSG00000180871", "FPR1": "ENSG00000171051", "MNDA": "ENSG00000163563"},
    "lymphocyte": {"CD3E": "ENSG00000198851", "CD3D": "ENSG00000167286", "CD79A": "ENSG00000105369", "MS4A1": "ENSG00000156738", "CD2": "ENSG00000116824"},
}
SEX_GENES = {"XIST": XIST, "RPS4Y1": RPS4Y1}


def diagnostic_genes() -> List[str]:
    """Genes the sex and marker metrics need, to pass as ``load_config(extra_features=...)``."""
    return list(SEX_GENES.values()) + [g for m in MARKERS.values() for g in m.values()]


def _standardise(X: pd.DataFrame) -> pd.DataFrame:
    sd = X.std(axis=0, ddof=1)
    return (X - X.mean(axis=0)) / sd.where(sd > 0, 1.0)


def pool_pcs(X: pd.DataFrame, n_pcs: int = 50, seed: int = 0) -> np.ndarray:
    """Top principal components of the pool after per gene standardisation across the pool."""
    Z = _standardise(X).values
    n = min(n_pcs, Z.shape[0] - 1, Z.shape[1])
    return PCA(n_components=n, random_state=seed).fit_transform(Z)


def mixing_metrics(X: pd.DataFrame, studies: pd.Series, n_pcs: int = 50, k: int = 20, seed: int = 0) -> pd.DataFrame:
    """Silhouette by study and the other-study fraction of the k nearest neighbours, pooled and per study."""
    P = pool_pcs(X, n_pcs, seed)
    s = np.asarray(studies.loc[X.index].astype(str), dtype=object)
    sil = silhouette_samples(P, s)
    kk = min(k, len(s) - 1)
    idx = NearestNeighbors(n_neighbors=kk + 1).fit(P).kneighbors(P, return_distance=False)
    # the first neighbour of a point is the point itself unless there are exact duplicates
    nbr = np.array([[j for j in row if j != i][:kk] for i, row in enumerate(idx)])
    other = (s[nbr] != s[:, None]).mean(axis=1)
    sizes = pd.Series(s).value_counts()
    expected = (len(s) - sizes.reindex(s).values) / (len(s) - 1)
    rows = []
    for scope, m in [("pool", np.ones(len(s), bool))] + [(st, s == st) for st in sorted(set(s))]:
        rows += [
            {"metric": "silhouette_study", "scope": scope, "value": float(sil[m].mean()), "n": int(m.sum())},
            {"metric": "knn_other_study", "scope": scope, "value": float(other[m].mean()), "n": int(m.sum())},
            {"metric": "knn_other_study_expected_random", "scope": scope, "value": float(expected[m].mean()), "n": int(m.sum())},
        ]
    return pd.DataFrame(rows).assign(n_pcs=P.shape[1], k=kk)


def lc_effect_conservation(X_before: pd.DataFrame, X_after: pd.DataFrame, y: pd.Series, studies: pd.Series) -> pd.DataFrame:
    """Spearman correlation across genes of the within-study LC effect before versus after."""
    eb = per_study_effects(X_before, y, studies)
    ea = per_study_effects(X_after.reindex(columns=X_before.columns), y, studies)
    rows = []
    for st in sorted(eb):
        b, a = eb[st]["g"], ea[st]["g"]
        ok = b.notna() & a.notna()
        rho = stats.spearmanr(b[ok], a[ok])[0] if ok.sum() > 2 else np.nan
        rows.append({"metric": "lc_effect_spearman", "scope": st, "before_vs_after": float(rho), "n": int(ok.sum())})
    out = pd.DataFrame(rows)
    if len(out):
        out = pd.concat([out, pd.DataFrame([{"metric": "lc_effect_spearman", "scope": "mean_over_studies", "before_vs_after": float(out["before_vs_after"].mean()), "n": int(len(out))}])], ignore_index=True)
    return out


def sex_separation(X: pd.DataFrame, sex: pd.Series) -> pd.DataFrame:
    """AUROC of XIST (female) and RPS4Y1 (male) and the silhouette by inferred sex on those genes."""
    sx = sex.reindex(X.index)
    m = sx.isin(["female", "male"])
    female = (sx[m] == "female").astype(int).values
    rows = []
    genes = [g for g in SEX_GENES.values() if g in X.columns]
    for name, g, pos in (("sex_auroc_XIST", XIST, 1), ("sex_auroc_RPS4Y1", RPS4Y1, 0)):
        ok = g in X.columns and 0 < female.sum() < len(female)
        val = roc_auc_score(female == pos, X.loc[m[m].index, g]) if ok else np.nan
        rows.append({"metric": name, "scope": "pool", "value": float(val), "n": int(m.sum()), "n_genes": int(g in X.columns)})
    if genes and 0 < female.sum() < len(female):
        Z = _standardise(X.loc[m[m].index, genes]).values
        sil = float(silhouette_samples(Z, female).mean())
    else:
        sil = np.nan
    rows.append({"metric": "silhouette_sex", "scope": "pool", "value": sil, "n": int(m.sum()), "n_genes": len(genes)})
    return pd.DataFrame(rows)


def marker_scores(X: pd.DataFrame, markers: Optional[Dict[str, Dict[str, str]]] = None) -> pd.DataFrame:
    """Mean of the marker genes standardised across the pool, one column per cell type."""
    markers = markers or MARKERS
    Z = _standardise(X)
    out = {}
    for cell, genes in markers.items():
        g = [e for e in genes.values() if e in Z.columns]
        out[cell] = Z[g].mean(axis=1) if g else pd.Series(np.nan, index=X.index)
    return pd.DataFrame(out)


def marker_preservation(X_before: pd.DataFrame, X_after: pd.DataFrame, studies: pd.Series, markers: Optional[Dict[str, Dict[str, str]]] = None) -> pd.DataFrame:
    markers = markers or MARKERS
    rows = []
    for cell, genes in markers.items():
        g = {k: e for k, e in genes.items() if e in X_before.columns and e in X_after.columns}
        if not g:
            rows.append({"metric": f"marker_score_spearman_{cell}", "scope": "pool", "before_vs_after": np.nan, "n": 0, "n_genes": 0})
            continue
        b = marker_scores(X_before, {cell: g})[cell]
        a = marker_scores(X_after, {cell: g})[cell]
        rows.append({"metric": f"marker_score_spearman_{cell}", "scope": "pool", "before_vs_after": float(stats.spearmanr(b, a)[0]), "n": int(len(b)), "n_genes": len(g)})
        within = []
        s = studies.loc[b.index].astype(str)
        for st in sorted(s.unique()):
            idx = s.index[s == st]
            r = float(stats.spearmanr(b[idx], a[idx])[0]) if len(idx) > 2 else np.nan
            within.append(r)
            rows.append({"metric": f"marker_score_spearman_{cell}", "scope": st, "before_vs_after": r, "n": int(len(idx)), "n_genes": len(g)})
        rows.append({"metric": f"marker_score_spearman_{cell}", "scope": "mean_over_studies", "before_vs_after": float(np.nanmean(within)), "n": int(len(within)), "n_genes": len(g)})
    return pd.DataFrame(rows)


def harmonisation_table(
    X_before: pd.DataFrame,
    X_after: pd.DataFrame,
    studies: pd.Series,
    y: pd.Series,
    features: Sequence[str],
    sex: Optional[pd.Series] = None,
    n_pcs: int = 50,
    k: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """All metrics for one harmonised pool. Columns: family, metric, scope, before, after,
    before_vs_after, n, n_genes, n_pcs, k. Mixing and sex separation fill before and after;
    the LC effect and marker comparisons fill before_vs_after."""
    feats = list(features)
    mb = mixing_metrics(X_before[feats], studies, n_pcs, k, seed)
    ma = mixing_metrics(X_after[feats], studies, n_pcs, k, seed)
    mix = mb.rename(columns={"value": "before"}).merge(ma[["metric", "scope", "value"]].rename(columns={"value": "after"}), on=["metric", "scope"])
    mix["n_genes"] = len(feats)
    parts = [mix.assign(family="technical_mixing")]
    parts.append(lc_effect_conservation(X_before[feats], X_after[feats], y, studies).assign(family="biological_conservation", n_genes=len(feats)))
    if sex is not None:
        sb, sa = sex_separation(X_before, sex), sex_separation(X_after, sex)
        parts.append(sb.rename(columns={"value": "before"}).merge(sa[["metric", "value"]].rename(columns={"value": "after"}), on="metric").assign(family="biological_conservation"))
    parts.append(marker_preservation(X_before, X_after, studies).assign(family="biological_conservation"))
    cols = ["family", "metric", "scope", "before", "after", "before_vs_after", "n", "n_genes", "n_pcs", "k"]
    return pd.concat(parts, ignore_index=True).reindex(columns=cols)


def run_harmonisation_diagnostics(
    cs: CohortSet,
    cfg: LosoConfig,
    features: Sequence[str],
    extra: Iterable[str] = (),
    harmonisers: Optional[List[str]] = None,
    sex_col: str = "inferred_sex",
    out: Optional[Path] = None,
    n_pcs: int = 50,
    k: int = 20,
    verbose: bool = True,
) -> pd.DataFrame:
    """Metrics for every LOSO fold and harmoniser. With ``out``, one table per fold and
    harmoniser is written to ``out/<held_out>/<harmoniser>.csv`` and the long table of all of
    them to ``out/all_metrics.csv``."""
    feats = list(features)
    extra = [g for g in extra if g not in set(feats)]
    X_all, M_all = cs.pooled(feats + extra)
    y_all = binary_label(M_all, cfg.positive, cfg.negative)
    s_all = M_all["study_id"].astype(str)
    lab = y_all.notna()
    X_all, M_all, y_all, s_all = X_all[lab], M_all[lab], y_all[lab], s_all[lab]
    prot_all = _covariate_matrix(M_all, cfg.protect) if cfg.protect else None
    sex_all = M_all[sex_col] if sex_col in M_all.columns else None
    tables = []
    for held_out in cs.study_ids:
        if cs[held_out].source.get("train_only") == "true":
            continue
        te = s_all == held_out
        tr = ~te
        if (y_all[te] == 1).sum() < cfg.min_test_per_class or (y_all[te] == 0).sum() < cfg.min_test_per_class:
            continue
        if s_all[tr].nunique() < cfg.min_train_studies:
            continue
        protect_tr = None if prot_all is None else prot_all.loc[tr[tr].index]
        for h in harmonisers or cfg.harmonisers:
            # extra genes enter only when every sample the harmoniser is fitted on has a value
            fit_rows = X_all if h == "combat_pooled" else X_all[tr]
            cols = feats + [g for g in extra if fit_rows[g].notna().all()]
            X_tr, X_te = X_all.loc[tr, cols], X_all.loc[te, cols]
            try:
                H, leak = fit_harmoniser(h, X_tr, s_all[tr], y_all[tr], X_te, s_all[te], protect_tr)
                X_h = H.transform(X_tr, s_all[tr])
                t = harmonisation_table(X_tr, X_h, s_all[tr], y_all[tr], feats, sex_all, n_pcs, k)
            except Exception as e:  # keep going, record the failure
                t = pd.DataFrame([{"family": "error", "metric": "error", "scope": str(e)[:200]}])
                leak = np.nan
            t.insert(0, "harmoniser", h)
            t.insert(0, "held_out", held_out)
            t["leak_risk"] = leak
            t["extra_genes"] = ",".join(g for g in cols if g not in set(feats))
            if out is not None:
                d = Path(out) / held_out
                d.mkdir(parents=True, exist_ok=True)
                t.to_csv(d / f"{h}.csv", index=False)
            tables.append(t)
            if verbose:
                sil = t[(t["metric"] == "silhouette_study") & (t["scope"] == "pool")]
                msg = f"silhouette by study {sil['before'].iloc[0]:.3f} -> {sil['after'].iloc[0]:.3f}" if len(sil) else t["scope"].iloc[0]
                print(f"[{held_out}] {h:>13} | {msg}", flush=True)
    res = pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()
    if out is not None:
        res.to_csv(Path(out) / "all_metrics.csv", index=False)
    return res
