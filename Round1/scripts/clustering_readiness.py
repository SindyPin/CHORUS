"""Clustering readiness: do LC cases cluster by cohort after harmonisation?

LC cases only (lc_status LC after each config's sample selection); the labels are not used otherwise. Blocks: the RNA
cohorts of configs/lc_transcriptomics_blood_v4.yaml (the genes shared by all cohorts after the expression filter) and
the affinity protein cohorts of configs/lc_proteomics_affinity_v4.yaml (the proteins shared by all, imputed as in the
config). No gene list is used (in particular not the MRCONTROL list).

Per block and label free harmoniser (identity, center_scale, quantile_ref, fsqn, rank; fitted on the LC samples of all
cohorts without labels and applied to them):
  1. features: for RNA the --top-var genes with the largest variance of the harmonised pooled matrix, for proteins all;
     standardised over the pooled samples; principal components (--n-pcs, at most samples - 1 and features)
  2. clustering of the components for k = 2 to 6: k-means (50 starts), hierarchical (Ward), consensus (--n-resample
     k-means runs on 80% subsamples of the samples; consensus matrix = times clustered together / times sampled
     together; final clusters by average linkage on 1 - consensus)
  3. adjusted Rand index (ARI) between clusters and cohort, normalised mutual information, silhouette on the
     components, smallest cluster, and for consensus the PAC (fraction of sample pairs with consensus in (0.1, 0.9))
  4. leave one cohort out: steps 1 and 2 repeated without that cohort (harmoniser, features and components refitted);
     ARI between these clusters and the clusters of the full run on the remaining samples (stability), and ARI with
     cohort among the remaining cohorts

  clusters.csv    per block, harmoniser, method and k: n_samples, n_cohorts, ari_cohort, nmi_cohort, silhouette,
                  min_cluster_size, pac (consensus only)
  stability.csv   per block, harmoniser, method, k and left out cohort: n_samples, ari_full (stability),
                  ari_cohort
  stability_summary.csv   per block, harmoniser, method and k: mean, min and max of ari_full over left out cohorts
  composition.csv per block, harmoniser, method, k, cluster and cohort: number of samples (full run)
  cohorts.csv     per block and cohort: LC samples used

  python scripts/clustering_readiness.py --out results/clustering_readiness
  python scripts/clustering_readiness.py --out results/clustering_readiness_v5 --blocks rna=configs/lc_transcriptomics_blood_v4.yaml,protein_affinity=configs/lc_proteomics_affinity_v5.yaml,protein_explore=configs/lc_proteomics_explore_v5.yaml

--blocks replaces the blocks (name=config, comma separated); the --top-var gene filter applies to the block named rna.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

from chorus.cli import load_config
from chorus.harmonise import make_harmoniser

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = {"rna": "configs/lc_transcriptomics_blood_v4.yaml", "protein_affinity": "configs/lc_proteomics_affinity_v4.yaml"}
LABEL_FREE = ["identity", "center_scale", "quantile_ref", "fsqn", "rank"]
KS = range(2, 7)
METHODS = ["kmeans", "hierarchical", "consensus"]


def lc_cases(config: str) -> tuple:
    cs, lc, cfg = load_config(str(ROOT / config))
    X, M = cs.pooled(cfg["filtered_features"])
    keep = M["lc_status"].astype(str) == lc.positive
    return X[keep], M.loc[keep, "study_id"].astype(str)


def embed(X: pd.DataFrame, s: pd.Series, harmoniser: str, top_var: int, n_pcs: int) -> np.ndarray:
    H = make_harmoniser(harmoniser).fit(X, s)
    Z = H.transform(X, s)
    Z = Z.loc[:, Z.notna().all(axis=0)]
    if top_var and Z.shape[1] > top_var:
        Z = Z[Z.var(axis=0).sort_values(ascending=False).index[:top_var]]
    sd = Z.std(axis=0, ddof=1)
    Z = Z.loc[:, sd > 0]
    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
    n = min(n_pcs, Z.shape[0] - 1, Z.shape[1])
    return PCA(n_components=n, random_state=0).fit_transform(Z.values)


def consensus(P: np.ndarray, k: int, n_resample: int, frac: float, seed: int) -> tuple:
    rng = np.random.default_rng(seed)
    n = P.shape[0]
    together, sampled = np.zeros((n, n)), np.zeros((n, n))
    for b in range(n_resample):
        idx = np.sort(rng.choice(n, int(round(frac * n)), replace=False))
        lab = KMeans(n_clusters=k, n_init=10, random_state=int(rng.integers(1 << 31))).fit_predict(P[idx])
        same = (lab[:, None] == lab[None, :]).astype(float)
        sampled[np.ix_(idx, idx)] += 1
        together[np.ix_(idx, idx)] += same
    C = np.divide(together, sampled, out=np.zeros_like(together), where=sampled > 0)
    np.fill_diagonal(C, 1.0)
    D = squareform(1 - C, checks=False)
    labels = fcluster(linkage(D, method="average"), t=k, criterion="maxclust") - 1
    iu = np.triu_indices(n, 1)
    pac = float(((C[iu] > 0.1) & (C[iu] < 0.9)).mean())
    return labels, pac


def cluster_all(P: np.ndarray, n_resample: int, seed: int) -> dict:
    """(method, k) -> (labels, pac or nan)."""
    out = {}
    for k in KS:
        out[("kmeans", k)] = (KMeans(n_clusters=k, n_init=50, random_state=seed).fit_predict(P), np.nan)
        out[("hierarchical", k)] = (AgglomerativeClustering(n_clusters=k, linkage="ward").fit_predict(P), np.nan)
        out[("consensus", k)] = consensus(P, k, n_resample, 0.8, seed + k)
    return out


def run_block(block: str, config: str, top_var: int, n_pcs: int, n_resample: int, seed: int) -> tuple:
    X, s = lc_cases(config)
    cohorts = s.value_counts().reindex(pd.unique(s))
    print(f"{block}: {X.shape[0]} LC samples, {X.shape[1]} features, " + ", ".join(f"{c} {n}" for c, n in cohorts.items()), flush=True)
    tv = top_var if block == "rna" else 0
    rows, stab, comp = [], [], []
    for h in LABEL_FREE:
        t0 = time.time()
        P = embed(X, s, h, tv, n_pcs)
        full = cluster_all(P, n_resample, seed)
        for (method, k), (lab, pac) in full.items():
            rows.append(
                {
                    "block": block, "harmoniser": h, "method": method, "k": k, "n_samples": len(lab), "n_cohorts": s.nunique(), "n_pcs": P.shape[1],
                    "ari_cohort": adjusted_rand_score(s, lab), "nmi_cohort": normalized_mutual_info_score(s, lab),
                    "silhouette": silhouette_score(P, lab) if len(set(lab)) > 1 else np.nan,
                    "min_cluster_size": int(np.bincount(lab).min()), "pac": pac,
                }
            )
            ct = pd.crosstab(pd.Series(lab, name="cluster"), pd.Series(s.values, name="cohort")).stack().rename("n").reset_index()
            comp.append(ct.assign(block=block, harmoniser=h, method=method, k=k))
        for left in cohorts.index:
            keep = (s != left).to_numpy()
            Pl = embed(X[keep], s[keep], h, tv, n_pcs)
            sub = cluster_all(Pl, n_resample, seed)
            for (method, k), (lab, _) in sub.items():
                stab.append(
                    {
                        "block": block, "harmoniser": h, "method": method, "k": k, "left_out": left, "n_samples": int(keep.sum()),
                        "ari_full": adjusted_rand_score(full[(method, k)][0][keep], lab), "ari_cohort": adjusted_rand_score(s[keep], lab),
                    }
                )
        print(f"  {h}: {time.time() - t0:.0f} s", flush=True)
    counts = cohorts.rename("n_lc").rename_axis("study_id").reset_index().assign(block=block)
    return rows, stab, comp, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-var", type=int, default=2000, help="RNA genes kept by variance after harmonisation")
    ap.add_argument("--n-pcs", type=int, default=20)
    ap.add_argument("--n-resample", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--blocks", help="name=config,... (default: the version 4 blocks)")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    blocks = dict(b.split("=", 1) for b in a.blocks.split(",")) if a.blocks else BLOCKS
    R, S, C, N = [], [], [], []
    for block, config in blocks.items():
        r, st, cp, n = run_block(block, config, a.top_var, a.n_pcs, a.n_resample, a.seed)
        R += r
        S += st
        C += cp
        N.append(n)
    R, S = pd.DataFrame(R), pd.DataFrame(S)
    R.to_csv(out / "clusters.csv", index=False)
    S.to_csv(out / "stability.csv", index=False)
    keys = ["block", "harmoniser", "method", "k"]
    S.groupby(keys, sort=False)["ari_full"].agg(ari_full_mean="mean", ari_full_min="min", ari_full_max="max", n_left_out="count").reset_index().to_csv(out / "stability_summary.csv", index=False)
    pd.concat(C, ignore_index=True)[[*keys, "cluster", "cohort", "n"]].to_csv(out / "composition.csv", index=False)
    pd.concat(N, ignore_index=True)[["block", "study_id", "n_lc"]].to_csv(out / "cohorts.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"blocks": blocks, "harmonisers": LABEL_FREE, "k": list(KS), "methods": METHODS, "top_var_rna": a.top_var, "n_pcs": a.n_pcs, "n_resample": a.n_resample, "subsample_fraction": 0.8, "seed": a.seed}, indent=1))
    print(R.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
