"""Noise ceiling of between-cohort concordance: how well can two cohorts agree, given their sizes?

For every cohort of a concordance config (each cohort on its own features, the samples of the config):

  split half     --n-splits random stratified half splits (each class split in two halves of floor(n/2) and the rest);
                 Hedges g of LC minus COVID_recovered in each half; Spearman between the halves. Computed on the cohort's
                 own features (cohorts.csv) and, for every pair, on the features the pair shares (the features of the
                 between-cohort correlation). Reliability: Spearman-Brown 2r / (1 + r) of the mean split half r.
  pairs          observed Spearman of the two cohorts' g over their shared features (as scripts/concordance.py);
                 disattenuated = observed / sqrt(reliability_a * reliability_b) on the shared features (NaN when a
                 reliability is not positive); expected if both cohorts came from one population: the larger cohort
                 (more samples) is split, --n-splits times, into a stratified subsample with the class counts of the
                 smaller cohort (capped at half of each class) and the remaining samples, and the Spearman of the two
                 parts' g over the shared features is taken; its mean and 2.5th and 97.5th percentiles are reported,
                 and below_2p5 marks an observed correlation under the 2.5th percentile.
  summary        per layer: pairs, pairs below the 2.5th percentile, and whether more than half are (the bar).

Hedges g needs at least 2 samples per class in a part; a part with fewer gives NaN.

  python scripts/concordance_noise_ceiling.py --out results/concordance_noise_ceiling
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from chorus.cli import load_config
from chorus.features.invariant import _hedges_g
from chorus.schema import binary_label

ROOT = Path(__file__).resolve().parents[1]
LAYERS = {
    "transcriptomics": "configs/concordance_rna_v5.yaml",
    "proteomics": "configs/concordance_proteomics_all_v5.yaml",
    "metabolomics": "configs/concordance_metabolomics_all.yaml",
}


def cohorts_of(config: str) -> dict:
    cs, lc, cfg = load_config(str(ROOT / config))
    out = {}
    for c in cs.cohorts:
        X = c.X.loc[:, c.X.notna().any(axis=0)].astype(float)
        y = binary_label(c.meta.loc[X.index], lc.positive, lc.negative)
        keep = y.notna().to_numpy()
        out[c.study_id] = (X[keep], y[keep].astype(int).to_numpy())
    return out


def g_of(X: np.ndarray, y: np.ndarray, idx: np.ndarray) -> np.ndarray:
    yy = y[idx]
    if min((yy == 1).sum(), (yy == 0).sum()) < 2:
        return np.full(X.shape[1], np.nan)
    with np.errstate(all="ignore"):
        return _hedges_g(X[idx][yy == 1], X[idx][yy == 0])[0]


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    return float(stats.spearmanr(a[ok], b[ok]).statistic)


def half_split(y: np.ndarray, rng) -> tuple:
    a, b = [], []
    for k in (0, 1):
        idx = rng.permutation(np.where(y == k)[0])
        h = len(idx) // 2
        a.append(idx[:h])
        b.append(idx[h:])
    return np.concatenate(a), np.concatenate(b)


def sub_split(y: np.ndarray, counts: dict, rng) -> tuple:
    """A stratified subsample with ``counts`` per class (capped at half of the class) and the rest."""
    a, b = [], []
    for k in (0, 1):
        idx = rng.permutation(np.where(y == k)[0])
        m = min(counts[k], len(idx) // 2)
        a.append(idx[:m])
        b.append(idx[m:])
    return np.concatenate(a), np.concatenate(b)


def sb(r: float) -> float:
    return 2 * r / (1 + r) if np.isfinite(r) and r > -1 else np.nan


def run_layer(layer: str, config: str, n_splits: int, seed: int) -> tuple:
    C = cohorts_of(config)
    names = list(C)
    cols = {n: C[n][0].columns for n in names}
    # split-half g of every cohort, all features, kept to reuse on each pair's shared features
    halves = {}
    coh_rows = []
    for i, n in enumerate(names):
        X, y = C[n][0].to_numpy(), C[n][1]
        rng = np.random.default_rng([seed, 0, i])
        H = [(g_of(X, y, a), g_of(X, y, b)) for a, b in (half_split(y, rng) for _ in range(n_splits))]
        halves[n] = H
        r = np.array([spearman(ga, gb) for ga, gb in H])
        coh_rows.append({"layer": layer, "cohort": n, "n_LC": int((y == 1).sum()), "n_COVID_recovered": int((y == 0).sum()), "n_features": X.shape[1],
                         "split_half_r_mean": np.nanmean(r) if np.isfinite(r).any() else np.nan, "split_half_r_sd": np.nanstd(r, ddof=1) if np.isfinite(r).sum() > 1 else np.nan,
                         "n_valid_splits": int(np.isfinite(r).sum()), "reliability_spearman_brown": sb(np.nanmean(r)) if np.isfinite(r).any() else np.nan})
        print(f"{layer} {n}: split half {coh_rows[-1]['split_half_r_mean']:.3f}", flush=True)
    pair_rows = []
    for (i, a), (j, b) in itertools.combinations(enumerate(names), 2):
        shared = cols[a].intersection(cols[b])
        Xa, ya = C[a][0][shared].to_numpy(), C[a][1]
        Xb, yb = C[b][0][shared].to_numpy(), C[b][1]
        row = {"layer": layer, "cohort_a": a, "cohort_b": b, "n_shared": len(shared), "n_a": len(ya), "n_b": len(yb)}
        if len(shared) < 3:
            pair_rows.append({**row, "observed": np.nan})
            continue
        obs = spearman(g_of(Xa, ya, np.arange(len(ya))), g_of(Xb, yb, np.arange(len(yb))))
        pos_a = cols[a].get_indexer(shared)
        pos_b = cols[b].get_indexer(shared)
        rel = {}
        for n, pos in ((a, pos_a), (b, pos_b)):
            r = np.array([spearman(ga[pos], gb[pos]) for ga, gb in halves[n]])
            rel[n] = (np.nanmean(r) if np.isfinite(r).any() else np.nan)
        rel_a, rel_b = sb(rel[a]), sb(rel[b])
        dis = obs / np.sqrt(rel_a * rel_b) if (np.isfinite(rel_a) and np.isfinite(rel_b) and rel_a > 0 and rel_b > 0) else np.nan
        big, small = (a, b) if len(ya) >= len(yb) else (b, a)
        Xg, yg = (Xa, ya) if big == a else (Xb, yb)
        ys = yb if small == b else ya
        counts = {0: int((ys == 0).sum()), 1: int((ys == 1).sum())}
        rng = np.random.default_rng([seed, 1, i, j])
        exp = np.array([spearman(g_of(Xg, yg, p), g_of(Xg, yg, q)) for p, q in (sub_split(yg, counts, rng) for _ in range(n_splits))])
        ok = exp[np.isfinite(exp)]
        lo = np.percentile(ok, 2.5) if len(ok) else np.nan
        pair_rows.append({**row, "observed": obs, "split_half_a_shared": rel[a], "split_half_b_shared": rel[b], "reliability_a_shared": rel_a, "reliability_b_shared": rel_b,
                          "disattenuated": dis, "larger_cohort": big, "expected_mean": ok.mean() if len(ok) else np.nan, "expected_p2_5": lo,
                          "expected_p97_5": np.percentile(ok, 97.5) if len(ok) else np.nan, "n_expected_draws": len(ok),
                          "below_2p5": bool(obs < lo) if np.isfinite(obs) and np.isfinite(lo) else np.nan})
    return pd.DataFrame(coh_rows), pd.DataFrame(pair_rows)


VERDICTS = ("agree", "truly differ", "too small to tell", "not evaluable")


def verdict(reliability_a: float, reliability_b: float, observed: float, expected_p2_5: float, threshold: float = 0.3) -> str:
    """One verdict per cohort pair (thresholds fixed before the validation of 18 November 2026):
    "too small to tell" when either cohort's Spearman-Brown reliability on the shared features is at or below
    ``threshold`` (or undefined); otherwise "agree" when the observed correlation is at or above the 2.5th percentile
    of the one population expectation, else "truly differ". "not evaluable" when the observed correlation or the
    expectation is undefined (fewer than 3 shared features, or no valid draw)."""
    if not (np.isfinite(observed) and np.isfinite(expected_p2_5)):
        return "not evaluable"
    if not (np.isfinite(reliability_a) and np.isfinite(reliability_b)) or min(reliability_a, reliability_b) <= threshold:
        return "too small to tell"
    return "agree" if observed >= expected_p2_5 else "truly differ"


def summary(P: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for layer, d in P.groupby("layer", sort=False):
        v = d["below_2p5"].dropna().astype(bool)
        rows.append({"layer": layer, "n_pairs": len(d), "n_pairs_evaluable": len(v), "n_below_2p5": int(v.sum()),
                     "fraction_below": v.mean() if len(v) else np.nan, "disagreement_real_by_bar": bool(v.sum() > len(v) / 2) if len(v) else None,
                     "median_observed": d["observed"].median(), "median_expected": d["expected_mean"].median(), "median_disattenuated": d["disattenuated"].median()})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-splits", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--layers", default=",".join(LAYERS))
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    Cs, Ps = [], []
    for layer in a.layers.split(","):
        c, p = run_layer(layer, LAYERS[layer], a.n_splits, a.seed)
        Cs.append(c)
        Ps.append(p)
    C, P = pd.concat(Cs, ignore_index=True), pd.concat(Ps, ignore_index=True)
    C.to_csv(out / "cohorts.csv", index=False)
    P.to_csv(out / "pairs.csv", index=False)
    S = summary(P)
    S.to_csv(out / "summary.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"layers": {l: LAYERS[l] for l in a.layers.split(",")}, "n_splits": a.n_splits, "seed": a.seed}, indent=1))
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
