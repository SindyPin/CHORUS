"""Transfer as a function of reliability and agreement: a synthesis of existing outputs (results/transfer_model).

1  stratified.csv     post hoc: from results/concordance_noise_ceiling/pairs.csv, per layer and reliability threshold
                      t (0.2, 0.3, 0.4, 0.5), the fraction of pairs below the one population 2.5th percentile among the
                      pairs where both cohorts have Spearman-Brown split half reliability above t on the shared features,
                      and among the rest
2  reliability.csv    split half reliability per held-out cohort and contrast, with the procedure of analysis A
                      (scripts/concordance_noise_ceiling.py: 200 stratified half splits, Hedges g in each half, Spearman
                      between halves, Spearman-Brown), on the features of the leave one study out run each point comes from:
                        LC          the LC blocks of concordance_vs_transfer_full (transcriptomics v4, affinity v5, MS v4,
                                    metabolomics), LC versus COVID_recovered
                        sex         configs/lc_transcriptomics_blood.yaml as in control_sex_label (all genes): inferred
                                    female versus male among the samples of the LC versus COVID_recovered contrast
                        acute       configs/control_acute_proteomics.yaml: acute versus never infected
                      plus mean_abs_g, the mean absolute Hedges g of the contrast over those features on the whole cohort
3  models.csv         to each LC and control point of results/concordance_vs_transfer_full/points.csv (no_xy points
                      excluded) the held-out cohort's reliability and mean_abs_g for its contrast are added. Ordinary least
                      squares of the held-out AUROC on: agreement (effect_spearman); reliability; agreement + reliability;
                      agreement * reliability; and the same with mean_abs_g in place of reliability. R2 and coefficients,
                      with 2,000 resample cohort bootstrap intervals (held-out cohorts resampled within contrast).
   gap.csv            for each model, the LC minus control gap: the coefficient of an LC indicator added to the model, with
                      a permutation p (the LC and control labels shuffled across points, 2,000 times, two sided)
   bar.csv            (a) R2(agreement + reliability) - R2(agreement) >= 0.15; (b) gap p > 0.05 in the agreement +
                      reliability model
4  cohorts.csv        per LC held-out cohort and block: samples per class, split half reliability, within cohort AUROC
                      (within_study_cv.csv, mean over its rows), held-out AUROC (mean over the reduced grid rows of
                      loso_results.csv), agreement (mean effect_spearman of the points), guided gain (arm 2 minus arm 1)
                      and local AUROC (arm 5) from results/concordance_guided/per_cohort.csv
   cohort_spearman.csv  Spearman of reliability with each of the others across those cohorts

  python scripts/transfer_model.py --out results/transfer_model
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from concordance_noise_ceiling import g_of, half_split, sb, spearman  # noqa: E402

from chorus.cli import load_config  # noqa: E402
from chorus.schema import binary_label  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = [0.2, 0.3, 0.4, 0.5]
POINT_SOURCES = {  # block of points.csv -> (config, contrast)
    "transcriptomics": ("configs/lc_transcriptomics_blood_v4.yaml", "lc"),
    "proteomics_affinity": ("configs/lc_proteomics_affinity_v5.yaml", "lc"),
    "proteomics_ms": ("configs/lc_proteomics_ms_v4.yaml", "lc"),
    "metabolomics": ("configs/lc_metabolomics.yaml", "lc"),
    "sex_transcriptomics": ("configs/lc_transcriptomics_blood.yaml", "sex"),
    "acute_proteomics": ("configs/control_acute_proteomics.yaml", "acute"),
}
LOSO = {"transcriptomics": "results/lc_transcriptomics_blood_v4", "proteomics_affinity": "results/lc_proteomics_affinity_v5",
        "proteomics_ms": "results/lc_proteomics_ms_v4", "metabolomics": "results/lc_metabolomics"}
GUIDED = {"transcriptomics": "transcriptomics_v4", "proteomics_affinity": "proteomics_affinity_v5", "proteomics_ms": "proteomics_ms_v4", "metabolomics": "metabolomics"}
RED_H = ["identity", "center_scale", "rank", "fsqn", "combat_addon", "chorus", "chorus_protected", "nested"]


# ------------------------------------------------------------------------------------------------ 1
def stratified(pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for t in THRESHOLDS:
        both = (pairs["reliability_a_shared"] > t) & (pairs["reliability_b_shared"] > t)
        for layer, d in pairs.groupby("layer", sort=False):
            for name, m in (("both_above", both.loc[d.index]), ("rest", ~both.loc[d.index])):
                v = d.loc[m, "below_2p5"].dropna().astype(bool)
                rows.append({"threshold": t, "layer": layer, "stratum": name, "n_pairs": len(v), "n_below_2p5": int(v.sum()), "fraction_below": v.mean() if len(v) else np.nan,
                             "more_than_half": bool(v.sum() > len(v) / 2) if len(v) else None})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ 2
def labelled(config: str, contrast: str) -> tuple:
    cs, lc, cfg = load_config(str(ROOT / config))
    X, M = cs.pooled(cfg["filtered_features"])
    if contrast == "sex":
        in_c = M["lc_status"].isin([lc.positive, lc.negative])
        y = pd.Series(np.nan, index=M.index)
        y[in_c & (M["inferred_sex"] == "female")] = 1.0
        y[in_c & (M["inferred_sex"] == "male")] = 0.0
    else:
        y = binary_label(M, lc.positive, lc.negative)
    keep = y.notna()
    return X[keep], y[keep].astype(int), M.loc[keep, "study_id"].astype(str)


def reliability(n_splits: int, seed: int) -> pd.DataFrame:
    rows = []
    for k, (block, (config, contrast)) in enumerate(POINT_SOURCES.items()):
        X, y, s = labelled(config, contrast)
        for j, c in enumerate(pd.unique(s)):
            m = (s == c).to_numpy()
            Xc, yc = X[m].to_numpy(float), y[m].to_numpy()
            rng = np.random.default_rng([seed, k, j])
            r = np.array([spearman(g_of(Xc, yc, a), g_of(Xc, yc, b)) for a, b in (half_split(yc, rng) for _ in range(n_splits))])
            g = g_of(Xc, yc, np.arange(len(yc)))
            rows.append({"block": block, "contrast": contrast, "cohort": c, "n_pos": int((yc == 1).sum()), "n_neg": int((yc == 0).sum()), "n_features": Xc.shape[1],
                         "split_half_r": np.nanmean(r) if np.isfinite(r).any() else np.nan, "reliability": sb(np.nanmean(r)) if np.isfinite(r).any() else np.nan,
                         "mean_abs_g": float(np.nanmean(np.abs(g)))})
            print(f"{block} {c}: reliability {rows[-1]['reliability']:.3f}", flush=True)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ 3
def ols(y: np.ndarray, X: np.ndarray) -> tuple:
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ beta
    r2 = 1 - (res**2).sum() / ((y - y.mean()) ** 2).sum()
    return beta, r2, res


def design(P: pd.DataFrame, terms: list) -> np.ndarray:
    cols = []
    for t in terms:
        if "*" in t:
            a, b = t.split("*")
            cols.append(P[a].to_numpy() * P[b].to_numpy())
        else:
            cols.append(P[t].to_numpy())
    return np.column_stack(cols)


MODELS = {
    "agreement": ["effect_spearman"],
    "reliability": ["reliability"],
    "agreement+reliability": ["effect_spearman", "reliability"],
    "agreement*reliability": ["effect_spearman", "reliability", "effect_spearman*reliability"],
    "mean_abs_g": ["mean_abs_g"],
    "agreement+mean_abs_g": ["effect_spearman", "mean_abs_g"],
    "agreement*mean_abs_g": ["effect_spearman", "mean_abs_g", "effect_spearman*mean_abs_g"],
}


def models(P: pd.DataFrame, n_boot: int, n_perm: int, seed: int) -> tuple:
    rng = np.random.default_rng(seed)
    y = P["auroc"].to_numpy()
    lc = (P["group"] == "LC").to_numpy().astype(float)
    clusters = P["contrast"] + "|" + P["block"] + "|" + P["held_out"]
    by_contrast = {c: [np.where(clusters.to_numpy() == k)[0] for k in pd.unique(clusters[P["contrast"] == c])] for c in pd.unique(P["contrast"])}
    rows, gaps = [], []
    for name, terms in MODELS.items():
        X = design(P, terms)
        beta, r2, res = ols(y, X)
        boots = []
        for _ in range(n_boot):
            idx = np.concatenate([np.concatenate([g[i] for i in rng.integers(0, len(g), len(g))]) for g in by_contrast.values()])
            b, r, _ = ols(y[idx], X[idx])
            boots.append(np.r_[b, r])
        B = np.array(boots)
        ci = np.nanpercentile(B, [2.5, 97.5], axis=0)
        row = {"model": name, "n_points": len(y), "r2": r2, "r2_lo": ci[0, -1], "r2_hi": ci[1, -1]}
        for i, t in enumerate(["intercept", *terms]):
            row[f"coef_{t}"], row[f"coef_{t}_lo"], row[f"coef_{t}_hi"] = beta[i], ci[0, i], ci[1, i]
        rows.append(row)
        # LC minus control gap: indicator added to the model, permutation of the labels
        bg, _, _ = ols(y, np.column_stack([X, lc]))
        null = np.array([ols(y, np.column_stack([X, rng.permutation(lc)]))[0][-1] for _ in range(n_perm)])
        gaps.append({"model": name, "gap_LC_minus_control": bg[-1], "mean_residual_LC": res[lc == 1].mean(), "mean_residual_control": res[lc == 0].mean(),
                     "p_permutation": float((1 + (np.abs(null) >= abs(bg[-1])).sum()) / (1 + n_perm)), "n_perm": n_perm})
    return pd.DataFrame(rows), pd.DataFrame(gaps)


# ------------------------------------------------------------------------------------------------ 4
def cohort_table(P: pd.DataFrame, rel: pd.DataFrame) -> pd.DataFrame:
    gp = pd.read_csv(ROOT / "results/concordance_guided/per_cohort.csv").pivot_table(index=["block", "held_out"], columns="arm", values="auroc")
    rows = []
    for block, folder in LOSO.items():
        r = pd.read_csv(ROOT / folder / "loso_results.csv")
        r = r[r["auroc"].notna() & r["harmoniser"].isin(RED_H) & r["selector"].isin(["variance", "differential", "invariant"]) & r["classifier"].isin(["logreg", "rf"])]
        w = pd.read_csv(ROOT / folder / "within_study_cv.csv")
        for c in sorted(r["held_out"].unique()):
            rr = rel[(rel["block"] == block) & (rel["cohort"] == c)].iloc[0]
            g = gp.loc[(GUIDED[block], c)] if (GUIDED[block], c) in gp.index else pd.Series(dtype=float)
            rows.append({"block": block, "cohort": c, "n_LC": rr["n_pos"], "n_COVID_recovered": rr["n_neg"], "reliability": rr["reliability"], "mean_abs_g": rr["mean_abs_g"],
                         "within_cohort_auroc": w.loc[w["study"] == c, "auroc"].mean(), "held_out_auroc": r.loc[r["held_out"] == c, "auroc"].mean(),
                         "agreement": P.loc[(P["block"] == block) & (P["held_out"] == c), "effect_spearman"].mean(),
                         "guided_gain": g.get("guided", np.nan) - g.get("all", np.nan), "local_auroc": g.get("calib", np.nan)})
    return pd.DataFrame(rows)


def cohort_spearman(T: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col in ["n_LC", "n_COVID_recovered", "mean_abs_g", "within_cohort_auroc", "held_out_auroc", "agreement", "guided_gain", "local_auroc"]:
        d = T[["reliability", col]].dropna()
        rho, p = stats.spearmanr(d["reliability"], d[col]) if len(d) >= 3 else (np.nan, np.nan)
        rows.append({"with": col, "n_cohorts": len(d), "spearman": rho, "p": p})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-splits", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stratified(pd.read_csv(ROOT / "results/concordance_noise_ceiling/pairs.csv")).to_csv(out / "stratified.csv", index=False)
    rel = reliability(a.n_splits, a.seed)
    rel.to_csv(out / "reliability.csv", index=False)
    P = pd.read_csv(ROOT / "results/concordance_vs_transfer_full/points.csv")
    P = P[P["group"].isin(["LC", "control"])].copy()
    P = P.merge(rel[["block", "cohort", "reliability", "mean_abs_g"]], left_on=["block", "held_out"], right_on=["block", "cohort"], how="left").drop(columns="cohort")
    if P["reliability"].isna().any():
        raise SystemExit(f"points without reliability: {P.loc[P['reliability'].isna(), ['block', 'held_out']].drop_duplicates().values.tolist()}")
    P = P.dropna(subset=["auroc", "effect_spearman"])
    P.to_csv(out / "points.csv", index=False)
    M, G = models(P, a.n_boot, a.n_perm, a.seed)
    M.to_csv(out / "models.csv", index=False)
    G.to_csv(out / "gap.csv", index=False)
    r2 = M.set_index("model")["r2"]
    gp = G.set_index("model").loc["agreement+reliability", "p_permutation"]
    bar = pd.DataFrame([
        {"criterion": "(a) R2(agreement + reliability) - R2(agreement) >= 0.15", "value": r2["agreement+reliability"] - r2["agreement"], "met": bool(r2["agreement+reliability"] - r2["agreement"] >= 0.15)},
        {"criterion": "(b) LC minus control gap not significant in the agreement + reliability model (p > 0.05)", "value": gp, "met": bool(gp > 0.05)},
    ])
    bar = pd.concat([bar, pd.DataFrame([{"criterion": "two factor account holds: (a) and (b)", "value": np.nan, "met": bool(bar["met"].all())}])], ignore_index=True)
    bar.to_csv(out / "bar.csv", index=False)
    T = cohort_table(P, rel)
    T.to_csv(out / "cohorts.csv", index=False)
    cohort_spearman(T).to_csv(out / "cohort_spearman.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"point_sources": POINT_SOURCES, "n_splits": a.n_splits, "n_boot": a.n_boot, "n_perm": a.n_perm, "seed": a.seed, "thresholds": THRESHOLDS}, indent=1))
    print(M[["model", "r2", "r2_lo", "r2_hi"]].round(3).to_string(index=False)); print(G.round(4).to_string(index=False)); print(bar.to_string(index=False))


if __name__ == "__main__":
    main()
