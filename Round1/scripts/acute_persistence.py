"""Persistence of the acute signature: does an acute COVID-19 versus never infected model score LC higher than
COVID_recovered after the acute phase?

Proteins (--out/protein):
  training  the cohorts and samples of configs/control_acute_proteomics.yaml (acute versus never_infected, proteins
            shared by all three cohorts). Per label free harmoniser (identity, center_scale, quantile_ref, fsqn, rank):
            the harmoniser is fitted on the training cohorts, the invariant selector (k = 50, the training cohorts as
            environments) chooses the proteins, then the harmoniser is refitted on the chosen proteins and logreg is
            fitted on them. The model is fitted once and never refitted.
  targets   the post acute LC and COVID_recovered samples of configs/concordance_proteomics_all_v4.yaml (one per
            subject, the earliest; each cohort imputed as there). Subjects in the training samples are removed (the
            same subject_id). A cohort shares the proteins when it measures at least --min-shared of the training
            proteins; the training proteins are first restricted to those measured in every sharing cohort that keeps
            at least --min-per-class subjects per group, so one model scores all of them. A cohort is scored only when
            it measures every protein of the model. The harmoniser is applied to each cohort on its own (addon).
RNA programmes (--out/rna_programmes):
  training  the acute contrast of scripts/mosaic_programmes.py (configs/mosaic_acute_rna.yaml and
            configs/mosaic_acute_protein.yaml): programme scores of every acute environment, harmonisers identity
            and center_scale, logreg on the programmes scored in every training environment and every target cohort.
  targets   the RNA cohorts of configs/mosaic_lc_rna_v4.yaml (the RNA side of the lc_v4 contrast).

Acute score: the logreg log odds of acute. Per target cohort: n per group, Hedges g of the score (LC minus
COVID_recovered) with its variance and a stratified bootstrap interval, AUROC with a stratified bootstrap interval;
where days_since_infection is recorded with more than one value in the cohort, the OLS coefficient of LC in
score ~ LC + days on the samples with days, divided by the pooled within group SD of the score on those samples
(g_adjusted, its variance from the coefficient's standard error). Random effects meta-analysis (DerSimonian and
Laird) across the scored cohorts of g, and of g_adjusted where available with g elsewhere (g_days).

  per_cohort.csv   harmoniser, cohort and the statistics above
  meta.csv         harmoniser, estimate (g or g_days), pooled effect, se, interval, z, p, Q, I2, tau2, k
  model.csv        harmoniser and the model's features with their coefficients
  sharing.csv      per target cohort: samples and subjects per group before and after removing training subjects,
                   training features measured, shares, scored
  scores.csv       per harmoniser and target sample: cohort, group, score

Acute check (--acute-check COHORT, version 5: Woodruff2023): the protein model scores the cohort's acute samples
(sample_phase acute) together with its post acute LC and COVID_recovered samples, harmonised as one cohort (addon);
acute_check.csv gives per harmoniser the mean score of each group, Hedges g (acute minus post acute) and the AUROC of the
score for acute versus post acute. The acute samples are used only here.

  python scripts/acute_persistence.py --out results/acute_persistence
  python scripts/acute_persistence.py --out results/acute_persistence_v5 --target-protein configs/concordance_proteomics_all_v5.yaml --acute-check Woodruff2023
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mosaic_programmes import gene_set_table, load_layer  # noqa: E402

from chorus.cli import load_config  # noqa: E402
from chorus.eval.metrics import bootstrap_ci  # noqa: E402
from chorus.features import make_selector  # noqa: E402
from chorus.features.invariant import _hedges_g, meta_analyse  # noqa: E402
from chorus.harmonise import make_harmoniser  # noqa: E402
from chorus.integrate.programmes import common_programmes  # noqa: E402
from chorus.models import make_classifier  # noqa: E402
from chorus.schema import binary_label  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TRAIN_PROTEIN = "configs/control_acute_proteomics.yaml"
TARGET_PROTEIN = "configs/concordance_proteomics_all_v4.yaml"
TARGET_RNA = "configs/mosaic_lc_rna_v4.yaml"
LABEL_FREE = ["identity", "center_scale", "quantile_ref", "fsqn", "rank"]
PROGRAMME_HARMONISERS = ["identity", "center_scale"]  # those of the mosaic acute model


def measured(X: pd.DataFrame) -> pd.Index:
    sd = X.std(axis=0, ddof=1)
    return X.columns[X.notna().all(axis=0).values & (sd > 0).values]


def log_odds(clf, X: pd.DataFrame) -> np.ndarray:
    p = np.clip(clf.predict_proba(X)[:, 1], 1e-12, 1 - 1e-12)
    return np.log(p / (1 - p))


def hedges_ci(a: np.ndarray, b: np.ndarray, n_boot: int, seed: int) -> tuple:
    rng = np.random.default_rng(seed)
    vals = [_hedges_g(rng.choice(a, len(a))[:, None], rng.choice(b, len(b))[:, None])[0][0] for _ in range(n_boot)]
    return tuple(np.nanpercentile(vals, [2.5, 97.5]))


def days_adjusted(score: np.ndarray, y: np.ndarray, days: np.ndarray) -> dict:
    """Coefficient of y in score ~ y + days (OLS) over the pooled within group SD of the score."""
    import statsmodels.api as sm

    ok = np.isfinite(days)
    if ok.sum() < 6 or len(np.unique(days[ok])) < 2 or min((y[ok] == 1).sum(), (y[ok] == 0).sum()) < 3:
        return {"n_days": int(ok.sum()), "g_adjusted": np.nan, "var_adjusted": np.nan, "days_coef": np.nan}
    s, yy, d = score[ok], y[ok], days[ok]
    fit = sm.OLS(s, sm.add_constant(np.column_stack([yy, d]))).fit()
    a, b = s[yy == 1], s[yy == 0]
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return {"n_days": int(ok.sum()), "g_adjusted": fit.params[1] / sp, "var_adjusted": (fit.bse[1] / sp) ** 2, "days_coef": fit.params[2]}


def cohort_stats(score: pd.Series, y: pd.Series, days: pd.Series, n_boot: int, seed: int) -> dict:
    s, yy = score.to_numpy(), y.to_numpy()
    a, b = s[yy == 1], s[yy == 0]
    g, v = (x[0] for x in _hedges_g(a[:, None], b[:, None]))
    glo, ghi = hedges_ci(a, b, n_boot, seed)
    auc, lo, hi = bootstrap_ci(yy, s, roc_auc_score, n_boot=n_boot, seed=seed)
    return {
        "n_lc": len(a), "n_recovered": len(b), "mean_lc": a.mean(), "mean_recovered": b.mean(),
        "g": g, "var_g": v, "g_lo": glo, "g_hi": ghi, "auroc": auc, "auroc_lo": lo, "auroc_hi": hi,
        **days_adjusted(s, yy, pd.to_numeric(days, errors="coerce").to_numpy(dtype=float)),
    }


def meta_table(per: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for h, d in per[per["scored"]].groupby("harmoniser", sort=False):
        for name in ("g", "g_days"):
            if name == "g":
                gg, vv = d["g"], d["var_g"]
            else:
                use = d["g_adjusted"].notna()
                gg, vv = d["g_adjusted"].where(use, d["g"]), d["var_adjusted"].where(use, d["var_g"])
            eff = {c: pd.DataFrame({"g": [g], "var": [v]}, index=["acute_score"]) for c, g, v in zip(d["cohort"], gg, vv)}
            m = meta_analyse(eff).iloc[0]
            rows.append(
                {"harmoniser": h, "estimate": name, "n_adjusted": int(d["g_adjusted"].notna().sum()) if name == "g_days" else 0, "effect": m["effect"], "se": m["se"],
                 "ci_lo": m["effect"] - 1.959964 * m["se"], "ci_hi": m["effect"] + 1.959964 * m["se"], "z": m["z"], "p": m["p"], "Q": m["Q"], "Q_p": m["Q_p"], "I2": m["I2"], "tau2": m["tau2"], "k": int(m["k"]),
                 "cohorts": ";".join(d["cohort"])}
            )
    return pd.DataFrame(rows)


def score_targets(models: dict, X: pd.DataFrame, M: pd.DataFrame, y: pd.Series, cohorts: list, n_boot: int, seed: int) -> tuple:
    """models: harmoniser -> (fitted harmoniser, classifier, features). Returns (per cohort rows, scores)."""
    rows, scores = [], []
    for h, (H, clf, feats) in models.items():
        for c in cohorts:
            m = (M["study_id"] == c).to_numpy()
            Xc = X.loc[m, feats]
            Z = H.transform(Xc, M.loc[m, "study_id"])
            s = pd.Series(log_odds(clf, Z), index=Xc.index)
            rows.append({"harmoniser": h, "cohort": c, "scored": True, **cohort_stats(s, y[m], M.loc[m].get("days_since_infection", pd.Series(np.nan, index=Xc.index)), n_boot, seed)})
            scores.append(pd.DataFrame({"harmoniser": h, "cohort": c, "sample": Xc.index, "group": np.where(y[m] == 1, "LC", "COVID_recovered"), "score": s.values}))
    return rows, scores


def sharing_row(c: str, M: pd.DataFrame, y: pd.Series, keep: pd.Series) -> dict:
    m = M["study_id"] == c
    sub = lambda mask, v: int(M.loc[mask & (y == v), "subject_id"].nunique())  # noqa: E731
    return {
        "cohort": c, "n_lc_before": sub(m, 1), "n_recovered_before": sub(m, 0),
        "n_lc_removed_training_subject": sub(m & ~keep, 1), "n_recovered_removed_training_subject": sub(m & ~keep, 0),
        "n_lc": sub(m & keep, 1), "n_recovered": sub(m & keep, 0),
    }


def acute_check(models: dict, config: str, cohort: str, n_boot: int, seed: int) -> pd.DataFrame:
    """Score ``cohort``'s acute and post acute (LC, COVID_recovered) samples with each protein model, harmonised together."""
    import yaml

    from chorus.io import Cohort

    entry = next(c for c in yaml.safe_load((ROOT / config).read_text())["cohorts"] if c["study_id"] == cohort)
    c = Cohort.from_files(cohort, "proteomics", str(ROOT / entry["matrix"]), str(ROOT / entry["meta"]), feature_space="uniprot")
    phase, st = c.meta["sample_phase"].astype(str), c.meta["lc_status"].astype(str)
    acute = phase == "acute"
    post = (phase == "post_acute") & st.isin(["LC", "COVID_recovered"])
    keep = (acute | post).to_numpy()
    X = c.X.loc[keep]
    X = X.fillna(X.median(axis=0))
    y = pd.Series(acute[keep].astype(int).to_numpy(), index=X.index)
    rows = []
    for h, (H, clf, feats) in models.items():
        Z = H.transform(X[feats], pd.Series(cohort, index=X.index))
        sc = log_odds(clf, Z)
        a, b = sc[y == 1], sc[y == 0]
        g, v = (x[0] for x in _hedges_g(a[:, None], b[:, None]))
        auc, lo, hi = bootstrap_ci(y.to_numpy(), sc, roc_auc_score, n_boot=n_boot, seed=seed)
        rows.append({"harmoniser": h, "cohort": cohort, "n_acute": len(a), "n_post_acute": len(b), "mean_acute": a.mean(), "mean_post_acute": b.mean(),
                     "mean_post_acute_lc": sc[(y == 0).to_numpy() & (st[keep] == "LC").to_numpy()].mean(),
                     "mean_post_acute_recovered": sc[(y == 0).to_numpy() & (st[keep] == "COVID_recovered").to_numpy()].mean(),
                     "g_acute_minus_post": g, "var_g": v, "auroc_acute_vs_post": auc, "auroc_lo": lo, "auroc_hi": hi})
    return pd.DataFrame(rows)


def protein(out: Path, min_shared: float, min_per_class: int, n_boot: int, seed: int, target: str = TARGET_PROTEIN, check: str = "") -> dict:
    cs, lc, cfg = load_config(str(ROOT / TRAIN_PROTEIN))
    Xtr, Mtr = cs.pooled(cfg["filtered_features"])
    ytr = binary_label(Mtr, lc.positive, lc.negative)
    Xtr, Mtr, ytr = Xtr[ytr.notna()], Mtr[ytr.notna()], ytr[ytr.notna()].astype(int)
    str_ = Mtr["study_id"].astype(str)
    train_subjects = set(Mtr["subject_id"].astype(str))
    train_feats = pd.Index(cfg["filtered_features"])

    ct, lt, ccfg = load_config(str(ROOT / target))
    X, M = ct.pooled(ccfg["filtered_features"])
    y = binary_label(M, lt.positive, lt.negative)
    X, M, y = X[y.notna()], M[y.notna()], y[y.notna()].astype(int)
    keep = ~M["subject_id"].astype(str).isin(train_subjects)
    share = []
    for c in ct.study_ids:
        r = sharing_row(c, M, y, keep)
        meas = measured(X.loc[(M["study_id"] == c).to_numpy()].reindex(columns=train_feats))
        r.update({"n_training_features": len(train_feats), "n_training_features_measured": len(meas), "shares": len(meas) >= min_shared * len(train_feats)})
        r["eligible"] = r["shares"] and min(r["n_lc"], r["n_recovered"]) >= min_per_class
        share.append(r)
    share = pd.DataFrame(share)
    eligible = list(share.loc[share["eligible"], "cohort"])
    X, M, y = X[keep.to_numpy()], M[keep], y[keep]
    feats = train_feats
    for c in eligible:
        feats = feats.intersection(measured(X.loc[(M["study_id"] == c).to_numpy()].reindex(columns=train_feats)))
    models, coef = {}, []
    for h in LABEL_FREE:
        H = make_harmoniser(h).fit(Xtr[feats], str_)
        sel = make_selector(name="invariant", k=50).fit(H.transform(Xtr[feats], str_), ytr, str_)
        chosen = list(sel.selected_)
        H2 = make_harmoniser(h).fit(Xtr[chosen], str_)
        clf = make_classifier("logreg", seed).fit(H2.transform(Xtr[chosen], str_), ytr)
        models[h] = (H2, clf, chosen)
        lr = clf.steps[-1][1] if hasattr(clf, "steps") else clf
        coef.append(pd.DataFrame({"harmoniser": h, "feature": chosen, "coefficient": np.ravel(lr.coef_)}))
    scored = [c for c in eligible if all(set(m[2]) <= set(measured(X.loc[(M["study_id"] == c).to_numpy(), m[2]])) for m in models.values())]
    share["n_model_features_pool"] = len(feats)
    share["scored"] = share["cohort"].isin(scored)
    rows, scores = score_targets(models, X, M, y, scored, n_boot, seed)
    if check:
        out.mkdir(parents=True, exist_ok=True)
        A = acute_check(models, target, check, n_boot, seed)
        A.to_csv(out / "acute_check.csv", index=False)
        print(A.round(3).to_string(index=False))
    return write(out, rows, scores, share, pd.concat(coef, ignore_index=True), {"training": TRAIN_PROTEIN, "targets": target, "acute_check": check or None, "harmonisers": LABEL_FREE, "selector": {"name": "invariant", "k": 50}, "n_training_features": len(train_feats), "n_features_pool": len(feats), "eligible": eligible, "scored": scored, "min_shared": min_shared})


def rna_programmes(out: Path, min_members: int, min_per_class: int, n_boot: int, seed: int) -> dict:
    T = gene_set_table()
    train = {layer: load_layer(cfg, layer, T, min_members) for layer, cfg in {"rna": "configs/mosaic_acute_rna.yaml", "protein": "configs/mosaic_acute_protein.yaml"}.items()}
    S_tr = pd.concat([train[l]["S"] for l in train])
    y_tr = pd.concat([train[l]["y"] for l in train]).astype(int)
    e_tr = pd.concat([train[l]["env"] for l in train])
    train_subjects = set(pd.concat([train[l]["M"]["subject_id"] for l in train]).astype(str))
    tgt = load_layer(TARGET_RNA, "rna", T, min_members)
    S, M, y, env = tgt["S"], tgt["M"], tgt["y"].astype(int), tgt["env"]
    keep = ~M["subject_id"].astype(str).isin(train_subjects)
    M2 = M.assign(study_id=env)
    share = pd.DataFrame([sharing_row(c, M2, y, keep) for c in pd.unique(env)])
    share["eligible"] = np.minimum(share["n_lc"], share["n_recovered"]) >= min_per_class
    targets = list(share.loc[share["eligible"], "cohort"])
    S, M2, y = S[keep], M2[keep], y[keep]
    progs = [p for p in common_programmes(S_tr, e_tr, list(pd.unique(e_tr))) if p in set(common_programmes(S, M2["study_id"], targets))]
    models, coef = {}, []
    for h in PROGRAMME_HARMONISERS:
        H = make_harmoniser(h).fit(S_tr[progs], e_tr)
        clf = make_classifier("logreg", seed).fit(H.transform(S_tr[progs], e_tr), y_tr)
        models[h] = (H, clf, progs)
        lr = clf.steps[-1][1] if hasattr(clf, "steps") else clf
        coef.append(pd.DataFrame({"harmoniser": h, "feature": progs, "coefficient": np.ravel(lr.coef_)}))
    share["n_programmes"] = len(progs)
    share["scored"] = share["eligible"]
    rows, scores = score_targets(models, S, M2, y, targets, n_boot, seed)
    return write(out, rows, scores, share, pd.concat(coef, ignore_index=True), {"training_environments": list(pd.unique(e_tr)), "targets": TARGET_RNA, "harmonisers": PROGRAMME_HARMONISERS, "n_programmes": len(progs), "min_members": min_members, "scored": targets})


def write(out: Path, rows, scores, share, coef, info) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    per = pd.DataFrame(rows)
    per.to_csv(out / "per_cohort.csv", index=False)
    meta = meta_table(per)
    meta.to_csv(out / "meta.csv", index=False)
    share.to_csv(out / "sharing.csv", index=False)
    coef.to_csv(out / "model.csv", index=False)
    pd.concat(scores, ignore_index=True).to_csv(out / "scores.csv", index=False)
    (out / "run_info.json").write_text(json.dumps(info, indent=1))
    print(share.to_string(index=False))
    print(per.round(3).to_string(index=False))
    print(meta.round(3).to_string(index=False), flush=True)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-shared", type=float, default=0.5, help="fraction of the training proteins a target cohort must measure to share them")
    ap.add_argument("--min-per-class", type=int, default=5)
    ap.add_argument("--min-members", type=int, default=5)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--parts", default="protein,rna_programmes")
    ap.add_argument("--target-protein", default=TARGET_PROTEIN, help="config of the post acute protein cohorts to score")
    ap.add_argument("--acute-check", default="", help="protein cohort whose acute samples are scored as a check")
    a = ap.parse_args()
    out = Path(a.out)
    parts = a.parts.split(",")
    if "protein" in parts:
        protein(out / "protein", a.min_shared, a.min_per_class, a.n_boot, a.seed, a.target_protein, a.acute_check)
    if "rna_programmes" in parts:
        rna_programmes(out / "rna_programmes", a.min_members, a.min_per_class, a.n_boot, a.seed)


if __name__ == "__main__":
    main()
