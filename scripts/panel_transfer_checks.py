"""Checks of the CORSIP2026 panel transfer (results/panel_transfer_CORSIP2026); existing results are read, never rewritten.

The configurations checked are the five with the highest CORSIP2026 AUROC in results/panel_transfer_CORSIP2026/loso_results.csv
and the median configuration (the lower middle one when the configurations that ran are sorted by AUROC). A configuration
is refitted exactly as chorus.eval.run_loso fits it (harmoniser on the training cohorts, selector, sample_matrix
covariate, classifier; nested chooses its harmoniser by the inner leave one study out); ``refit`` checks that the refit
reproduces the stored AUROC.

  confounding   sex and age of CORSIP2026 by LC status; AUROC of sex alone and of age alone; transfer AUROC within each sex
                where both classes have at least 3 samples (every configuration, from the stored predictions); logistic
                regression of LC on the standardised model score, sex and standardised age (odds ratio per SD of score,
                Wald 95% interval); the largest model weights of the checked configurations (logreg: coefficient on
                standardised features; rf: impurity importance) with chromosome and the sex difference of each gene in
                CORSIP2026 (Welch t test of log CPM, female minus male, Benjamini-Hochberg)
  permute       lc_status permuted within each training cohort (class counts kept) --n-perm times, the checked
                configurations refitted and CORSIP2026 scored; null mean, SD and quantiles, p = (1 + null >= observed) /
                (1 + n_perm)
  concordance   Hedges g of LC versus COVID_recovered per cohort on the 495 genes; Spearman correlation and sign agreement
                of CORSIP2026 with each training cohort; leave one training cohort out: the checked configurations
                refitted without that cohort and CORSIP2026 scored
  random        the version 4 leave one study out (configs/lc_transcriptomics_blood_v4_panel495.yaml grid) on --n-sets
                random sets of 495 genes drawn from the 10,053 genes expressed in every v4 cohort

  python scripts/panel_transfer_checks.py confounding --out results/panel_transfer_CORSIP2026_checks
  python scripts/panel_transfer_checks.py permute --out results/panel_transfer_CORSIP2026_checks --n-perm 200 [--start 0 --stop 50]
  python scripts/panel_transfer_checks.py permutation_summary --out results/panel_transfer_CORSIP2026_checks
  python scripts/panel_transfer_checks.py concordance --out results/panel_transfer_CORSIP2026_checks
  python scripts/panel_transfer_checks.py random --out results/lc_transcriptomics_blood_v4_panel495/random_genes --n-sets 20 [--start 0 --stop 5]
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from chorus.cli import load_config
from chorus.eval import run_loso
from chorus.eval.loso import _covariate_matrix, _fit_harmoniser, _forced_extra, _inner_select_harmoniser, selector_label
from chorus.features import make_selector
from chorus.features.invariant import per_study_effects
from chorus.models import make_classifier
from chorus.schema import binary_label
from chorus.sex import gene_chromosomes

ROOT = Path(__file__).resolve().parents[1]
CONFIG = "configs/panel_transfer_CORSIP2026.yaml"
RESULTS = ROOT / "results" / "panel_transfer_CORSIP2026"
TEST = "CORSIP2026"
GENE_TABLE = "data/quant/GSE224615/tx2gene_gencode_v47.tsv"
HGNC = "data/ref/hgnc_complete_set.txt"
KEYS = ["harmoniser", "selector", "classifier"]


def load():
    cs, lc, _ = load_config(str(ROOT / CONFIG))
    X, M = cs.pooled()
    y = binary_label(M, lc.positive, lc.negative)
    keep = y.notna()
    X, M, y = X[keep], M[keep], y[keep].astype(int)
    s = M["study_id"].astype(str)
    return {"X": X, "M": M, "y": y, "s": s, "lc": lc, "cov": _covariate_matrix(M, lc.covariates), "prot": _covariate_matrix(M, lc.protect) if lc.protect else None}


def checked_configs(results: pd.DataFrame) -> pd.DataFrame:
    ok = results[results["auroc"].notna()].sort_values("auroc", ascending=False, kind="stable").reset_index(drop=True)
    top = ok.head(5).assign(role=[f"top{i + 1}" for i in range(5)])
    asc = ok.sort_values("auroc", kind="stable").reset_index(drop=True)
    med = asc.iloc[[(len(asc) - 1) // 2]].assign(role="median")
    return pd.concat([top, med], ignore_index=True)[["role", *KEYS, "auroc", "harmoniser_chosen"]]


def selector_cfg(lc, label: str) -> dict:
    return next(s for s in lc.selectors if selector_label(s) == label)


def fit_score(D: dict, h: str, sel_label: str, clf_name: str, y_train: pd.Series = None, train_mask: pd.Series = None, return_model: bool = False):
    """Fit one configuration on the training cohorts (``train_mask``, default every cohort but CORSIP2026) with labels
    ``y_train`` (default the curated labels) and return the CORSIP2026 scores (and the fitted parts)."""
    X, s, lc = D["X"], D["s"], D["lc"]
    te = (s == TEST).to_numpy()
    tr = ~te if train_mask is None else (train_mask.to_numpy() & ~te)
    y = D["y"] if y_train is None else y_train
    X_tr, X_te, y_tr, s_tr, s_te = X[tr], X[te], y[tr], s[tr], s[te]
    prot = None if D["prot"] is None else D["prot"].loc[X_tr.index]
    sc = selector_cfg(lc, sel_label)
    chosen = h
    if h == "nested":
        chosen, _ = _inner_select_harmoniser(lc.nested_candidates or lc.harmonisers, X_tr, s_tr, y_tr, sc, lc, prot)
    Xh_tr, Xh_te, _, _ = _fit_harmoniser(chosen, X_tr, s_tr, y_tr, X_te, s_te, prot, lc.harmoniser_params.get(chosen))
    sel = make_selector(**sc).fit(Xh_tr, y_tr, s_tr)
    forced = _forced_extra(lc, sel.selected_, Xh_tr.columns)
    Xs_tr = pd.concat([sel.transform(Xh_tr), Xh_tr[forced], D["cov"].loc[X_tr.index]], axis=1)
    Xs_te = pd.concat([sel.transform(Xh_te), Xh_te.reindex(columns=forced), D["cov"].loc[X_te.index]], axis=1)
    clf = make_classifier(clf_name, lc.random_state).fit(Xs_tr, y_tr)
    p = pd.Series(clf.predict_proba(Xs_te)[:, 1], index=X_te.index)
    if return_model:
        return p, {"clf": clf, "columns": list(Xs_tr.columns), "chosen": chosen}
    return p


def permuted_labels(y: pd.Series, s: pd.Series, rng) -> pd.Series:
    out = y.copy()
    for c in pd.unique(s):
        if c == TEST:
            continue
        idx = s.index[s == c]
        out.loc[idx] = rng.permutation(y.loc[idx].to_numpy())
    return out


def auroc(D, p: pd.Series) -> float:
    return float(roc_auc_score(D["y"].loc[p.index], p))


# ------------------------------------------------------------------------------------------------ confounding
def weights(model: dict, clf_name: str) -> pd.Series:
    est = model["clf"].steps[-1][1] if hasattr(model["clf"], "steps") else model["clf"]
    w = np.ravel(est.coef_) if clf_name.startswith("logreg") else est.feature_importances_
    return pd.Series(w, index=model["columns"])


def sex_difference(D: dict) -> pd.DataFrame:
    m = (D["s"] == TEST).to_numpy()
    X, sex = D["X"][m], D["M"].loc[m, "sex"].astype(str)
    f, ml = X[(sex == "female").to_numpy()], X[(sex == "male").to_numpy()]
    t, p = stats.ttest_ind(f, ml, equal_var=False)
    q = pd.Series(p, index=X.columns).rank(method="first")
    fdr = pd.Series(p, index=X.columns) * len(p) / q
    fdr = fdr.sort_values(ascending=False).cummin().reindex(X.columns).clip(upper=1)
    return pd.DataFrame({"female_minus_male": f.mean() - ml.mean(), "sex_t": t, "sex_p": p, "sex_fdr": fdr})


def confounding(out: Path, top_n: int) -> None:
    D = load()
    res = pd.read_csv(RESULTS / "loso_results.csv")
    pred = pd.read_csv(RESULTS / "loso_predictions.csv")
    C = checked_configs(res)
    m = (D["s"] == TEST).to_numpy()
    meta, y = D["M"][m], D["y"][m]
    sex, age = meta["sex"].astype(str), pd.to_numeric(meta["age"], errors="coerce")
    rows = []
    for g, lab in [(1, "LC"), (0, "COVID_recovered")]:
        k = (y == g).to_numpy()
        rows.append({"group": lab, "n": int(k.sum()), "female": int((sex[k] == "female").sum()), "male": int((sex[k] == "male").sum()),
                     "age_mean": age[k].mean(), "age_sd": age[k].std(), "age_min": age[k].min(), "age_max": age[k].max()})
    demo = pd.DataFrame(rows)
    ct = pd.crosstab(y, sex)
    demo.attrs = {}
    fisher_p = stats.fisher_exact(ct.values)[1] if ct.shape == (2, 2) else np.nan
    mw_p = stats.mannwhitneyu(age[(y == 1).to_numpy()], age[(y == 0).to_numpy()]).pvalue
    single = pd.DataFrame([
        {"predictor": "sex (female = 1)", "auroc": roc_auc_score(y, (sex == "female").astype(int)), "test": "Fisher exact", "p": fisher_p},
        {"predictor": "age", "auroc": roc_auc_score(y, age), "test": "Mann-Whitney", "p": mw_p},
    ])
    demo.to_csv(out / "corsip_demographics.csv", index=False)
    single.to_csv(out / "corsip_single_predictors.csv", index=False)
    print(demo.round(2).to_string(index=False)); print(single.round(4).to_string(index=False))

    # transfer AUROC within each sex, every configuration
    sx = pred.merge(meta[["sex"]].rename_axis("sample").reset_index(), on="sample")
    within = []
    for key, g in sx.groupby(KEYS):
        row = dict(zip(KEYS, key))
        row["auroc_all"] = roc_auc_score(g["y"], g["p"])
        for sv in ("female", "male"):
            h = g[g["sex"] == sv]
            n1, n0 = int((h["y"] == 1).sum()), int((h["y"] == 0).sum())
            row[f"n_lc_{sv}"], row[f"n_recovered_{sv}"] = n1, n0
            row[f"auroc_{sv}"] = roc_auc_score(h["y"], h["p"]) if min(n1, n0) >= 3 else np.nan
        within.append(row)
    W = pd.DataFrame(within).sort_values("auroc_all", ascending=False)
    W.to_csv(out / "within_sex_auroc.csv", index=False)

    # score adjusted for sex and age, checked configurations (and every configuration)
    import statsmodels.api as sm

    adj = []
    for key, g in sx.groupby(KEYS):
        g = g.set_index("sample")
        Z = pd.DataFrame({"score_z": (g["p"] - g["p"].mean()) / g["p"].std(ddof=1), "female": (meta.loc[g.index, "sex"] == "female").astype(float),
                          "age_z": (age.loc[g.index] - age.mean()) / age.std(ddof=1)})
        try:
            fit = sm.Logit(g["y"].astype(float), sm.add_constant(Z)).fit(disp=0, maxiter=200)
            b, se = fit.params["score_z"], fit.bse["score_z"]
            r = {"or_per_sd": np.exp(b), "or_lo": np.exp(b - 1.959964 * se), "or_hi": np.exp(b + 1.959964 * se), "p": fit.pvalues["score_z"], "converged": bool(fit.mle_retvals.get("converged", True))}
        except Exception as e:  # perfect separation and the like
            r = {"or_per_sd": np.nan, "or_lo": np.nan, "or_hi": np.nan, "p": np.nan, "converged": False, "error": str(e)[:120]}
        adj.append({**dict(zip(KEYS, key)), "auroc": roc_auc_score(g["y"], g["p"]), **r})
    A = pd.DataFrame(adj).sort_values("auroc", ascending=False)
    A.to_csv(out / "adjusted_association.csv", index=False)

    # weights of the checked configurations, refitted on the curated labels
    chrom = gene_chromosomes(str(ROOT / GENE_TABLE))
    sym = pd.read_csv(ROOT / HGNC, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna().drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
    sd = sex_difference(D)
    wt, check = [], []
    for r in C.itertuples():
        p, model = fit_score(D, r.harmoniser, r.selector, r.classifier, return_model=True)
        a = auroc(D, p)
        check.append({"role": r.role, "harmoniser": r.harmoniser, "selector": r.selector, "classifier": r.classifier, "auroc_stored": r.auroc, "auroc_refit": a, "harmoniser_chosen": model["chosen"]})
        w = weights(model, r.classifier).drop([c for c in model["columns"] if c.startswith("cov_")], errors="ignore")
        top = w.reindex(w.abs().sort_values(ascending=False).index).head(top_n)
        for rank, (gid, v) in enumerate(top.items(), 1):
            wt.append({"role": r.role, "harmoniser": r.harmoniser, "selector": r.selector, "classifier": r.classifier, "rank": rank, "gene_id": gid, "symbol": sym.get(gid),
                       "chromosome": chrom.get(gid), "weight": v, "weight_type": "coefficient" if r.classifier.startswith("logreg") else "impurity_importance", **sd.loc[gid].to_dict()})
    pd.DataFrame(check).to_csv(out / "refit_check.csv", index=False)
    WT = pd.DataFrame(wt)
    WT.to_csv(out / "top_weights.csv", index=False)
    sd.assign(symbol=sd.index.map(sym), chromosome=sd.index.map(chrom)).rename_axis("gene_id").to_csv(out / "corsip_sex_difference.csv")
    print(pd.DataFrame(check).round(4).to_string(index=False))
    print(A.head(8).round(3).to_string(index=False))


# ------------------------------------------------------------------------------------------------ permutation
def permute(out: Path, n_perm: int, seed: int, start: int = 0, stop: int = None) -> None:
    """Permutations start to stop (default all) into permutation_null_<start>.csv; permutation i uses the seed
    (seed, i), so chunks run in parallel give the same draws as one run. ``summarise`` merges the chunks."""
    D = load()
    C = checked_configs(pd.read_csv(RESULTS / "loso_results.csv"))
    rows = []
    t0 = time.time()
    stop = n_perm if stop is None else min(stop, n_perm)
    part = out / f"permutation_null_{start:03d}.csv"
    for i in range(start, stop):
        yp = permuted_labels(D["y"], D["s"], np.random.default_rng([seed, i]))
        for r in C.itertuples():
            try:
                a = auroc(D, fit_score(D, r.harmoniser, r.selector, r.classifier, y_train=yp))
            except Exception as e:
                a = np.nan
            rows.append({"permutation": i, "role": r.role, "harmoniser": r.harmoniser, "selector": r.selector, "classifier": r.classifier, "auroc": a})
        if i % 10 == 9:
            print(f"permutation {i + 1}: {time.time() - t0:.0f} s", flush=True)
            pd.DataFrame(rows).to_csv(part, index=False)
    pd.DataFrame(rows).to_csv(part, index=False)


def permutation_summary(out: Path) -> None:
    C = checked_configs(pd.read_csv(RESULTS / "loso_results.csv"))
    parts = sorted(out.glob("permutation_null_*.csv"))
    P = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True).sort_values(["permutation", "role"], kind="stable")
    P.to_csv(out / "permutation_null.csv", index=False)
    summ = []
    for r in C.itertuples():
        a = P[(P["role"] == r.role)]["auroc"].dropna().to_numpy()
        summ.append({"role": r.role, "harmoniser": r.harmoniser, "selector": r.selector, "classifier": r.classifier, "auroc_observed": r.auroc,
                     "n_permutations": len(a), "null_mean": a.mean(), "null_sd": a.std(ddof=1), "null_q05": np.quantile(a, 0.05), "null_q50": np.quantile(a, 0.5),
                     "null_q95": np.quantile(a, 0.95), "null_max": a.max(), "p_permutation": (1 + (a >= r.auroc).sum()) / (1 + len(a))})
    S = pd.DataFrame(summ)
    S.to_csv(out / "permutation_summary.csv", index=False)
    print(S.round(3).to_string(index=False))


# ------------------------------------------------------------------------------------------------ concordance
def concordance(out: Path) -> None:
    D = load()
    C = checked_configs(pd.read_csv(RESULTS / "loso_results.csv"))
    eff = per_study_effects(D["X"], D["y"], D["s"])
    G = pd.DataFrame({c: e["g"] for c, e in eff.items()})
    G.rename_axis("gene_id").to_csv(out / "effects_495.csv")
    rows = []
    for c in G.columns:
        if c == TEST:
            continue
        ok = G[[c, TEST]].dropna()
        rho, p = stats.spearmanr(ok[c], ok[TEST])
        rows.append({"training_cohort": c, "n_genes": len(ok), "spearman_with_corsip": rho, "p": p, "sign_agreement": float((np.sign(ok[c]) == np.sign(ok[TEST])).mean())})
    conc = pd.DataFrame(rows)
    loo = []
    for c in conc["training_cohort"]:
        mask = D["s"] != c
        for r in C.itertuples():
            try:
                a = auroc(D, fit_score(D, r.harmoniser, r.selector, r.classifier, train_mask=mask))
                err = None
            except Exception as e:
                a, err = np.nan, str(e)[:120]
            loo.append({"left_out": c, "role": r.role, "harmoniser": r.harmoniser, "selector": r.selector, "classifier": r.classifier, "auroc_full": r.auroc, "auroc_without": a, "change": a - r.auroc, "error": err})
        print(f"left out {c}", flush=True)
    L = pd.DataFrame(loo)
    L.to_csv(out / "leave_one_training_cohort_out.csv", index=False)
    mean_change = L.groupby("left_out")["change"].mean().rename("mean_auroc_change")
    conc = conc.merge(mean_change, left_on="training_cohort", right_index=True)
    rho, p = stats.spearmanr(conc["spearman_with_corsip"], conc["mean_auroc_change"])
    conc.to_csv(out / "concordance_with_corsip.csv", index=False)
    (out / "concordance_summary.json").write_text(json.dumps({"spearman_concordance_vs_mean_auroc_change": rho, "p": p, "n_training_cohorts": len(conc)}, indent=1))
    print(conc.round(3).to_string(index=False)); print("Spearman concordance vs mean AUROC change:", round(rho, 3), round(p, 3))


# ------------------------------------------------------------------------------------------------ random genes
def random_genes(out: Path, n_sets: int, seed: int, start: int = 0, stop: int = None) -> None:
    """Gene set i is drawn with the seed (seed, i), so chunks run in parallel give the same sets as one run."""
    cs, lc, cfg = load_config(str(ROOT / "configs/lc_transcriptomics_blood_v4.yaml"))
    _, lc5, _ = load_config(str(ROOT / "configs/lc_transcriptomics_blood_v4_panel495.yaml"))
    universe = sorted(cfg["filtered_features"])
    out.mkdir(parents=True, exist_ok=True)
    for i in range(start, n_sets if stop is None else min(stop, n_sets)):
        genes = sorted(str(g) for g in np.random.default_rng([seed, i]).choice(universe, 495, replace=False))  # plain str: sklearn rejects numpy str_ names
        d = out / f"set{i:02d}"
        d.mkdir(exist_ok=True)
        (d / "genes.txt").write_text("\n".join(genes) + "\n")
        t0 = time.time()
        res, _, _ = run_loso(cs, lc5, features=pd.Index(genes), verbose=False)
        res.to_csv(d / "loso_results.csv", index=False)
        print(f"set {i}: {time.time() - t0:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("check", choices=["confounding", "permute", "permutation_summary", "concordance", "random"])
    ap.add_argument("--start", type=int, default=0, help="first permutation or gene set of this chunk")
    ap.add_argument("--stop", type=int, help="end (exclusive) of this chunk")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--n-sets", type=int, default=20)
    ap.add_argument("--top", type=int, default=20, help="largest weights reported per checked configuration")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.check == "confounding":
        confounding(out, a.top)
    elif a.check == "permute":
        permute(out, a.n_perm, a.seed, a.start, a.stop)
    elif a.check == "permutation_summary":
        permutation_summary(out)
    elif a.check == "concordance":
        concordance(out)
    else:
        random_genes(out, a.n_sets, a.seed, a.start, a.stop)


if __name__ == "__main__":
    main()
