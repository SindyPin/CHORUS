"""Validation of the noise ceiling as a method (results/noise_ceiling_validation).

The pair procedure is that of scripts/concordance_noise_ceiling.py, unchanged:
  - Hedges g per feature in each cohort.
  - Split half reliability: 200 stratified half splits, Spearman-Brown corrected, on the features the pair shares.
  - Observed correlation of the two cohorts' g over the shared features.
  - One population expectation: the larger cohort is split into a stratified part with the class counts of the smaller
    cohort (capped at half of each class) and the rest, 200 times. The mean and the 2.5th and 97.5th percentiles of the
    two parts' correlation are kept.
  - Random streams as in that script, so its Spearman results are reproduced.
  - Every correlation is computed twice: Spearman (primary, as before) and Pearson (secondary), on the same splits.

Verdict (concordance_noise_ceiling.verdict): "too small to tell" if either reliability is at or below the threshold
(0.3; 0.2, 0.4 and 0.5 as sensitivity), else "agree" if observed >= the 2.5th percentile, else "truly differ".
Every pair is marked same or different platform (PLATFORM).

Sets:
  1. verdicts.csv     every pair of every set below, with verdicts for each measure and threshold
  2. pseudo_cohorts   large cohorts split into two pseudo cohorts with no shared participant (n per class each, n 10,
                      20 and 40 where each class has at least 2n; 200 repeats); LC: Buggert2025, Fineschi2026,
                      Woodruff2023, Cervia2024_SomaScan, Su2022 (concordance_proteomics_all_v5); acute against never
                      infected: Su2022 (control_acute_proteomics)
  3. controls         inferred sex on the RNA cohorts of results/control_sex_label (all genes; without chrX and chrY)
                      and acute against never infected (control_acute_proteomics), with the cohorts, samples and
                      features of the adaptive_borrowing tasks (concordance_guided.load_task)
  4. shuffled         LC labels shuffled within each cohort (class counts kept), 50 repeats per layer, layers of
                      concordance_noise_ceiling
  5. lc_verdicts      the verdict applied to results/concordance_noise_ceiling/pairs.csv (Spearman), with the Pearson
                      columns of the recomputation
  6. benefit          for every task and held-out cohort of results/adaptive_borrowing: arm 6 minus arm 5 at 8 per
                      class and half, against predictors from that task's pairs (same cohorts and features as the arms):
                      (i) share of training cohorts with verdict agree, (ii) mean observed / expected mean, (iii) mean
                      observed correlation; Spearman with a cohort bootstrap (held-out cohorts resampled)

  python scripts/noise_ceiling_validation.py --out results/noise_ceiling_validation --workers 10
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concordance_guided as cg  # noqa: E402
import concordance_noise_ceiling as nc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MEASURES = ("spearman", "pearson")
THRESHOLDS = (0.2, 0.3, 0.4, 0.5)
PRIMARY = 0.3
N_SPLITS, N_PSEUDO, N_SHUFFLE, N_BOOT = 200, 200, 50, 2000
PSEUDO_SIZES = (10, 20, 40)
PSEUDO = {"LC": ("configs/concordance_proteomics_all_v5.yaml", ["Buggert2025", "Fineschi2026", "Woodruff2023", "Cervia2024_SomaScan", "Su2022"]),
          "acute": ("configs/control_acute_proteomics.yaml", ["Su2022"])}
CONTROL_TASKS = ["sex_all_genes", "sex_no_xy", "acute_proteomics"]
LC_TASKS = list(cg.BLOCKS)
# platform, as the technology that measured the cohort: Olink Explore (NGS readout) and Olink Target 96 (qPCR readout)
# are kept apart, mass spectrometry proteomics is one platform, and each metabolomics cohort is its own platform
PLATFORM = {
    **{c: "bulk RNA-seq" for c in ["GSE226260_WB", "GSE226260_PBMC", "GSE224615", "GSE251849", "GSE270045", "PRJNA1149196_SA", "PRJNA1149196_QLD", "GSE297499", "GSE228320"]},
    "CORSIP2026": "NanoString",
    **{c: "Olink Explore" for c in ["Buggert2025", "Fineschi2026", "Woodruff2023"]},
    **{c: "Olink Target 96" for c in ["Su2022", "Bansal2025", "Zoodsma2022_Hannover"]},
    "Cervia2024_SomaScan": "SomaScan",
    **{c: "MS proteomics" for c in ["Cervia2024_MS", "Chowdhury2024", "Chowdhury2026", "Kovarik2023_plasma"]},
    "Martinez2024": "CE/GC/LC-MS", "Lopez2023": "targeted LC-MS", "Saito2024": "CIL LC-MS", "Su2022_metabolomics": "untargeted LC-MS",
}


# ------------------------------------------------------------------------------------------------ pair engine
def corr(a: np.ndarray, b: np.ndarray, kind: str) -> float:
    if kind == "spearman":
        return nc.spearman(a, b)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(stats.pearsonr(a[ok], b[ok]).statistic)


def _mean(v: np.ndarray) -> float:
    return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan


def run_set(C: dict, n_splits: int = N_SPLITS, seed: int = 0, label: str = "") -> pd.DataFrame:
    """Every pair of the cohorts ``C`` ({name: (X DataFrame, y array)}), as concordance_noise_ceiling.run_layer, for both
    measures."""
    names = list(C)
    cols = {n: C[n][0].columns for n in names}
    halves = {}
    for i, n in enumerate(names):
        X, y = C[n][0].to_numpy(), C[n][1]
        rng = np.random.default_rng([seed, 0, i])
        halves[n] = [(nc.g_of(X, y, a), nc.g_of(X, y, b)) for a, b in (nc.half_split(y, rng) for _ in range(n_splits))]
    rows = []
    for (i, a), (j, b) in itertools.combinations(enumerate(names), 2):
        shared = cols[a].intersection(cols[b])
        Xa, ya = C[a][0][shared].to_numpy(), C[a][1]
        Xb, yb = C[b][0][shared].to_numpy(), C[b][1]
        row = {"set": label, "cohort_a": a, "cohort_b": b, "platform_a": PLATFORM.get(a, "?"), "platform_b": PLATFORM.get(b, "?"),
               "same_platform": PLATFORM.get(a, "?a") == PLATFORM.get(b, "?b"), "n_shared": len(shared), "n_a": len(ya), "n_b": len(yb)}
        if len(shared) < 3:
            rows.append({**row, **{f"observed_{m}": np.nan for m in MEASURES}})
            continue
        ga, gb = nc.g_of(Xa, ya, np.arange(len(ya))), nc.g_of(Xb, yb, np.arange(len(yb)))
        pos = {a: cols[a].get_indexer(shared), b: cols[b].get_indexer(shared)}
        big = a if len(ya) >= len(yb) else b
        Xg, yg = (Xa, ya) if big == a else (Xb, yb)
        ys = yb if big == a else ya
        counts = {0: int((ys == 0).sum()), 1: int((ys == 1).sum())}
        rng = np.random.default_rng([seed, 1, i, j])
        parts = [nc.sub_split(yg, counts, rng) for _ in range(n_splits)]
        draws = [(nc.g_of(Xg, yg, p), nc.g_of(Xg, yg, q)) for p, q in parts]
        row["larger_cohort"] = big
        for m in MEASURES:
            obs = corr(ga, gb, m)
            sh = {n: _mean(np.array([corr(h1[pos[n]], h2[pos[n]], m) for h1, h2 in halves[n]])) for n in (a, b)}
            rel_a, rel_b = nc.sb(sh[a]), nc.sb(sh[b])
            exp = np.array([corr(p, q, m) for p, q in draws])
            ok = exp[np.isfinite(exp)]
            row.update({f"observed_{m}": obs, f"split_half_a_{m}": sh[a], f"split_half_b_{m}": sh[b], f"reliability_a_{m}": rel_a, f"reliability_b_{m}": rel_b,
                        f"expected_mean_{m}": ok.mean() if len(ok) else np.nan, f"expected_p2_5_{m}": np.percentile(ok, 2.5) if len(ok) else np.nan,
                        f"expected_p97_5_{m}": np.percentile(ok, 97.5) if len(ok) else np.nan, f"n_expected_draws_{m}": len(ok)})
        rows.append(row)
    return add_verdicts(pd.DataFrame(rows))


def add_verdicts(P: pd.DataFrame) -> pd.DataFrame:
    for m in MEASURES:
        for c in (f"reliability_a_{m}", f"reliability_b_{m}", f"observed_{m}", f"expected_p2_5_{m}", f"expected_mean_{m}"):
            if c not in P:
                P[c] = np.nan
        for t in THRESHOLDS:
            P[f"verdict_{m}_r{t}"] = [nc.verdict(ra, rb, o, e, t) for ra, rb, o, e in
                                      P[[f"reliability_a_{m}", f"reliability_b_{m}", f"observed_{m}", f"expected_p2_5_{m}"]].itertuples(index=False)]
    P["verdict"] = P[f"verdict_spearman_r{PRIMARY}"]
    return P


# ------------------------------------------------------------------------------------------------ data
def cohorts_of_task(task: str) -> dict:
    """The cohorts of an adaptive_borrowing task, each on its non-missing features (the arms' data)."""
    D = cg.load_task(task)
    out = {}
    for c in pd.unique(D["s"]):
        m = (D["s"] == c).to_numpy()
        X = D["X"][m].astype(float)
        out[c] = (X.loc[:, X.notna().any(axis=0)], D["y"][m].to_numpy().astype(int))
    return out


def pseudo_split(y: np.ndarray, n: int, rng) -> tuple | None:
    """Two disjoint stratified pseudo cohorts with ``n`` samples of each class (one sample per participant in these
    cohorts, so no participant is shared); None when a class has fewer than 2n samples."""
    if min((y == 0).sum(), (y == 1).sum()) < 2 * n:
        return None
    a, b = [], []
    for k in (0, 1):
        idx = rng.permutation(np.where(y == k)[0])
        a.append(idx[:n])
        b.append(idx[n:2 * n])
    return np.sort(np.concatenate(a)), np.sort(np.concatenate(b))


def shuffle_within(C: dict, rng) -> dict:
    return {n: (X, rng.permutation(y)) for n, (X, y) in C.items()}


# ------------------------------------------------------------------------------------------------ jobs
_DATA: dict = {}


def _init(data):
    _DATA.update(data)


def job(kind: str, *args) -> pd.DataFrame:
    t0 = time.time()
    if kind == "layer":  # LC layer of concordance_noise_ceiling, same seed
        layer, = args
        out = run_set(_DATA[f"layer:{layer}"], label=f"lc:{layer}").assign(layer=layer)
    elif kind == "task":
        task, = args
        out = run_set(cohorts_of_task(task), label=f"task:{task}").assign(task=task)
    elif kind == "shuffle":
        layer, r = args
        C = shuffle_within(_DATA[f"layer:{layer}"], np.random.default_rng([7, r]))
        out = run_set(C, seed=1000 + r, label=f"shuffled:{layer}").assign(layer=layer, repeat=r)
    elif kind == "pseudo":
        contrast, cohort, n = args
        X, y = _DATA[f"pseudo:{contrast}"][cohort]
        rows = []
        for r in range(N_PSEUDO):
            a, b = pseudo_split(y, n, np.random.default_rng([11, n, r]))
            P = run_set({f"{cohort}_A": (X.iloc[a], y[a]), f"{cohort}_B": (X.iloc[b], y[b])}, seed=2000 + r, label=f"pseudo:{contrast}")
            rows.append(P.assign(contrast=contrast, cohort=cohort, size=n, repeat=r, same_platform=True))
        out = pd.concat(rows, ignore_index=True)
    print(f"done {kind} {args} in {time.time() - t0:.0f} s", flush=True)
    return out


# ------------------------------------------------------------------------------------------------ summaries
def tally(P: pd.DataFrame, by: list, col: str = "verdict") -> pd.DataFrame:
    t = P.groupby(by)[col].value_counts().unstack(fill_value=0).reindex(columns=list(nc.VERDICTS), fill_value=0)
    t["n"] = t.sum(axis=1)
    t["evaluable"] = t["agree"] + t["truly differ"]
    for v in ("agree", "truly differ", "too small to tell"):
        t[f"fraction_{v}"] = t[v] / t["n"]
    t["fraction_agree_of_evaluable"] = t["agree"] / t["evaluable"].replace(0, np.nan)
    t["fraction_differ_of_evaluable"] = t["truly differ"] / t["evaluable"].replace(0, np.nan)
    return t.reset_index()


def tally_all(P: pd.DataFrame, by: list, platform: bool = True) -> pd.DataFrame:
    """Tallies for each measure and threshold, overall and (``platform``) split by same and different platform."""
    out = []
    for m in MEASURES:
        for t in THRESHOLDS:
            col = f"verdict_{m}_r{t}"
            out.append(tally(P, by, col).assign(measure=m, threshold=t, platform_group="all"))
            if platform:
                for g, d in P.groupby("same_platform"):
                    out.append(tally(d, by, col).assign(measure=m, threshold=t, platform_group="same platform" if g else "different platform"))
    return pd.concat(out, ignore_index=True)


def spearman_boot(d: pd.DataFrame, x: str, y: str, n_boot: int, seed: int, x2: str = None) -> dict:
    """Spearman of x with y over units, with a cohort bootstrap (held-out cohorts resampled with all their units);
    with ``x2``, the difference rho(x) - rho(x2) on the same resamples."""
    d = d.dropna(subset=[x, y] + ([x2] if x2 else []))

    def rho(t, col):
        return stats.spearmanr(t[col], t[y]).statistic if len(t) >= 3 and t[col].nunique() > 1 and t[y].nunique() > 1 else np.nan

    est = rho(d, x) - (rho(d, x2) if x2 else 0)
    groups = {c: g for c, g in d.groupby("held_out")}
    names = list(groups)
    rng = np.random.default_rng(seed)
    b = []
    for _ in range(n_boot):
        t = pd.concat([groups[c] for c in rng.choice(names, len(names))])
        b.append(rho(t, x) - (rho(t, x2) if x2 else 0))
    b = np.array(b, float)
    b = b[np.isfinite(b)]
    lo, hi = (np.percentile(b, 2.5), np.percentile(b, 97.5)) if len(b) else (np.nan, np.nan)
    return {"n_units": len(d), "n_cohorts": len(names), "estimate": est, "ci_lo": lo, "ci_hi": hi, "interval_above_zero": bool(lo > 0) if np.isfinite(lo) else None, "n_valid_resamples": len(b)}


def predictors(T: pd.DataFrame) -> pd.DataFrame:
    """Per task and held-out cohort, predictors from its pairs with the training cohorts (all other cohorts of the task)."""
    rows = []
    for task, d in T.groupby("task"):
        for h in pd.unique(pd.concat([d["cohort_a"], d["cohort_b"]])):
            p = d[(d["cohort_a"] == h) | (d["cohort_b"] == h)]
            ratio = p["observed_spearman"] / p["expected_mean_spearman"]
            rows.append({"task": task, "held_out": h, "n_training": len(p), "n_agree": int((p["verdict"] == "agree").sum()),
                         "n_truly_differ": int((p["verdict"] == "truly differ").sum()), "n_too_small": int((p["verdict"] == "too small to tell").sum()),
                         "share_agree": (p["verdict"] == "agree").mean(), "mean_ratio_observed_expected": ratio[np.isfinite(ratio)].mean(),
                         "mean_observed": p["observed_spearman"].mean()})
    return pd.DataFrame(rows)


def benefit_table(T: pd.DataFrame, per_cohort: pd.DataFrame) -> tuple:
    pc = per_cohort.pivot_table(index=["task", "size", "held_out"], columns="arm", values="auroc").reset_index()
    pc["benefit"] = pc["calib_all"] - pc["calib"]
    B = pc[pc["size"].isin(["8", "half"])][["task", "size", "held_out", "calib", "calib_all", "benefit"]].merge(predictors(T), on=["task", "held_out"], how="left")
    B["group"] = np.where(B["task"].isin(LC_TASKS), "LC", "control")
    cor = []
    for size in ("8", "half"):
        for scope, d in (("all tasks", B[B["size"] == size]), ("LC", B[(B["size"] == size) & (B["group"] == "LC")])):
            for name, x in (("i share agree", "share_agree"), ("ii mean observed / expected", "mean_ratio_observed_expected"), ("iii mean observed (no ceiling)", "mean_observed")):
                cor.append({"size": size, "scope": scope, "predictor": name, **spearman_boot(d, x, "benefit", N_BOOT, 0)})
            for name, x in (("i - iii", "share_agree"), ("ii - iii", "mean_ratio_observed_expected")):
                cor.append({"size": size, "scope": scope, "predictor": name, **spearman_boot(d, x, "benefit", N_BOOT, 0, x2="mean_observed")})
    # mean benefit of the held-out cohort over its pairs, by the pair's verdict
    long = []
    for task, d in T.groupby("task"):
        for _, r in d.iterrows():
            for h, o in ((r["cohort_a"], r["cohort_b"]), (r["cohort_b"], r["cohort_a"])):
                long.append({"task": task, "held_out": h, "training_cohort": o, "verdict": r["verdict"]})
    L = pd.DataFrame(long).merge(B[["task", "size", "held_out", "benefit", "group"]], on=["task", "held_out"])
    byv = L.groupby(["size", "group", "verdict"])["benefit"].agg(["mean", "count"]).reset_index()
    byv = pd.concat([byv, L.groupby(["size", "verdict"])["benefit"].agg(["mean", "count"]).reset_index().assign(group="all tasks")], ignore_index=True)
    return B, pd.DataFrame(cor), byv


def bars(pseudo_t: pd.DataFrame, controls_t: pd.DataFrame, shuffled_t: pd.DataFrame, cor: pd.DataFrame) -> pd.DataFrame:
    sel = lambda t: t[(t["measure"] == "spearman") & (t["threshold"] == PRIMARY) & (t["platform_group"] == "all")]  # noqa: E731
    p = sel(pseudo_t)
    a_ok = bool((p["fraction_truly differ"] <= 0.05).all())
    c = sel(controls_t).set_index("set")
    b_rows = {s: c.loc[f"task:{s}", "fraction_agree_of_evaluable"] for s in ("sex_all_genes", "acute_proteomics")}
    b_ok = all(np.isfinite(v) and v >= 0.80 for v in b_rows.values())
    s = sel(shuffled_t)
    c_ok = bool((s["fraction_truly differ"] <= 0.05).all())
    k = cor[(cor["size"] == "8") & (cor["scope"] == "all tasks")].set_index("predictor")
    d_ok = any(k.loc[x, "estimate"] > 0 and cg._true(k.loc[x, "interval_above_zero"]) for x in ("i share agree", "ii mean observed / expected"))
    e_ok = any(k.loc[x, "estimate"] > k.loc["iii mean observed (no ceiling)", "estimate"] for x in ("i share agree", "ii mean observed / expected"))
    fmt = lambda x: f"{x:.3f}"  # noqa: E731
    return pd.DataFrame([
        {"bar": "a: pseudo cohorts, truly differ in at most 5% of splits at every size (each contrast and size)", "met": a_ok,
         "detail": "; ".join(f"{r.contrast} n={r['size']}: {fmt(r['fraction_truly differ'])}" for _, r in p.iterrows())},
        {"bar": "b: positive controls, at least 80% agree among pairs not too small to tell (sex on all genes, acute)", "met": b_ok,
         "detail": "; ".join(f"{k_}: {fmt(v)}" for k_, v in b_rows.items())},
        {"bar": "c: shuffled labels, truly differ in at most 5% of pairs (each layer)", "met": c_ok,
         "detail": "; ".join(f"{r.layer}: {fmt(r['fraction_truly differ'])}" for _, r in s.iterrows())},
        {"bar": "d: predictor i or ii has a positive Spearman with the benefit of combining, interval above zero, all tasks, 8 per class", "met": d_ok,
         "detail": "; ".join(f"{x}: {fmt(k.loc[x, 'estimate'])} ({fmt(k.loc[x, 'ci_lo'])} to {fmt(k.loc[x, 'ci_hi'])})" for x in ("i share agree", "ii mean observed / expected", "iii mean observed (no ceiling)"))},
        {"bar": "e: predictor i or ii has a higher Spearman than iii, all tasks, 8 per class (descriptive if the interval includes zero)", "met": e_ok,
         "detail": "; ".join(f"{x}: {fmt(k.loc[x, 'estimate'])} ({fmt(k.loc[x, 'ci_lo'])} to {fmt(k.loc[x, 'ci_hi'])})" for x in ("i - iii", "ii - iii"))},
    ])


def main():
    from concurrent.futures import ProcessPoolExecutor

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--bars-only", action="store_true", help="recompute bars.csv from the saved summary files")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.bars_only:
        pt = pd.read_csv(out / "pseudo_cohorts.csv")
        Bars = bars(pt[pt["cohort"] == "all"], pd.read_csv(out / "controls.csv"), pd.read_csv(out / "shuffled.csv"), pd.read_csv(out / "benefit_correlation.csv", dtype={"size": str}))
        Bars.to_csv(out / "bars.csv", index=False)
        print(Bars.to_string(index=False))
        return
    t0 = time.time()
    data = {f"layer:{l}": nc.cohorts_of(c) for l, c in nc.LAYERS.items()}
    for k, (cfg, names) in PSEUDO.items():
        C = data["layer:proteomics"] if k == "LC" else nc.cohorts_of(cfg)
        data[f"pseudo:{k}"] = {n: C[n] for n in names}
    for k, C in data.items():
        assert "PRJNA1184005" not in C, k  # analysis set frozen at version 5
    print(f"data loaded in {time.time() - t0:.0f} s", flush=True)
    jobs = [("task", t) for t in ["sex_all_genes", "sex_no_xy", "transcriptomics_v4"]] + [("layer", "transcriptomics")]
    jobs += [("shuffle", "transcriptomics", r) for r in range(N_SHUFFLE)]
    jobs += [("task", t) for t in cg.TASKS if t not in ("sex_all_genes", "sex_no_xy", "transcriptomics_v4")] + [("layer", "proteomics"), ("layer", "metabolomics")]
    jobs += [("pseudo", k, c, n) for k, (_, names) in PSEUDO.items() for c in names for n in PSEUDO_SIZES
             if pseudo_split(data[f"pseudo:{k}"][c][1], n, np.random.default_rng(0)) is not None]
    jobs += [("shuffle", l, r) for l in ("proteomics", "metabolomics") for r in range(N_SHUFFLE)]
    with ProcessPoolExecutor(max_workers=a.workers, initializer=_init, initargs=(data,)) as ex:
        futs = {j: ex.submit(job, *j) for j in jobs}
        res = {j: f.result() for j, f in futs.items()}
    by = lambda kind: pd.concat([v for k, v in res.items() if k[0] == kind], ignore_index=True)  # noqa: E731
    LCr, T, S, Ps = by("layer"), by("task"), by("shuffle"), by("pseudo")
    V = pd.concat([LCr, T, S, Ps], ignore_index=True)
    V.to_csv(out / "verdicts.csv", index=False)
    # 2 pseudo cohorts
    Ps.to_csv(out / "pseudo_cohorts_splits.csv", index=False)
    pt = tally_all(Ps, ["contrast", "size"], platform=False)
    pd.concat([pt.assign(cohort="all"), tally_all(Ps, ["contrast", "cohort", "size"], platform=False)], ignore_index=True).to_csv(out / "pseudo_cohorts.csv", index=False)
    # 3 controls
    Ct = T[T["task"].isin(CONTROL_TASKS)]
    ct = tally_all(Ct, ["set"])
    Ct.to_csv(out / "controls_pairs.csv", index=False)
    ct.to_csv(out / "controls.csv", index=False)
    # 4 shuffled
    S.to_csv(out / "shuffled_pairs.csv", index=False)
    st = tally_all(S, ["layer"])
    st.to_csv(out / "shuffled.csv", index=False)
    # 5 LC: the verdict on the committed pairs.csv; Pearson from the recomputation (same splits)
    old = pd.read_csv(ROOT / "results/concordance_noise_ceiling/pairs.csv")
    old = old.rename(columns={"observed": "observed_spearman", "reliability_a_shared": "reliability_a_spearman", "reliability_b_shared": "reliability_b_spearman",
                              "expected_p2_5": "expected_p2_5_spearman", "expected_mean": "expected_mean_spearman"})
    rec = LCr.set_index(["layer", "cohort_a", "cohort_b"])
    keep = [c for c in rec.columns if c.endswith("_pearson") or c in ("platform_a", "platform_b", "same_platform")]
    L = old.merge(rec[keep].reset_index(), on=["layer", "cohort_a", "cohort_b"], how="left")
    L = add_verdicts(L)
    chk = old.set_index(["layer", "cohort_a", "cohort_b"])["observed_spearman"].sub(rec["observed_spearman"]).abs().max()
    chk_e = old.set_index(["layer", "cohort_a", "cohort_b"])["expected_p2_5_spearman"].sub(rec["expected_p2_5_spearman"]).abs().max()
    L.to_csv(out / "lc_verdicts.csv", index=False)
    lt = tally_all(L, ["layer"])
    lt.to_csv(out / "lc_verdicts_summary.csv", index=False)
    # 6 benefit
    pc = pd.read_csv(ROOT / "results/adaptive_borrowing/per_cohort.csv", dtype={"size": str})
    B, cor, byv = benefit_table(T, pc)
    B.to_csv(out / "benefit.csv", index=False)
    cor.to_csv(out / "benefit_correlation.csv", index=False)
    byv.to_csv(out / "benefit_by_verdict.csv", index=False)
    Bars = bars(pt, ct, st, cor)
    Bars.to_csv(out / "bars.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({
        "pair_procedure": "scripts/concordance_noise_ceiling.py (split half 200, expectation 200 draws), Spearman primary, Pearson secondary",
        "thresholds": THRESHOLDS, "primary_threshold": PRIMARY, "n_splits": N_SPLITS, "pseudo": {k: v[1] for k, v in PSEUDO.items()}, "pseudo_sizes": PSEUDO_SIZES,
        "n_pseudo_repeats": N_PSEUDO, "n_shuffle_repeats": N_SHUFFLE, "n_boot": N_BOOT, "layers": nc.LAYERS, "tasks": cg.TASKS, "platform": PLATFORM,
        "lc_recomputation_max_abs_difference": {"observed_spearman": float(chk), "expected_p2_5_spearman": float(chk_e)},
        "jobs": [list(map(str, j)) for j in jobs], "workers": a.workers, "runtime_s": round(time.time() - t0)}, indent=1))
    pd.set_option("display.width", 250)
    print("LC recomputation max abs difference: observed", chk, "p2.5", chk_e)
    sel = lambda t: t[(t["measure"] == "spearman") & (t["threshold"] == PRIMARY) & (t["platform_group"] == "all")]  # noqa: E731
    for name, t in (("pseudo", pt), ("controls", ct), ("shuffled", st), ("LC", lt)):
        print(name)
        print(sel(t).round(3).to_string(index=False))
    print(cor.round(3).to_string(index=False))
    print(Bars.to_string(index=False))


if __name__ == "__main__":
    main()
