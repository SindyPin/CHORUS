"""Concordance-guided training: does choosing (or weighting) training cohorts by their agreement with a small labelled
calibration part of the target cohort help?

For each block (BLOCKS) and each held-out cohort of its leave one study out (train_only cohorts are never held out):
20 repeated stratified 50/50 splits of the held-out cohort into a calibration part and a test part (sklearn
train_test_split, stratify, random_state = repeat). A cohort whose calibration part has fewer than 4 samples per class is
skipped and listed in skipped.csv. On the calibration part only, Hedges g of LC minus COVID_recovered is computed per
feature and its Spearman agreement with the g of each training cohort (all of its samples) over the block's features.
Arms, each fitted with the grid of this analysis and scored on the test part. Grid (reduced from the version 5 grid for
run time, as decided on 18 November 2026): harmonisers identity, center_scale and combat_addon (no nested, fsqn, rank or
chorus); selectors invariant and differential with the block config's k; logreg and rf (200 trees, one core per fit; the
benchmark uses 500 trees).

  1 all          all training cohorts (the leave one study out)
  2 guided       the training cohorts with agreement above 0
  3 weighted     the cohorts of 2, classifier sample weights proportional to max(0, agreement) of the sample's cohort
                 (harmoniser and selector unweighted)
  4 random       as many training cohorts as in 2, chosen at random (5 draws per split, seed (repeat, draw))
  5 calib        the calibration part alone
  6 calib_all    the calibration part (as its own study) plus all training cohorts
  7 calib_guided the calibration part plus the cohorts of 2

The harmoniser is fitted on the training studies of the arm (with labels only where the operator uses them, as in the
benchmark) and applied to the test part as a new study without labels (addon). Arms 2, 3, 4 and 7 are undefined when no
training cohort agrees (empty guided selection), and an arm with a single training study cannot fit the harmonisers that
need two (recorded as errors, as in the benchmark). Sample_matrix and other covariates of the block config are appended
as in the benchmark.

Outputs (--out): units/<block>__<cohort>.csv (one row per split, arm, draw and configuration: auroc or error),
agreement.csv (per split and training cohort), skipped.csv, and with `summarise`: per_cohort.csv (mean AUROC per arm over
splits and configurations, arm 4 first averaged over its draws), arms.csv (per block and arm: mean over held-out cohorts
with a 2,000 resample cohort bootstrap interval), differences.csv (paired differences 2-1, 2-4, 7-5 per block: per cohort
the mean over the splits where both arms are defined, then the mean over cohorts with the cohort bootstrap interval),
empty.csv (how often the guided selection is empty) and bars.csv.

  python scripts/concordance_guided.py all --out results/concordance_guided --workers 10   (every unit, smallest blocks first,
      in parallel processes; each block's summary is written as soon as its units finish, and the overall summary at the end)
  python scripts/concordance_guided.py run --block proteomics_affinity_v5 --held-out Su2022 --out results/concordance_guided
  python scripts/concordance_guided.py summarise --out results/concordance_guided [--block metabolomics]
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

from chorus.cli import load_config
from chorus.eval.loso import _covariate_matrix, _fit_harmoniser, _forced_extra, selector_label
from chorus.features import make_selector
from chorus.features.invariant import _hedges_g
from chorus.models import make_classifier
from chorus.schema import binary_label

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = {
    "transcriptomics_v4": "configs/lc_transcriptomics_blood_v4.yaml",
    "proteomics_affinity_v4": "configs/lc_proteomics_affinity_v4.yaml",
    "proteomics_ms_v4": "configs/lc_proteomics_ms_v4.yaml",
    "metabolomics": "configs/lc_metabolomics.yaml",
    "proteomics_affinity_v5": "configs/lc_proteomics_affinity_v5.yaml",
    "proteomics_explore_v5": "configs/lc_proteomics_explore_v5.yaml",
}
BLOCK_ORDER = ["metabolomics", "proteomics_ms_v4", "proteomics_affinity_v5", "proteomics_affinity_v4", "proteomics_explore_v5", "transcriptomics_v4"]  # smallest first
HARMONISERS = ["identity", "center_scale", "combat_addon"]
NESTED = False
CLASSIFIERS = ["logreg", "rf"]
RF_KWARGS = {"n_estimators": 200, "n_jobs": 1}
ARMS = ["all", "guided", "weighted", "random", "calib", "calib_all", "calib_guided"]
N_SPLITS, N_DRAWS, MIN_CAL_PER_CLASS = 20, 5, 4


def load_block(block: str) -> dict:
    cs, lc, cfg = load_config(str(ROOT / BLOCKS[block]))
    k = {s["name"]: s.get("k") for s in lc.selectors}
    lc.selectors = [{"name": n, "k": k[n]} for n in ("invariant", "differential")]
    lc.harmonisers, lc.classifiers = HARMONISERS, CLASSIFIERS
    X, M = cs.pooled()
    y = binary_label(M, lc.positive, lc.negative)
    keep = y.notna()
    X, M, y = X[keep], M[keep], y[keep].astype(int)
    held = [c.study_id for c in cs.cohorts if c.source.get("train_only") != "true"]
    return {"X": X, "M": M, "y": y, "s": M["study_id"].astype(str), "lc": lc, "held": held,
            "cov": _covariate_matrix(M, lc.covariates), "prot": _covariate_matrix(M, lc.protect) if lc.protect else None}


def hedges(X: pd.DataFrame, y: pd.Series) -> np.ndarray:
    with np.errstate(all="ignore"):
        return _hedges_g(X[(y == 1).to_numpy()].to_numpy(float), X[(y == 0).to_numpy()].to_numpy(float))[0]


def agreement(D: dict, cal: pd.Index, train_cohorts: list) -> pd.Series:
    gc = hedges(D["X"].loc[cal], D["y"].loc[cal])
    out = {}
    for c in train_cohorts:
        m = (D["s"] == c).to_numpy()
        gt = hedges(D["X"][m], D["y"][m])
        ok = np.isfinite(gc) & np.isfinite(gt)
        out[c] = float(stats.spearmanr(gc[ok], gt[ok]).statistic) if ok.sum() >= 3 else np.nan
    return pd.Series(out)


def fit_grid(D: dict, train_idx: pd.Index, studies: pd.Series, test_idx: pd.Index, test_study: str, weights: pd.Series = None) -> list:
    """Fit every configuration of the reduced grid on ``train_idx`` (study ids ``studies``) and score ``test_idx``."""
    lc = D["lc"]
    X_tr, X_te = D["X"].loc[train_idx], D["X"].loc[test_idx]
    y_tr, y_te = D["y"].loc[train_idx], D["y"].loc[test_idx]
    s_tr, s_te = studies.loc[train_idx], pd.Series(test_study, index=test_idx)
    prot = None if D["prot"] is None else D["prot"].loc[train_idx]
    rows = []
    for sel_cfg in lc.selectors:
        for h in HARMONISERS:
            chosen = h
            try:
                Xh_tr, Xh_te, _, _ = _fit_harmoniser(chosen, X_tr, s_tr, y_tr, X_te, s_te, prot, lc.harmoniser_params.get(chosen))
                sel = make_selector(**sel_cfg).fit(Xh_tr, y_tr, s_tr)
            except Exception as e:
                rows += [{"harmoniser": h, "harmoniser_chosen": chosen, "selector": selector_label(sel_cfg), "classifier": c, "auroc": np.nan, "error": str(e)[:150]} for c in CLASSIFIERS]
                continue
            forced = _forced_extra(lc, sel.selected_, Xh_tr.columns)
            A = pd.concat([sel.transform(Xh_tr), Xh_tr[forced], D["cov"].loc[train_idx]], axis=1)
            B = pd.concat([sel.transform(Xh_te), Xh_te.reindex(columns=forced), D["cov"].loc[test_idx]], axis=1)
            for c in CLASSIFIERS:
                try:
                    clf = make_classifier(c, lc.random_state, **(RF_KWARGS if c == "rf" else {}))
                    if weights is None:
                        clf.fit(A, y_tr)
                    else:
                        clf.fit(A, y_tr, clf__sample_weight=weights.loc[train_idx].to_numpy())
                    p = clf.predict_proba(B)[:, 1]
                    rows.append({"harmoniser": h, "harmoniser_chosen": chosen, "selector": selector_label(sel_cfg), "classifier": c, "auroc": roc_auc_score(y_te, p), "error": None})
                except Exception as e:
                    rows.append({"harmoniser": h, "harmoniser_chosen": chosen, "selector": selector_label(sel_cfg), "classifier": c, "auroc": np.nan, "error": str(e)[:150]})
    return rows


def run_unit(block: str, held: str, out: Path) -> None:
    D = load_block(block)
    s, y = D["s"], D["y"]
    target = s.index[s == held]
    train_cohorts = [c for c in pd.unique(s) if c != held]
    rows, agree = [], []
    t0 = time.time()
    for r in range(N_SPLITS):
        cal, test = train_test_split(target, test_size=0.5, stratify=y.loc[target], random_state=r)
        cal, test = pd.Index(cal), pd.Index(test)
        n_cal = y.loc[cal].value_counts()
        if r == 0 and min(n_cal.get(0, 0), n_cal.get(1, 0)) < MIN_CAL_PER_CLASS:
            pd.DataFrame([{"block": block, "held_out": held, "n_cal_LC": int(n_cal.get(1, 0)), "n_cal_recovered": int(n_cal.get(0, 0)),
                           "reason": f"calibration part has fewer than {MIN_CAL_PER_CLASS} samples per class"}]).to_csv(out / "units" / f"{block}__{held}__skipped.csv", index=False)
            print(f"{block} {held}: skipped", flush=True)
            return
        ag = agreement(D, cal, train_cohorts)
        agree += [{"block": block, "held_out": held, "split": r, "training_cohort": c, "agreement": v} for c, v in ag.items()]
        guided = [c for c in train_cohorts if ag.get(c, np.nan) > 0]
        studies = s.copy()
        studies.loc[cal] = f"{held}__calibration"
        idx_of = lambda cs_: s.index[s.isin(cs_)]  # noqa: E731
        arms = [("all", 0, idx_of(train_cohorts), None, len(train_cohorts))]
        if guided:
            w = s.map(lambda c: max(0.0, ag.get(c, 0.0)))
            arms += [("guided", 0, idx_of(guided), None, len(guided)), ("weighted", 0, idx_of(guided), w, len(guided))]
            for d in range(N_DRAWS):
                pick = list(np.random.default_rng([r, d]).choice(train_cohorts, len(guided), replace=False))
                arms.append(("random", d, idx_of(pick), None, len(pick)))
            arms.append(("calib_guided", 0, idx_of(guided).append(cal), None, len(guided) + 1))
        arms += [("calib", 0, cal, None, 1), ("calib_all", 0, idx_of(train_cohorts).append(cal), None, len(train_cohorts) + 1)]
        for arm, draw, tr, w, n_st in arms:
            for row in fit_grid(D, tr, studies, test, f"{held}__test", w):
                rows.append({"block": block, "held_out": held, "split": r, "arm": arm, "draw": draw, "n_guided": len(guided), "n_training_studies": n_st, **row})
        if not guided:
            rows.append({"block": block, "held_out": held, "split": r, "arm": "guided_empty", "draw": 0, "n_guided": 0})
        print(f"{block} {held} split {r}: {len(guided)} guided, {time.time() - t0:.0f} s", flush=True)
    pd.DataFrame(rows).to_csv(out / "units" / f"{block}__{held}.csv", index=False)
    pd.DataFrame(agree).to_csv(out / "units" / f"{block}__{held}__agreement.csv", index=False)


# ------------------------------------------------------------------------------------------------ summary
def boot_ci(v: np.ndarray, n_boot: int = 2000, seed: int = 0) -> tuple:
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    b = [rng.choice(v, len(v)).mean() for _ in range(n_boot)]
    return float(v.mean()), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def _true(v) -> bool:
    """True only for a true boolean (numpy or Python); None and NaN (one cohort, no interval) are not true."""
    return isinstance(v, (bool, np.bool_)) and bool(v)


def summarise(out: Path, block: str = None) -> None:
    """Summary files of every finished unit (or of one block, into --out/blocks/<block>/)."""
    pat = f"{block}__*" if block else "*"
    units = [f for f in sorted((out / "units").glob(f"{pat}.csv")) if not f.name.endswith(("__agreement.csv", "__skipped.csv"))]
    if block:
        out = out / "blocks" / block
        out.mkdir(parents=True, exist_ok=True)
    if not units:
        return
    U = pd.concat([pd.read_csv(f) for f in units], ignore_index=True)
    src = out if not block else out.parents[1]
    sk = [pd.read_csv(f) for f in sorted((src / "units").glob(f"{pat}__skipped.csv"))]
    (pd.concat(sk, ignore_index=True) if sk else pd.DataFrame(columns=["block", "held_out", "reason"])).to_csv(out / "skipped.csv", index=False)
    ag = [pd.read_csv(f) for f in sorted((src / "units").glob(f"{pat}__agreement.csv"))]
    pd.concat(ag, ignore_index=True).to_csv(out / "agreement.csv", index=False)
    E = U[U["arm"] == "guided_empty"]
    S = U[U["arm"].isin(ARMS)]
    split_counts = S.groupby(["block", "held_out"])["split"].nunique()
    empty = E.groupby(["block", "held_out"]).size().reindex(split_counts.index, fill_value=0)
    pd.DataFrame({"n_splits": split_counts, "n_empty_guided": empty, "fraction_empty": empty / split_counts}).reset_index().to_csv(out / "empty.csv", index=False)
    # per split: mean AUROC over configurations (arm 4: over draws and configurations)
    per_split = S.groupby(["block", "held_out", "arm", "split"])["auroc"].mean().rename("auroc").reset_index()
    per_split.to_csv(out / "per_split.csv", index=False)
    pc = per_split.groupby(["block", "held_out", "arm"])["auroc"].agg(["mean", "count"]).rename(columns={"mean": "auroc", "count": "n_splits"}).reset_index()
    pc.to_csv(out / "per_cohort.csv", index=False)
    arms = []
    for (b, arm), d in pc.groupby(["block", "arm"]):
        m, lo, hi = boot_ci(d["auroc"].to_numpy())
        arms.append({"block": b, "arm": arm, "n_cohorts": int(d["auroc"].notna().sum()), "mean_auroc": m, "ci_lo": lo, "ci_hi": hi})
    A = pd.DataFrame(arms)
    A["arm"] = pd.Categorical(A["arm"], ARMS)
    A.sort_values(["block", "arm"]).to_csv(out / "arms.csv", index=False)
    W = per_split.pivot_table(index=["block", "held_out", "split"], columns="arm", values="auroc")
    diffs = []
    for a, b in [("guided", "all"), ("guided", "random"), ("calib_guided", "calib")]:
        if a not in W or b not in W:
            continue
        d = (W[a] - W[b]).dropna().groupby(level=["block", "held_out"]).mean()
        for blk, v in d.groupby(level="block"):
            m, lo, hi = boot_ci(v.to_numpy())
            diffs.append({"block": blk, "comparison": f"{a} - {b}", "n_cohorts": len(v), "mean_difference": m, "ci_lo": lo, "ci_hi": hi, "interval_above_zero": bool(lo > 0) if len(v) >= 2 else None})  # one cohort: no interval
    Dd = pd.DataFrame(diffs)
    Dd.to_csv(out / "differences.csv", index=False)
    piv = Dd.set_index(["block", "comparison"])
    bars = []
    blocks = sorted(Dd["block"].unique())
    useful = [b for b in blocks if all(_true(piv.loc[(b, c), "interval_above_zero"]) for c in ("guided - all", "guided - random") if (b, c) in piv.index) and all((b, c) in piv.index for c in ("guided - all", "guided - random"))]
    local = [b for b in blocks if (b, "calib_guided - calib") in piv.index and piv.loc[(b, "calib_guided - calib"), "mean_difference"] > 0]
    local_ci = [b for b in blocks if (b, "calib_guided - calib") in piv.index and _true(piv.loc[(b, "calib_guided - calib"), "interval_above_zero"])]
    bars.append({"bar": "guided useful: arm 2 beats arms 1 and 4 with the paired interval above zero, in at least two blocks", "blocks_meeting": ";".join(useful), "n_blocks": len(useful), "met": len(useful) >= 2})
    bars.append({"bar": "adds value over local data: arm 7 beats arm 5 (mean difference above zero), in at least two blocks", "blocks_meeting": ";".join(local), "n_blocks": len(local), "met": len(local) >= 2})
    bars.append({"bar": "same, with the paired interval above zero", "blocks_meeting": ";".join(local_ci), "n_blocks": len(local_ci), "met": len(local_ci) >= 2})
    pd.DataFrame(bars).to_csv(out / "bars.csv", index=False)
    print(A.round(3).to_string(index=False)); print(Dd.round(3).to_string(index=False)); print(pd.DataFrame(bars).to_string(index=False))


def _run_unit(args):
    block, held, out = args
    run_unit(block, held, Path(out))
    return block, held


def run_all(out: Path, workers: int) -> None:
    """Every (block, held-out cohort) unit in parallel processes, smallest blocks first; a block's summary is written
    when its last unit finishes, the overall summary at the end."""
    import json
    from concurrent.futures import ProcessPoolExecutor, as_completed

    units = [(b, h, str(out)) for b in BLOCK_ORDER for h in load_block(b)["held"]]
    left = {b: sum(1 for u in units if u[0] == b) for b in BLOCK_ORDER}
    (out / "run_info.json").write_text(json.dumps({"blocks": {b: BLOCKS[b] for b in BLOCK_ORDER}, "harmonisers": HARMONISERS, "nested": NESTED, "selectors": ["invariant", "differential"],
                                                    "classifiers": CLASSIFIERS, "rf": RF_KWARGS, "n_splits": N_SPLITS, "n_random_draws": N_DRAWS, "min_calibration_per_class": MIN_CAL_PER_CLASS,
                                                    "units": [u[:2] for u in units], "workers": workers}, indent=1))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_run_unit, u) for u in units]
        for f in as_completed(futs):
            b, h = f.result()
            left[b] -= 1
            print(f"done {b} {h} ({left[b]} left in block)", flush=True)
            if left[b] == 0:
                summarise(out, b)
                print(f"block {b} summarised", flush=True)
    summarise(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["run", "all", "list", "summarise"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--out", required=True)
    ap.add_argument("--block")
    ap.add_argument("--held-out")
    a = ap.parse_args()
    out = Path(a.out)
    (out / "units").mkdir(parents=True, exist_ok=True)
    if a.command == "run":
        run_unit(a.block, a.held_out, out)
    elif a.command == "all":
        run_all(out, a.workers)
    elif a.command == "list":  # block and held-out cohort of every unit, for a parallel runner
        for b in BLOCKS:
            for h in load_block(b)["held"]:
                print(b, h)
    else:
        summarise(out, a.block)


if __name__ == "__main__":
    main()
