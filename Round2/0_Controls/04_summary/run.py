"""Step 04: summary of the known answer controls (steps 01 to 03).

Every summary rests on the first 15 folds (5 folds by the first 3 of the 10 seeded repeats, Amendment 3), also for
jobs of step 01 that ran more folds before the amendment.

Writes, in output/:
  metrics_by_arm.csv      every control, layer, arm and model: mean and SD of every metric over the folds
  paired_differences.csv  every model and metric: paired difference between arms on the same folds (positive =
                          better), 95% interval from the corrected resampled t test, verdict
  counts.csv              per comparison and metric: models that improve, do not change, get worse
  failures.csv            every failure of steps 01 to 03
  answers.csv             the plain answer per control and layer (primary metric: AUROC for sex, MAE for age)
  models_not_run.csv      per model and reason: the controls, layers and arms where it has no result (Amendment 2, item 3)

  python Round2/0_Controls/04_summary/run.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as K  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
STEPS = ["01_one_cohort", "02_horizontal_one_plus_one", "03_vertical_one_plus_one"]
# (better candidate, reference): difference = candidate - reference, oriented so that positive is better
HORIZONTAL = [("A_shared", "A_all"), ("B_naive", "A_shared"), ("B_scaled", "A_shared"), ("B_scaled", "B_naive"),
              ("B_naive", "A_all"), ("B_scaled", "A_all")]
VERTICAL = [("joined", "layer_1"), ("joined", "layer_2"), ("late_fusion", "layer_1"), ("late_fusion", "layer_2"),
            ("late_fusion", "joined"), ("layer_2", "layer_1")]


def load() -> pd.DataFrame:
    rows = []
    for s in STEPS:
        for f in sorted((HERE.parent / s / "output" / "raw").glob("*.csv")):
            if not f.name.endswith("_pred.csv"):
                d = pd.read_csv(f)
                rows.append(d[d["fold"] < K.N_FOLDS])  # Amendment 3: the first 3 repeats only, also where 10 were run
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    R = load()
    R["family"] = np.where(R["step"].str.startswith("03"), "vertical", "horizontal")
    metric_cols = [m for m in K.CLASSIFICATION_METRICS + K.REGRESSION_METRICS if m in R]
    g = R.groupby(["family", "control", "layer", "arm", "model"])
    M = g[metric_cols].agg(["mean", "std"])
    M.columns = [f"{a}_{b}" for a, b in M.columns]
    M = M.join(g["fold"].nunique().rename("n_folds")).reset_index()
    M.to_csv(OUT / "metrics_by_arm.csv", index=False)

    D = []
    for (fam, control, layer), d in R.groupby(["family", "control", "layer"]):
        task = K.TASK[control]
        mets = K.CLASSIFICATION_METRICS if task == "classification" else K.REGRESSION_METRICS
        ref = d[d["arm"].isin(["A_all", "layer_1"])]
        ratio = (ref["n_test"] / ref["n_train"]).median()  # the anchor's (or participants') own split
        for cand, base in (HORIZONTAL if fam == "horizontal" else VERTICAL):
            for model, dm in d.groupby("model"):
                a = dm[dm["arm"] == cand].set_index("fold")
                b = dm[dm["arm"] == base].set_index("fold")
                common = a.index.intersection(b.index)
                for met in mets:
                    if len(common) < 2:
                        D.append(dict(family=fam, control=control, layer=layer, comparison=f"{cand} vs {base}", model=model, metric=met,
                                      n_folds=len(common), verdict="not evaluable"))
                        continue
                    diff = a.loc[common, met].to_numpy() - b.loc[common, met].to_numpy()
                    if met in K.LOWER_IS_BETTER:
                        diff = -diff
                    mean, lo, hi = K.corrected_t(diff, ratio, 1.0)
                    D.append(dict(family=fam, control=control, layer=layer, comparison=f"{cand} vs {base}", model=model, metric=met,
                                  n_folds=len(common), difference=mean, ci_low=lo, ci_high=hi, verdict=K.verdict(lo, hi)))
    D = pd.DataFrame(D)
    D.to_csv(OUT / "paired_differences.csv", index=False)
    Cn = (D.groupby(["family", "control", "layer", "comparison", "metric", "verdict"]).size().unstack("verdict", fill_value=0).reset_index())
    for v in ("improves", "does not change", "gets worse", "not evaluable"):
        if v not in Cn:
            Cn[v] = 0
    Cn = Cn[["family", "control", "layer", "comparison", "metric", "improves", "does not change", "gets worse", "not evaluable"]]
    Cn.to_csv(OUT / "counts.csv", index=False)

    F = [pd.read_csv(f) for s in STEPS if (f := HERE.parent / s / "output" / "failures.csv").exists()]
    F = pd.concat(F, ignore_index=True) if F else pd.DataFrame(columns=["task", "arm", "model", "reason", "message", "device", "memory"])
    if len(F):
        F["message"] = F["message"].astype(str).str.replace(r"\s*\n\s*", " ", regex=True)
    F.to_csv(OUT / "failures.csv", index=False)
    if len(F):
        F["where"] = F["task"].str.split("/").str[1:].str.join(" ") + " " + F["arm"]
        N = (F.groupby(["model", "reason"]).agg(jobs=("job", "nunique"), where=("where", lambda w: "; ".join(sorted(set(w)))),
                                                example_message=("message", "first")).reset_index())
        N.to_csv(OUT / "models_not_run.csv", index=False)

    P2 = pd.read_csv(HERE.parent / "02_horizontal_one_plus_one" / "output" / "pairs.csv")
    ans = []
    for _, c in Cn.iterrows():
        prim = K.PRIMARY[K.TASK[c["control"]]]
        if c["metric"] != prim:
            continue
        if c["family"] == "horizontal" and c["comparison"] not in ("B_scaled vs A_shared", "B_naive vs A_shared"):
            continue
        if c["family"] == "vertical" and c["comparison"] not in ("joined vs layer_1", "joined vs layer_2", "late_fusion vs layer_1", "late_fusion vs layer_2"):
            continue
        note = ""
        if c["comparison"] == "B_naive vs A_shared":
            ss = P2.loc[(P2["control"] == c["control"]) & (P2["layer"] == c["layer"]), "same_scale"]
            note = "interpretable (same scale)" if len(ss) and bool(ss.iloc[0]) else "reported only: different platforms or scales (Amendment 1, item 5)"
        n = c["improves"] + c["does not change"] + c["gets worse"]
        ans.append(dict(control=c["control"], layer=c["layer"], question="adding one cohort" if c["family"] == "horizontal" else "adding one layer",
                        comparison=c["comparison"], metric=prim, models_evaluated=n, improves=c["improves"],
                        does_not_change=c["does not change"], gets_worse=c["gets worse"], note=note))
    pd.DataFrame(ans).to_csv(OUT / "answers.csv", index=False)
    print(pd.DataFrame(ans).to_string(index=False))
    print(F[["task", "arm", "model", "reason", "message"]].to_string(index=False) if len(F) else "no failures")


if __name__ == "__main__":
    main()
