"""Step 05: summary of Round2/5_Clinical (steps 02 to 04; step 01 is summarised from its own tables).

Every summary rests on 15 folds (5 folds by the first 3 of the 10 seeded repeats, Amendment 3).

Writes, in output/:
  metrics_by_arm.csv      every step, task, unit or study, arm and model: mean and SD of every metric over the folds
  paired_differences.csv  every model and metric: paired difference between arms on the same folds (positive = better),
                          95% interval from the corrected resampled t test, verdict
  counts.csv              per comparison and metric: models that improve, do not change, get worse
  answers.csv             the comparisons that answer the README questions, on AUROC
  failures.csv            every failure of steps 02 to 04
  models_not_run.csv      per model and reason: jobs and where (Amendment 2, item 3)

  python Round2/5_Clinical/05_summary/run.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import clinical as Q  # noqa: E402

K = Q.K
HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
STEPS = ["02_clinical_alone", "03_horizontal_one_plus_one", "04_vertical"]
COMPARISONS = {  # step: (candidate, reference, question)
    "02_clinical_alone": [("clinical", "age_sex", "do the other clinical variables add to age and sex?"),
                          ("clinical_plus_timing", "clinical", "sensitivity: does time since infection add to the clinical model?")],
    "03_horizontal_one_plus_one": [("A_shared", "A_all", "shared variables against all of the anchor's"),
                                   ("B_naive", "A_shared", "does combining cohorts help the clinical model? (no adjustment)"),
                                   ("B_scaled", "A_shared", "does combining cohorts help the clinical model? (within cohort scaling)"),
                                   ("B_scaled", "B_naive", "scaling against none")],
    "04_vertical": [("joined", "clinical", "does omics add value beyond clinical data?"),
                    ("joined", "omics", "does clinical data add value beyond omics?"),
                    ("clinical", "omics", "does a clinical model beat the omics model?"),
                    ("late_fusion", "clinical", "late fusion against clinical alone"),
                    ("late_fusion", "omics", "late fusion against omics alone")],
}
REFERENCE_ARM = {"02_clinical_alone": "clinical", "03_horizontal_one_plus_one": "A_all", "04_vertical": "omics"}


def load() -> pd.DataFrame:
    rows = []
    for s in STEPS:
        for f in sorted((Q.HERE / s / "output" / "raw").glob("*.csv")):
            if not f.name.endswith("_pred.csv"):
                d = pd.read_csv(f)
                rows.append(d[d["fold"] < K.N_FOLDS])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    R = load()
    mets = K.CLASSIFICATION_METRICS
    g = R.groupby(["step", "control", "layer", "arm", "model"])
    M = g[mets].agg(["mean", "std"])
    M.columns = [f"{a}_{b}" for a, b in M.columns]
    M = M.join(g["fold"].nunique().rename("n_folds")).reset_index().rename(columns={"control": "task", "layer": "unit_or_study"})
    M.to_csv(OUT / "metrics_by_arm.csv", index=False)

    Dl = []
    for (step, task, layer), d in R.groupby(["step", "control", "layer"]):
        ref = d[d["arm"] == REFERENCE_ARM[step]]
        ratio = (ref["n_test"] / ref["n_train"]).median()
        for cand, base, question in COMPARISONS[step]:
            for model, dm in d.groupby("model"):
                a = dm[dm["arm"] == cand].set_index("fold")
                b = dm[dm["arm"] == base].set_index("fold")
                common = a.index.intersection(b.index)
                for met in mets:
                    r = dict(step=step, task=task, unit_or_study=layer, comparison=f"{cand} vs {base}", question=question, model=model, metric=met,
                             n_folds=len(common))
                    if len(common) < 2:
                        if len(a) or len(b):
                            Dl.append({**r, "verdict": "not evaluable"})
                        continue
                    diff = a.loc[common, met].to_numpy() - b.loc[common, met].to_numpy()
                    if met in K.LOWER_IS_BETTER:
                        diff = -diff
                    mean, lo, hi = K.corrected_t(diff, ratio, 1.0)
                    Dl.append({**r, "difference": mean, "ci_low": lo, "ci_high": hi, "verdict": K.verdict(lo, hi)})
    D = pd.DataFrame(Dl)
    D.to_csv(OUT / "paired_differences.csv", index=False)
    Cn = D.groupby(["step", "task", "unit_or_study", "comparison", "question", "metric", "verdict"]).size().unstack("verdict", fill_value=0).reset_index()
    for v in ("improves", "does not change", "gets worse", "not evaluable"):
        if v not in Cn:
            Cn[v] = 0
    Cn = Cn[["step", "task", "unit_or_study", "comparison", "question", "metric", "improves", "does not change", "gets worse", "not evaluable"]]
    Cn.to_csv(OUT / "counts.csv", index=False)
    A = Cn[Cn["metric"] == "AUROC"].copy()
    A["models_evaluated"] = A["improves"] + A["does not change"] + A["gets worse"]
    A.drop(columns=["metric"]).to_csv(OUT / "answers.csv", index=False)

    F = [pd.read_csv(f) for s in STEPS if (f := Q.HERE / s / "output" / "failures.csv").exists()]
    F = pd.concat(F, ignore_index=True) if F else pd.DataFrame(columns=["task", "arm", "model", "reason", "message", "job", "fold"])
    if len(F):
        F["message"] = F["message"].astype(str).str.replace(r"\s*\n\s*", " ", regex=True)
    F.to_csv(OUT / "failures.csv", index=False)
    if len(F):
        F["where"] = F["task"].str.split("/").str[1:].str.join(" ") + " " + F["arm"]
        N = (F.groupby(["model", "reason"]).agg(jobs=("job", "nunique"), where=("where", lambda w: "; ".join(sorted(set(w)))),
                                                example_message=("message", "first")).reset_index())
        N.to_csv(OUT / "models_not_run.csv", index=False)
    pd.set_option("display.width", 250)
    print(A.drop(columns=["question"]).to_string(index=False))
    print(F.groupby(["model", "reason"]).size().to_string() if len(F) else "no failures")


if __name__ == "__main__":
    main()
