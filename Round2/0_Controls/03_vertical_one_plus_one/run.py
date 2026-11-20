"""Step 03: one layer plus one, on the same participants (Su2022: Olink and plasma metabolomics; Cervia2024:
SomaScan and MS proteomics on the same sera).

Participants: those with the control variable recorded and a sample in both layers (the same sample id), one sample
per participant (the earliest). Folds on these participants, shared by every arm:

  layer_1      layer 1 alone (100 features)
  layer_2      layer 2 alone (100 features)
  joined       both joined, 100 features per layer selected separately (early fusion)
  late_fusion  mean of the layer_1 and layer_2 predictions of the same model (no extra fit)

  python Round2/0_Controls/03_vertical_one_plus_one/run.py --time
  python Round2/0_Controls/03_vertical_one_plus_one/run.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as K  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
PAIRS = {"Su2022": ("Su2022", "Su2022_metabolomics"), "Cervia2024": ("Cervia2024_SomaScan", "Cervia2024_MS")}


def build(cohort: str, control: str):
    l1, l2 = PAIRS[cohort]
    X1, _, m1 = K.control_data(l1, control)  # one per participant within each layer first
    _, m1all, c1, _ = K.load(l1)
    _, m2all, c2, _ = K.load(l2)
    common_ids = sorted(set(m1all["sample_id"].astype(str)) & set(m2all["sample_id"].astype(str)))
    X1, y, m = K.control_data(l1, control, samples=common_ids)
    X2, y2, m2 = K.control_data(l2, control, samples=list(m["sample_id"].astype(str)))
    ids = [s for s in m["sample_id"].astype(str) if s in set(m2["sample_id"].astype(str))]
    X1, X2 = X1.loc[ids], X2.loc[ids]
    y = K.outcome(m.set_index(m["sample_id"].astype(str)).loc[ids], control).to_numpy()
    fl = K.folds(y, K.TASK[control])
    p = K.CACHE / f"03_{control}_{cohort}.npz"
    K.save_cache(p, X1.to_numpy(), y, fl, X2=X2.to_numpy())
    info = dict(control=control, cohort=cohort, layer_1=l1, layer_2=l2, n_participants=len(ids),
                n_features_layer_1=X1.shape[1], n_features_layer_2=X2.shape[1])
    if control == "sex":
        info.update(n_female=int((y == 1).sum()), n_male=int((y == 0).sum()))
    else:
        info.update(age_min=float(y.min()), age_median=float(np.median(y)), age_max=float(y.max()))
    return p, info, m.set_index(m["sample_id"].astype(str)).loc[ids]


def late_fusion(runner, jobs):
    """Mean of the layer_1 and layer_2 held out predictions of the same model, fold by fold."""
    by = {(j["control"], j["cohort"], j["model"], j["arm"]): j for j in jobs}
    for (control, cohort, model, arm), j in by.items():
        if arm != "layer_1":
            continue
        j2 = by[(control, cohort, model, "layer_2")]
        f1, f2 = runner.out / "raw" / f"{j['id']}_pred.csv", runner.out / "raw" / f"{j2['id']}_pred.csv"
        jl = {**j, "id": j["id"].replace("__layer_1__", "__late_fusion__"), "arm": "late_fusion"}
        if runner.failed(j["id"]) or runner.failed(j2["id"]):
            if not runner.failed(jl["id"]):
                if model in K.REMOVED:
                    r = dict(reason=K.REMOVED[model], message="Amendment 2: removed from the model grid, so there is no late fusion")
                else:
                    r = dict(reason="error", message="a single layer model failed, so there is no late fusion")
                runner.record(jl, dict(ok=False, fold=-1, seconds=0.0, **r), 0, 0)
            continue
        if not (f1.exists() and f2.exists()):
            continue
        P1, P2 = pd.read_csv(f1), pd.read_csv(f2)
        d = np.load(j["cache"])
        have = runner.done(jl["id"])
        for fold in sorted(f for f in set(P1["fold"]) & set(P2["fold"]) - have if f < K.N_FOLDS):
            p = (P1.loc[P1["fold"] == fold, "pred"].to_numpy() + P2.loc[P2["fold"] == fold, "pred"].to_numpy()) / 2
            te = d[f"te{fold}"]
            runner.record(jl, dict(ok=True, fold=fold, seconds=0.0, pred=p, info={}, metrics=K.metrics(d["y"][te], p, j["task"])),
                          len(d[f"tr{fold}"]), len(te))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    runner = K.Runner(OUT, "03_vertical_one_plus_one")
    jobs, info, used = [], [], []
    for control in ("sex", "age"):
        for cohort in PAIRS:
            p, inf, m = build(cohort, control)
            info.append(inf)
            used.append(m[["study_id", "subject_id", "sample_id"]].assign(control=control, cohort=cohort))
            for arm in K.ARMS_V:
                for model in K.CLASSICAL + K.ENSEMBLE + K.FOUNDATION:
                    jobs.append(dict(id=f"03__{control}__{cohort}__{arm}__{model}", control=control, layer=cohort, cohort=cohort, arm=arm,
                                     model=model, task=K.TASK[control], k=K.K_VERTICAL, vertical=True, cache=str(p)))
    pd.DataFrame(info).to_csv(OUT / "cohorts_used.csv", index=False)
    pd.concat(used).to_csv(K.CACHE / "03_participants.csv", index=False)
    K.record_static_failures(runner, jobs, {})
    K.run_step(runner, jobs, OUT, "03", a.time)
    late_fusion(runner, jobs)


if __name__ == "__main__":
    main()
