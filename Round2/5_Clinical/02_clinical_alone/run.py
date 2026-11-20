"""Step 02: every model on the clinical variables alone (Amendment 5, C), per unit with at least 10 per class.

Arms (same participants and folds within a unit):
  clinical              every variable kept in the unit (variables.csv)
  age_sex               age and sex only (Lopez2023: its sex code)
  clinical_plus_timing  clinical plus time since infection: sensitivity analysis only (Amendment 5, B2)

Encoding inside the training fold, no feature selection (clinical.encode). 15 folds (Amendment 3).

  python Round2/5_Clinical/02_clinical_alone/run.py [--time]
"""

import argparse
import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import clinical as Q  # noqa: E402

K = Q.K
OUT = Path(__file__).resolve().parent / "output"
STEP = "02_clinical_alone"


def arm_sets(kept: list, timing: str) -> dict:
    a = {"clinical": kept, "age_sex": [v for v in ("age", "sex", "sex_code") if v in kept]}
    if timing:
        a["clinical_plus_timing"] = kept + [timing]
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(Q.HERE / "output/units.csv").fillna("")
    D = pd.read_csv(Q.HERE / "variables.csv")
    runner = K.Runner(OUT, STEP, folder="5_Clinical")
    jobs, used, skipped = [], [], []
    for u in U.itertuples():
        kept = Q.kept_vars(D, u.unit, u.task)
        if not u.ten_per_class or not kept:
            skipped.append(dict(task=u.task, unit=u.unit, arm="all", reason="fewer than 10 per class" if not u.ten_per_class else "no clinical variable kept"))
            continue
        with open(Q.CACHE / f"{u.task}__{u.unit}.pkl", "rb") as f:
            R = pickle.load(f)["R"]
        y = R["y"].to_numpy()
        fl = K.folds(y, "classification")
        seen = {}
        for arm, vs in arm_sets(kept, u.timing_variable).items():
            if not vs:
                skipped.append(dict(task=u.task, unit=u.unit, arm=arm, reason="neither age nor sex is kept in this unit"))
                continue
            if tuple(vs) in seen:
                skipped.append(dict(task=u.task, unit=u.unit, arm=arm, reason=f"same variables as the {seen[tuple(vs)]} arm"))
                continue
            seen[tuple(vs)] = arm
            N, C, nn, cn = Q.arrays(R, vs)
            p = Q.CACHE / f"02__{u.task}__{u.unit}__{arm}.npz"
            K.save_cache(p, N, y, fl, N=N, C=C)
            used.append(dict(task=u.task, unit=u.unit, arm=arm, sensitivity=arm == "clinical_plus_timing", n_LC=int(y.sum()),
                             n_recovered=int((y == 0).sum()), variables="; ".join(vs), numeric="; ".join(nn), categorical="; ".join(cn)))
            for model in Q.models():
                jobs.append(dict(id=f"02__{u.task}__{u.unit}__{arm}__{model}", control=u.task, layer=u.unit, arm=arm, model=model,
                                 task="classification", k=0, cache=str(p), module="clinical"))
    pd.DataFrame(used).to_csv(OUT / "arms.csv", index=False)
    pd.DataFrame(skipped).to_csv(OUT / "arms_not_run.csv", index=False)
    K.record_static_failures(runner, jobs, {})
    K.run_step(runner, jobs, OUT, STEP, a.time)


if __name__ == "__main__":
    main()
