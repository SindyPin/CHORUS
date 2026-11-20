"""Step 04: omics and clinical data on the same participants and folds (Amendment 5, C).

Every omics study of a unit with at least 10 per class and at least one clinical variable kept (Su2022 and
Cervia2024 give one comparison per omics layer). Participants: the unit's participants with a sample of that study in
the task's phase (the earliest); their clinical values come from the unit's table.

  omics        the omics layer alone: training fold median imputation, Welch t top 100 (vertical rule), standardised
  clinical     the clinical variables alone: encoding inside the training fold, no feature selection
  joined       both (early fusion): the 100 omics features plus every clinical column
  late_fusion  mean of the omics and clinical predictions of the same model, fold by fold (no extra fit)

The expression models (BulkFormer, BulkRNABert) are not run here: they embed RNA and have no clinical or joined form,
and every comparison pairs the same model across arms. 15 folds (Amendment 3).

  python Round2/5_Clinical/04_vertical/run.py [--time]
"""

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import clinical as Q  # noqa: E402

K = Q.K
OUT = Path(__file__).resolve().parent / "output"
STEP = "04_vertical"
ARMS = ["omics", "clinical", "joined"]


def late_fusion(runner, jobs):
    by = {(j["control"], j["layer"], j["model"], j["arm"]): j for j in jobs}
    for (task, study, model, arm), j in by.items():
        if arm != "omics":
            continue
        j2 = by[(task, study, model, "clinical")]
        jl = {**j, "id": j["id"].replace("__omics__", "__late_fusion__"), "arm": "late_fusion"}
        if runner.failed(j["id"]) or runner.failed(j2["id"]):
            if not runner.failed(jl["id"]):
                r = (dict(reason=K.REMOVED[model], message="Amendment 2: removed from the model grid, so there is no late fusion")
                     if model in K.REMOVED else dict(reason="error", message="a single arm failed, so there is no late fusion"))
                runner.record(jl, dict(ok=False, fold=-1, seconds=0.0, **r), 0, 0)
            continue
        f1, f2 = runner.out / "raw" / f"{j['id']}_pred.csv", runner.out / "raw" / f"{j2['id']}_pred.csv"
        if not (f1.exists() and f2.exists()):
            continue
        P1, P2 = pd.read_csv(f1), pd.read_csv(f2)
        d = np.load(j["cache"])
        have = runner.done(jl["id"])
        for fold in sorted(f for f in set(P1["fold"]) & set(P2["fold"]) - have if f < K.N_FOLDS):
            p = (P1.loc[P1["fold"] == fold, "pred"].to_numpy() + P2.loc[P2["fold"] == fold, "pred"].to_numpy()) / 2
            te = d[f"te{fold}"]
            runner.record(jl, dict(ok=True, fold=fold, seconds=0.0, pred=p, info={}, metrics=K.metrics(d["y"][te], p, "classification")),
                          len(d[f"tr{fold}"]), len(te))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(Q.HERE / "output/units.csv").fillna("")
    D = pd.read_csv(Q.HERE / "variables.csv")
    E = K.config_entries()
    runner = K.Runner(OUT, STEP, folder="5_Clinical")
    jobs, used, skipped = [], [], []
    for u in U.itertuples():
        kept = Q.kept_vars(D, u.unit, u.task)
        if not u.ten_per_class or not kept:
            continue
        P = pickle.load(open(Q.CACHE / f"{u.task}__{u.unit}.pkl", "rb"))
        R, chosen = P["R"].set_index("subject_id"), P["chosen"]
        for study in u.studies.split(" + "):
            if K.PARENT.get(study, study) not in E:
                skipped.append(dict(task=u.task, unit=u.unit, study=study, reason="no omics matrix in the Round1 configs"))
                continue
            X, m, is_counts, _ = K.load(study)
            ch = chosen[(chosen["study"] == study) & chosen["sample_id"].isin(X.index)]
            y = ch["y"].to_numpy()
            if min(int(y.sum()), int((y == 0).sum())) < K.MIN_PER_CLASS:
                skipped.append(dict(task=u.task, unit=u.unit, study=study, reason="fewer than 10 per class with a sample of this layer"))
                continue
            Xo = K.transform(X.loc[ch["sample_id"]], is_counts).to_numpy()
            Rs = R.loc[ch["subject_id"]].reset_index()
            N, C, nn, cn = Q.arrays(Rs, kept)
            fl = K.folds(y, "classification")
            p = Q.CACHE / f"04__{u.task}__{study}.npz"
            K.save_cache(p, Xo, y, fl, N=N, C=C)
            used.append(dict(task=u.task, unit=u.unit, study=study, n_LC=int(y.sum()), n_recovered=int((y == 0).sum()),
                             n_omics_features=Xo.shape[1], clinical_variables="; ".join(kept)))
            for arm in ARMS:
                for model in Q.models():
                    jobs.append(dict(id=f"04__{u.task}__{study}__{arm}__{model}", control=u.task, layer=study, arm=arm, model=model,
                                     task="classification", k=K.K_VERTICAL, cache=str(p), module="clinical"))
    pd.DataFrame(used).to_csv(OUT / "cohorts_used.csv", index=False)
    pd.DataFrame(skipped).to_csv(OUT / "cohorts_not_used.csv", index=False)
    print(pd.DataFrame(used).to_string(index=False))
    K.record_static_failures(runner, jobs, {})
    K.run_step(runner, jobs, OUT, STEP, a.time)
    late_fusion(runner, jobs)


if __name__ == "__main__":
    main()
