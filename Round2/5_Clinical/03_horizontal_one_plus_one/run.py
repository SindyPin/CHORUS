"""Step 03: the clinical model on one unit plus one (Amendment 5, C), per task.

Anchor (Amendment 7, item 2): the largest unit that records age and sex (both kept) and has at least 10 per class.
Partner (Amendment 6): the unit sharing the most kept variables with the anchor, sharing at least age and sex; ties go
to the larger unit. If no unit qualifies, step 03 is skipped for that task and the reason is written. Units are groups
of people, so every person is counted once (the Su2022 and Cervia2024 layers are one unit each).

  A_all      anchor training fold, every variable kept in the anchor
  A_shared   anchor training fold, the variables kept in both
  B_naive    the same plus every partner participant
  B_scaled   the same after z scoring the numeric variables within each cohort (categorical variables as they are)

Tested on the anchor's held out folds; category codes are shared by anchor and partner. 15 folds (Amendment 3).

  python Round2/5_Clinical/03_horizontal_one_plus_one/run.py [--time]
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
STEP = "03_horizontal_one_plus_one"


def choose(U: pd.DataFrame, D: pd.DataFrame, task: str) -> tuple[dict | None, list]:
    """Anchor (Amendment 7, item 2): the largest unit that records age and sex and has 10 per class. Partner (Amendment 6): the unit
    sharing the most kept variables with the anchor, sharing at least age and sex; ties go to the larger unit. If no
    unit qualifies, the task is skipped and the reason returned."""
    ok = U[(U["task"] == task) & U["ten_per_class"] & (U["n_variables_kept"] > 0)].sort_values(["n", "unit"], ascending=[False, True])
    age_sex = ok[[{"age", "sex"} <= set(Q.kept_vars(D, u, task)) for u in ok["unit"]]]
    if age_sex.empty:
        return None, [dict(task=task, anchor="", reason="step 03 skipped: no unit records age and sex with 10 per class")]
    a = age_sex["unit"].iloc[0]
    ok = pd.concat([age_sex.iloc[:1], ok[ok["unit"] != a]])  # candidates: every other unit with clinical data and 10 per class
    ka = set(Q.kept_vars(D, a, task))
    cands = []
    for r in ok.iloc[1:].itertuples():
        shared = ka & set(Q.kept_vars(D, r.unit, task))
        cands.append(dict(task=task, anchor=a, candidate=r.unit, n=r.n, n_shared=len(shared), shares_age_and_sex={"age", "sex"} <= shared,
                          shared="; ".join(sorted(shared))))
    C = pd.DataFrame(cands)
    q = C[C["shares_age_and_sex"]] if len(C) else C
    if q.empty:
        why = (f"anchor {a} (n = {int(ok['n'].iloc[0])}; kept: {'; '.join(Q.kept_vars(D, a, task))}) "
               + ("does not record both age and sex" if not {"age", "sex"} <= ka else "has no unit sharing both age and sex"))
        return None, [dict(task=task, anchor=a, reason=f"step 03 skipped: {why}")] + cands
    b = q.sort_values(["n_shared", "n", "candidate"], ascending=[False, False, True]).iloc[0]
    shared = [v for v in Q.kept_vars(D, a, task) if v in set(b["shared"].split("; "))]
    overlap = "; ".join(v for k, v in Q.POSSIBLE_OVERLAP.items() if {a, b["candidate"]} <= k)
    return dict(task=task, anchor=a, partner=b["candidate"], shared="; ".join(shared), anchor_variables="; ".join(Q.kept_vars(D, a, task)),
                possible_overlap=overlap), cands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(Q.HERE / "output/units.csv").fillna("")
    D = pd.read_csv(Q.HERE / "variables.csv")
    runner = K.Runner(OUT, STEP, folder="5_Clinical")
    jobs, pairs, passed_all = [], [], []
    for task in Q.TASKS:
        pr, passed = choose(U, D, task)
        passed_all += passed
        if pr is None:
            continue
        RA = pickle.load(open(Q.CACHE / f"{task}__{pr['anchor']}.pkl", "rb"))["R"]
        RP = pickle.load(open(Q.CACHE / f"{task}__{pr['partner']}.pkl", "rb"))["R"]
        y, yp = RA["y"].to_numpy(), RP["y"].to_numpy()
        fl = K.folds(y, "classification")
        shared = pr["shared"].split("; ")
        pairs.append({**pr, "n_anchor_LC": int(y.sum()), "n_anchor_recovered": int((y == 0).sum()),
                      "n_partner_LC": int(yp.sum()), "n_partner_recovered": int((yp == 0).sum())})
        for arm in ("A_all", "A_shared", "B_naive", "B_scaled"):
            vs = pr["anchor_variables"].split("; ") if arm == "A_all" else shared
            both = pd.concat([RA, RP], ignore_index=True)
            N, C, _, _ = Q.arrays(both, vs)  # one coding for anchor and partner
            n = len(RA)
            extra = dict(N=N[:n], C=C[:n])
            if arm.startswith("B"):
                extra.update(Np=N[n:], Cp=C[n:], yp=yp)
            p = Q.CACHE / f"03__{task}__{arm}.npz"
            K.save_cache(p, N[:n], y, fl, **extra)
            for model in Q.models():
                jobs.append(dict(id=f"03__{task}__{pr['anchor']}__{arm}__{model}", control=task, layer=pr["anchor"], arm=arm, model=model,
                                 task="classification", k=0, cache=str(p), module="clinical"))
    pd.DataFrame(pairs).to_csv(OUT / "pairs.csv", index=False)
    pd.DataFrame(passed_all).to_csv(OUT / "partner_candidates.csv", index=False)
    print(pd.DataFrame(pairs).to_string(index=False))
    K.record_static_failures(runner, jobs, {})
    K.run_step(runner, jobs, OUT, STEP, a.time)


if __name__ == "__main__":
    main()
