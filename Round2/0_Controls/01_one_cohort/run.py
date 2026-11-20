"""Step 01: one cohort (A_all). Known answer controls: sex (classification) and age (regression).

Chooses the anchor and partner of each control and layer among the cohorts that record the variable (PROTOCOL.md
rule with Amendment 1), then runs A_all: the anchor alone with all its features, repeated stratified 5 fold cross
validation (10 repeats), every model.

  python Round2/0_Controls/01_one_cohort/run.py --time     fold 0 of every model only (timing pass: estimate only)
  python Round2/0_Controls/01_one_cohort/run.py            fold 0, the estimate, then all 50 folds if the step fits in 12 h (resumes)
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as K  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    C = pd.concat([K.candidates("sex"), K.candidates("age")], ignore_index=True)
    C.to_csv(OUT / "candidates.csv", index=False)
    A = K.choose(C)
    A.to_csv(OUT / "control_anchors.csv", index=False)
    print(A.to_string(index=False))
    runner = K.Runner(OUT, "01_one_cohort")
    jobs, info, used = [], [], []
    for r in A.itertuples():
        if not r.anchor:
            continue
        paths, emb_fail, inf, (ma, _) = K.build_horizontal(r.control, r.layer, r.anchor, r.partner)
        info.append(inf)
        used.append(ma[["study_id", "subject_id", "sample_id"]].assign(control=r.control, layer=r.layer, role="anchor"))
        j = K.jobs_for(r.control, r.layer, paths, ["A_all"], "01")
        K.record_static_failures(runner, j, emb_fail)
        jobs += j
    pd.DataFrame(info).to_csv(OUT / "cohorts_used.csv", index=False)
    pd.concat(used).to_csv(K.CACHE / "01_participants.csv", index=False)  # sample ids: not committed
    K.run_step(runner, jobs, OUT, "01", a.time)


if __name__ == "__main__":
    main()
