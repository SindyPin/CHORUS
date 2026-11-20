"""Step 02: one cohort plus one (A_shared, B_naive, B_scaled), on the anchors and partners chosen in step 01.

  A_shared   anchor alone, features shared with the partner
  B_naive    anchor training fold plus all partner participants, shared features, no adjustment
  B_scaled   the same after z scoring each feature within each cohort, without labels

Always tested on the anchor's held out fold (the folds of step 01).

  python Round2/0_Controls/02_horizontal_one_plus_one/run.py --time
  python Round2/0_Controls/02_horizontal_one_plus_one/run.py
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import common as K  # noqa: E402

OUT = Path(__file__).resolve().parent / "output"
STEP01 = Path(__file__).resolve().parents[1] / "01_one_cohort" / "output"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--time", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    A = pd.read_csv(STEP01 / "control_anchors.csv").fillna("")
    runner = K.Runner(OUT, "02_horizontal_one_plus_one")
    jobs, info, used = [], [], []
    for r in A.itertuples():
        if not (r.anchor and r.partner):
            continue
        paths, emb_fail, inf, (_, mp) = K.build_horizontal(r.control, r.layer, r.anchor, r.partner)
        info.append({**inf, "same_scale": bool(r.same_scale),
                     "interpretable_arms": "A_shared, B_naive, B_scaled" if r.same_scale else "A_shared, B_scaled (B_naive reported only: Amendment 1, item 5)"})
        used.append(mp[["study_id", "subject_id", "sample_id"]].assign(control=r.control, layer=r.layer, role="partner"))
        j = K.jobs_for(r.control, r.layer, paths, ["A_shared", "B_naive", "B_scaled"], "02")
        K.record_static_failures(runner, j, emb_fail)
        jobs += j
    pd.DataFrame(info).to_csv(OUT / "pairs.csv", index=False)
    pd.concat(used).to_csv(K.CACHE / "02_partner_participants.csv", index=False)
    K.run_step(runner, jobs, OUT, "02", a.time)


if __name__ == "__main__":
    main()
