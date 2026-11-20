"""Variable decisions of Round2/5_Clinical (Amendment 5, B6) and the per participant clinical tables.

Writes:
  variables.csv          every column recorded about the participants, per task and unit: status (eligible, early,
                         timing, excluded), kind, type, reason, participants, recorded, missing rate, kept
  who_mapping.csv        the severity mapping onto the WHO Clinical Progression Scale
  output/units.csv       per task and unit: studies, LC and recovered participants, variables kept, timing variable
  cache/<task>__<unit>.pkl   the per participant raw tables (participant level data: not committed)

  python Round2/5_Clinical/build_variables.py
"""

import pickle
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import clinical as Q  # noqa: E402


def synapse_rows(R: pd.DataFrame, chosen: pd.DataFrame, task: str, unit: str) -> list:
    """Every column of the GSE215865 Synapse table that is not already a variable row: its decision and the number of
    the unit's participants with a value on or before the day of their sample."""
    import re

    import numpy as np

    T = pd.read_csv(Q.CODE / Q.SYN.TABLE, low_memory=False)
    day = T.dropna(subset=["Blood_Sample_ID"]).set_index("Blood_Sample_ID")["Days_Since_First_Sample"]
    ch = chosen[chosen["study"] == "GSE215865"].copy()
    ch["subj"] = ch["subject_id"].str.replace("GSE215865:", "", regex=False)
    ch["cut"] = ch["sample_id"].map(lambda s: day.get(re.sub(r"_Plate_\d+$", "", s), np.nan))
    T = T.merge(ch[["subj", "cut"]], left_on="Subject_ID", right_on="subj")
    T = T[T["Days_Since_First_Sample"] <= T["cut"]]
    n = len(R)
    out = []
    for c, (status, kind, typ, why, used) in Q.SYNAPSE_DECISIONS.items():
        if status == "eligible" and used == c:
            continue  # a variable row of its own (decisions)
        rec = int(T.loc[T[c].notna(), "Subject_ID"].nunique())
        out.append(dict(task=task, unit=unit, variable=f"Synapse: {c}", status=status, kind=kind, type="",
                        reason=why + (f" (as {used})" if used and used != c and "used as" not in why else ""), n=n, n_recorded=rec,
                        missing_rate=round(1 - rec / n, 3), kept=False))
    return out


def units_for(task: str, S: pd.DataFrame) -> dict:
    out = {}
    for st in Q.studies_for(task, S):
        out.setdefault(Q.UNIT_OF.get(st, st), []).append(st)
    return out


def main():
    Q.CACHE.mkdir(exist_ok=True)
    (Q.HERE / "output").mkdir(exist_ok=True)
    S = Q.labels()
    D, U = [], []
    for task in Q.TASKS:
        for unit, studies in sorted(units_for(task, S).items()):
            R, chosen = Q.raw_table(studies, task, S)
            if R.empty:
                continue
            dec = Q.decisions(R, unit, task)
            if "GSE215865" in studies and Q.SYNAPSE_DECISIONS:
                dec += synapse_rows(R, chosen, task, unit)
            D += dec
            Dd = pd.DataFrame(dec)
            kept = list(Dd.loc[Dd["kept"], "variable"])
            tv = Dd[(Dd["status"] == "timing") & Dd["reason"].str.contains("used in the timing")]
            n1, n0 = int((R["y"] == 1).sum()), int((R["y"] == 0).sum())
            overlap = "; ".join(v for k, v in Q.POSSIBLE_OVERLAP.items() if unit in k)
            U.append(dict(task=task, unit=unit, studies=" + ".join(studies), n_LC=n1, n_recovered=n0, n=len(R),
                          ten_per_class=n1 >= 10 and n0 >= 10, n_variables_kept=len(kept), variables_kept="; ".join(kept),
                          timing_variable=tv["variable"].iloc[0] if len(tv) else "", possible_overlap=overlap))
            with open(Q.CACHE / f"{task}__{unit}.pkl", "wb") as f:
                pickle.dump(dict(R=R, chosen=chosen), f)
    D = pd.DataFrame(D)
    D.to_csv(Q.HERE / "variables.csv", index=False)
    W = [dict(source="acute_severity (Round1 harmonised)", value=k, who_cps=v[0], note=v[1]) for k, v in Q.WHO_CPS.items()]
    W += [dict(source="GSE224615 who_category_acute (WHO 8 point ordinal scale)", value=k, who_cps=v[0], note=v[1]) for k, v in Q.WHO8_CPS.items()]
    pd.DataFrame(W).to_csv(Q.HERE / "who_mapping.csv", index=False)
    U = pd.DataFrame(U)
    U.to_csv(Q.HERE / "output/units.csv", index=False)
    pd.set_option("display.width", 250)
    print(U.drop(columns=["variables_kept", "possible_overlap"]).to_string(index=False))


if __name__ == "__main__":
    main()
