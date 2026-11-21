"""GSE215865: how the TACO labels (github.com/SindyPin/TACO, commit cac2afb, Data/411_DCEG_Paths_average.csv, column
Long_COVID) were built, and which rule on the Synapse table (syn35876735) reproduces them. No label is changed.

The rule in TACO (ML_Models/DCE_CV_Abation.ipynb, data preparation cell, step 3), on original_clinical_data_table.csv
(the same content as the Synapse clinical_data_table.csv, checked cell by cell):
  per ROW of the table, 1 if any daily Post_COVID19_Symptom_* (not _Ever) is True, or SARSCoV2_Recovery_Status is
  "Not_Recovered", or Post_COVID19_Health_Ever is "Worse"; else 0. Then drop_duplicates("Subject_ID") keeps the FIRST
  row of each subject in file order.

Variants compared with TACO per subject (agree / disagree / TACO labels a subject the variant cannot label):
  a      Post_COVID19_Symptom_*_Ever only, any row
  b      a plus Health, Quality_Of_Life or Daily_Mood _Ever "Worse" (Thompson 2023)
  c      the daily Post_COVID19_* columns (symptom True, or health, quality of life or mood "Worse") over all rows,
         missing as false (every subject is labelled)
  d_a, d_b, d_c   each of a, b, c plus COVID19_Symptoms_At_Presentation_Reviewed == "COVID19_Symptoms" (the only
         other symptom column with a value per subject; the others are onset dates and days)
  e_taco          the TACO rule as written (first row of each subject)
  e_any_row       the TACO rule over every row of a subject (any row)
  e_first_row_*   the TACO rule on the first row without one of its three parts (to see which part decides)
  recovery_first  SARSCoV2_Recovery_Status on the first row alone ("Not_Recovered" 1, otherwise 0)

Writes per subject data/clinical/GSE215865/GSE215865_taco_rule_variants.csv (participant level, not committed) and
output/taco_rule_comparison.csv and output/taco_never_answered.csv (aggregate).

  python Round2/0_Labels/GSE215865_checklist/taco_rule.py <path to 411_DCEG_Paths_average.csv at cac2afb, first 2 columns>
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TABLE = CODE / "data/clinical/GSE215865/clinical_data_table.csv"
PRIVATE = CODE / "data/clinical/GSE215865/GSE215865_taco_rule_variants.csv"
DET = ["Health", "Quality_Of_Life", "Daily_Mood"]


def boolean(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().map({"TRUE": True, "FALSE": False})


def checklist(T: pd.DataFrame, sym: list, det: list, missing_false: bool) -> pd.Series:
    """Per subject: 1 if any symptom item True or any deterioration item "Worse" on any row; 0 if answered and none;
    NaN if never answered (unless missing_false)."""
    pos = T[sym].apply(boolean).fillna(False).any(axis=1) | T[det].eq("Worse").any(axis=1)
    ans = T[sym + det].notna().any(axis=1)
    g = pd.DataFrame({"s": T["Subject_ID"], "pos": pos, "ans": ans}).groupby("s")
    out = g["pos"].any().astype(float)
    if not missing_false:
        out[~g["ans"].any()] = np.nan
    return out


def taco_row(T: pd.DataFrame, daily_sym: list, sym=True, rec=True, health=True) -> pd.Series:
    """The TACO create_label, row by row (row[col] is True on the values as read by pandas)."""
    r = pd.Series(False, index=T.index)
    if sym:
        r |= T[daily_sym].apply(lambda s: s.map(lambda v: v is True)).any(axis=1)
    if rec:
        r |= T["SARSCoV2_Recovery_Status"].eq("Not_Recovered")
    if health:
        r |= T["Post_COVID19_Health_Ever"].eq("Worse")
    return r.astype(int)


def main(taco_path: str):
    (HERE / "output").mkdir(parents=True, exist_ok=True)
    T = pd.read_csv(TABLE, low_memory=False)
    taco = pd.read_csv(taco_path, usecols=["Subject_ID", "Long_COVID"]).set_index("Subject_ID")["Long_COVID"]
    ever_sym = [c for c in T.columns if c.startswith("Post_COVID19_Symptom_") and c.endswith("_Ever")]
    daily_sym = [c[:-5] for c in ever_sym]
    ever_det = [f"Post_COVID19_{d}_Ever" for d in DET]
    daily_det = [f"Post_COVID19_{d}" for d in DET]
    first = T.drop_duplicates("Subject_ID").index
    pres = (T["COVID19_Symptoms_At_Presentation_Reviewed"].eq("COVID19_Symptoms")).groupby(T["Subject_ID"]).any().astype(float)

    V = pd.DataFrame(index=sorted(T["Subject_ID"].unique()))
    V["a"] = checklist(T, ever_sym, [], False)
    V["b"] = checklist(T, ever_sym, ever_det, False)
    V["c"] = checklist(T, daily_sym, daily_det, True)
    for k in ("a", "b", "c"):
        # presentation symptoms add a 1; a subject the checklist cannot label gets 1 if symptomatic, else stays unlabelled
        V[f"d_{k}"] = np.where(pres.reindex(V.index) == 1, 1.0, V[k])
    rows = {"e_taco": dict(), "e_first_row_no_recovery": dict(rec=False), "e_first_row_no_symptoms": dict(sym=False),
            "e_first_row_no_health": dict(health=False)}
    for k, kw in rows.items():
        r = taco_row(T, daily_sym, **kw)
        V[k] = r.loc[first].set_axis(T.loc[first, "Subject_ID"]).reindex(V.index).astype(float)
    V["e_any_row"] = taco_row(T, daily_sym).groupby(T["Subject_ID"]).max().reindex(V.index).astype(float)
    V["recovery_first"] = T.loc[first].set_index("Subject_ID")["SARSCoV2_Recovery_Status"].eq("Not_Recovered").astype(float).reindex(V.index)
    V["answered_checklist"] = ~checklist(T, ever_sym + daily_sym, ever_det + daily_det, False).isna()
    V["recovery_status_first_row"] = T.loc[first].set_index("Subject_ID")["SARSCoV2_Recovery_Status"].reindex(V.index)
    V["taco"] = taco.reindex(V.index)
    V["taco_subject"] = V.index.isin(taco.index)
    V.index.name = "Subject_ID"
    V.to_csv(PRIVATE)

    S = V[V["taco_subject"]]
    out = []
    for k in ("a", "b", "c", "d_a", "d_b", "d_c", "e_taco", "e_any_row", "e_first_row_no_recovery", "e_first_row_no_symptoms",
              "e_first_row_no_health", "recovery_first"):
        x, t = S[k], S["taco"]
        lab = x.notna()
        out.append(dict(variant=k, taco_subjects=len(S), agree=int((x[lab] == t[lab]).sum()), disagree=int((x[lab] != t[lab]).sum()),
                        taco_1_variant_0=int(((t == 1) & (x == 0)).sum()), taco_0_variant_1=int(((t == 0) & (x == 1)).sum()),
                        taco_labels_variant_cannot=int((~lab).sum()), variant_1=int((x == 1).sum()), variant_0=int((x == 0).sum()),
                        all_767_subjects_1=int((V[k] == 1).sum()), all_767_subjects_0=int((V[k] == 0).sum())))
    O = pd.DataFrame(out)
    O.to_csv(HERE / "output/taco_rule_comparison.csv", index=False)
    N = S[~S["answered_checklist"]]
    W = N.groupby(["taco", "recovery_status_first_row"]).size().rename("subjects").reset_index()
    W.to_csv(HERE / "output/taco_never_answered.csv", index=False)
    pd.set_option("display.width", 250)
    print(O.to_string(index=False))
    print(f"\nTACO subjects who never answered the checklist: {len(N)}")
    print(W.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1])
