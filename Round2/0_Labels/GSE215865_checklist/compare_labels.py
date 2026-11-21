"""GSE215865: LC and recovered labels rebuilt from the post COVID checklist of the Mount Sinai COVID-19 Biobank
clinical table (Synapse syn35876735; Thompson et al., Nat Med 2023, 29:236-246), compared with the TACO labels used so
far. Under Amendment 7 (21 November 2026) these rebuilt labels replace the TACO labels in Round2: the checklist
definition below is the main label, and symptoms only (the *_Ever symptom items, without health, quality of life or
mood) is the sensitivity label; subjects who never answered are unknown under both (build_labels.py reads the
round2_main and round2_symptoms_only columns written here).

Definition (Thompson 2023, Methods "Clinical data generation" and Table 1):
  - checklist items are of two types, symptoms and clinical deterioration (health, quality of life and daily mood
    compared with before COVID-19); every item is coded true for a symptom present or for "Worse";
  - answers are collapsed over every checklist a participant completed: an item is true if it was true at any time;
  - PASC: at least one item true ("any PASC symptom"); no PASC: every answered item false ("no PASC symptoms");
  - participants with no checklist answer are unknown.
  The table's *_Ever columns are these collapsed answers. The paper analyses COVID-19 patients (controls answered
  about their hospitalisation instead), so LC and recovered are kept for SARS-CoV-2 positive participants only; both
  positivity sources are reported (GEO "covid-19 positive", used so far, and the table's COVID19_Ever_Positive).
  Variant: the paper removed items reported by 5 or fewer participants with RNA-seq before modelling; the labels are
  also rebuilt without those items, to see whether any label changes.

Writes (participant level, under data/: never committed):
  data/clinical/GSE215865/GSE215865_checklist_labels.csv   per subject: rebuilt label, TACO label, Round2 label, category
Writes (aggregate, committed):
  output/comparison_counts.csv, output/summary.csv

  python Round2/0_Labels/GSE215865_checklist/compare_labels.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TABLE = CODE / "data/clinical/GSE215865/clinical_data_table.csv"
TACO = CODE / "data/meta/GSE215865_labels_from_TACO.csv"
ROUND2 = CODE / "Round2/0_Labels/output/sample_labels.csv"
PRIVATE = CODE / "data/clinical/GSE215865/GSE215865_checklist_labels.csv"
DETERIORATION = ["Post_COVID19_Health_Ever", "Post_COVID19_Quality_Of_Life_Ever", "Post_COVID19_Daily_Mood_Ever"]


def boolean(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().map({"TRUE": True, "FALSE": False})


def rebuild(inv: pd.DataFrame, items: list) -> pd.Series:
    sym = [c for c in items if c not in DETERIORATION]
    det = [c for c in items if c in DETERIORATION]
    S = inv[sym].apply(boolean)
    answered = inv[sym + det].notna().any(axis=1)
    any_true = S.fillna(False).any(axis=1) | inv[det].eq("Worse").any(axis=1)
    return pd.Series(np.where(~answered, "unknown", np.where(any_true, "PASC", "no_PASC")), index=inv.index)


def category(r) -> str:
    """Agreement category of one GEO subject (LC side: TACO 1; recovered side: TACO 0 and COVID positive)."""
    t, c, pos = r["taco_Long_COVID"], r["rebuilt"], r["geo_covid_positive"]
    if pd.isna(t):
        return "no TACO label"
    if c == "unknown":
        return f"TACO {int(t)}, no checklist answer"
    if t == 1 and c == "PASC":
        return "agree: LC"
    if t == 0 and c == "no_PASC":
        return "agree: no LC" + ("" if pos == "Positive" else " (COVID negative in GEO)")
    return f"disagree: TACO {int(t)}, checklist {c}" + ("" if pos == "Positive" else " (COVID negative in GEO)")


def main():
    (HERE / "output").mkdir(parents=True, exist_ok=True)
    T = pd.read_csv(TABLE, low_memory=False)
    inv = T.drop_duplicates("Subject_ID").set_index("Subject_ID")
    items = [c for c in T.columns if c.startswith("Post_COVID19_Symptom_") and c.endswith("_Ever")] + DETERIORATION
    lab = rebuild(inv, items)

    G = pd.read_csv(TACO).drop_duplicates("ch_subject_id").set_index("ch_subject_id")
    rna = G.index
    # variant: items true for 5 or fewer participants with RNA-seq removed (as before the paper's models)
    n_true = {c: int((boolean(inv.loc[inv.index.isin(rna), c]) if c not in DETERIORATION else inv.loc[inv.index.isin(rna), c].eq("Worse")).sum()) for c in items}
    rare = [c for c, n in n_true.items() if n <= 5]
    lab_v = rebuild(inv, [c for c in items if c not in rare])
    lab_s = rebuild(inv, [c for c in items if c not in DETERIORATION])  # Amendment 7: symptoms only (sensitivity)

    S = pd.read_csv(ROUND2)
    S = S[S["study_id"] == "GSE215865"].copy()
    S["subject_id"] = S["subject_id"].astype(str).str.replace("GSE215865:", "", regex=False)  # Round2 ids carry the study
    r2 = S.groupby("subject_id")["lc_status"].agg(lambda s: "; ".join(sorted(set(map(str, s)))))
    r2_acute = S[S["sample_phase_round2"] == "acute"].groupby("subject_id")["lc_status"].first()

    P = pd.DataFrame({"rebuilt": lab, "rebuilt_without_rare_items": lab_v, "rebuilt_symptoms_only": lab_s,
                      "synapse_COVID19_Ever_Positive": inv["COVID19_Ever_Positive"].astype(str),
                      "synapse_patient_classification": inv["Patient_Classification_At_First_Sample"]})
    P = P.join(G[["Long_COVID", "ch_covid-19_positive"]].rename(columns={"Long_COVID": "taco_Long_COVID", "ch_covid-19_positive": "geo_covid_positive"}), how="left")
    P["has_rnaseq"] = P.index.isin(rna)
    P["round2_label"] = r2.reindex(P.index)
    P["round2_label_acute"] = r2_acute.reindex(P.index)
    P["category"] = [category(r) if r["has_rnaseq"] else "no RNA-seq (not in GEO)" for _, r in P.iterrows()]
    P["rebuilt_round2_form"] = np.where((P["rebuilt"] == "PASC") & (P["geo_covid_positive"] == "Positive"), "LC",
                                        np.where((P["rebuilt"] == "no_PASC") & (P["geo_covid_positive"] == "Positive"), "COVID_recovered",
                                                 np.where(P["rebuilt"] == "unknown", "unknown", "not used (COVID negative or not in GEO)")))
    # Amendment 7: Round2 labels. COVID positive in GEO: PASC is LC, no PASC is COVID_recovered, no answer is unknown;
    # COVID negative in GEO: never_infected (as before; these subjects answered about their hospitalisation, not about
    # COVID-19). Not in GEO: no RNA-seq, not used.
    for col, src in (("round2_main", "rebuilt"), ("round2_symptoms_only", "rebuilt_symptoms_only")):
        pos = P["geo_covid_positive"] == "Positive"
        P[col] = np.select([~P["has_rnaseq"], ~pos, P[src] == "PASC", P[src] == "no_PASC"],
                           ["not in GEO", "never_infected", "LC", "COVID_recovered"], "unknown")
    P.index.name = "subject_id"
    P.to_csv(PRIVATE)

    C = P[P["has_rnaseq"]].groupby("category").size().rename("subjects").reset_index()
    C.to_csv(HERE / "output/comparison_counts.csv", index=False)
    g = P[P["has_rnaseq"]]
    acute = g[g["round2_label_acute"].isin(["LC", "COVID_recovered"])]
    rows = [
        ("Synapse subjects", len(P)), ("answered the checklist", int((P["rebuilt"] != "unknown").sum())),
        ("rebuilt PASC (all subjects)", int((P["rebuilt"] == "PASC").sum())), ("rebuilt no PASC (all subjects)", int((P["rebuilt"] == "no_PASC").sum())),
        ("paper Table 1: answered / PASC / no PASC", "232 / 195 / 37 (RNA-seq 165 / 140 / 25)"),
        ("GEO (RNA-seq) subjects", len(g)), ("GEO subjects with a checklist answer", int((g["rebuilt"] != "unknown").sum())),
        ("GEO subjects: rebuilt PASC / no PASC", f"{int((g['rebuilt'] == 'PASC').sum())} / {int((g['rebuilt'] == 'no_PASC').sum())}"),
        ("GEO subjects, COVID positive in GEO: rebuilt PASC / no PASC",
         f"{int(((g['rebuilt'] == 'PASC') & (g['geo_covid_positive'] == 'Positive')).sum())} / {int(((g['rebuilt'] == 'no_PASC') & (g['geo_covid_positive'] == 'Positive')).sum())}"),
        ("items removed in the variant (true for <= 5 RNA-seq subjects)", "; ".join(rare) or "none"),
        ("labels changed by the variant (all subjects)", int((P["rebuilt"] != P["rebuilt_without_rare_items"]).sum())),
        ("Round2 early detection participants (acute sample, LC or recovered)", len(acute)),
        ("  of whom with a checklist answer", int((acute["rebuilt"] != "unknown").sum())),
        ("  Round2 LC and rebuilt PASC", int(((acute["round2_label_acute"] == "LC") & (acute["rebuilt"] == "PASC")).sum())),
        ("  Round2 LC and rebuilt no PASC", int(((acute["round2_label_acute"] == "LC") & (acute["rebuilt"] == "no_PASC")).sum())),
        ("  Round2 recovered and rebuilt PASC", int(((acute["round2_label_acute"] == "COVID_recovered") & (acute["rebuilt"] == "PASC")).sum())),
        ("  Round2 recovered and rebuilt no PASC", int(((acute["round2_label_acute"] == "COVID_recovered") & (acute["rebuilt"] == "no_PASC")).sum())),
        ("  Round2 LC, no checklist answer", int(((acute["round2_label_acute"] == "LC") & (acute["rebuilt"] == "unknown")).sum())),
        ("  Round2 recovered, no checklist answer", int(((acute["round2_label_acute"] == "COVID_recovered") & (acute["rebuilt"] == "unknown")).sum())),
        ("positivity: GEO positive vs Synapse COVID19_Ever_Positive (GEO subjects)",
         "; ".join(f"{a} / {b}: {n}" for (a, b), n in g.groupby(["geo_covid_positive", "synapse_COVID19_Ever_Positive"]).size().items())),
    ]
    pd.DataFrame(rows, columns=["quantity", "value"]).to_csv(HERE / "output/summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(pd.DataFrame(rows, columns=["quantity", "value"]).to_string(index=False))
    print()
    print(C.to_string(index=False))


if __name__ == "__main__":
    main()
