"""GSE225349: SomaScan (v4) serum proteomics of the Mount Sinai COVID-19 Biobank (Paranjpe et al., Commun Med 2023,
PMID 37308534), the same participants as GSE215865 (RNA-seq) and the Synapse clinical table. Built as an acute
proteomics cohort; no model is run on it here.

Input (GEO, downloaded 21 November 2026 into data/proteomics/GSE225349/raw/, not committed):
  GSE225349_non-normalized_counts_covid_only.txt.gz   aptamers x 1,105 samples, log2 RFU (the "covid only" set, every
                                                      time point); the "covid_AKI" file (437 samples, one per patient,
                                                      linear RFU) is a subset of these samples and is not used
  GSE225349_somalogic_protein_ids_validated.csv.gz    SeqId to UniProt
  GSE225349_series_matrix.txt.gz                      per sample: subject id, blood sample id, sex, age, days since first
                                                      sample, COVID-19 status
Sample id: the GEO blood sample id (Subj_xxxxxxxxT<n>), which is the Synapse Blood_Sample_ID; subject id
"GSE215865:Subj_xxxxxxxx", the same as in GSE215865, so that the two layers are one group of people.
Proteins: the Round1 somascan adapter (UniProt keys; aptamers of one protein averaged on the log2 scale) and the Round1
write_cohort (ClinicalSchema validation). lc_status (Amendment 7; the TACO labels are withdrawn): the post COVID
checklist labels rebuilt from the Synapse table (Round2/0_Labels/GSE215865_checklist), per subject: COVID positive in
this series and PASC is LC, no PASC is COVID_recovered, never answered is unknown; COVID negative is never_infected.
The symptoms only labels (sensitivity) are counted in output/overlap.csv; the metadata schema has no column for them.

Writes data/proteomics/GSE225349/GSE225349_{matrix,features}.csv and data/meta/GSE225349_meta.csv (not committed),
and output/build_summary.csv and output/overlap.csv (aggregate counts, committed).

  python Round2/0_Labels/GSE225349/build.py
"""

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
RAW = CODE / "data/proteomics/GSE225349/raw"
sys.path.insert(0, str(CODE / "Round1/scripts/proteomics"))
from common import write_cohort  # noqa: E402  (Round1 builder helper)

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "GSE225349"
PLATFORM = "SomaScan discovery platform v4 (GEO GSE225349, covid only set, log2 RFU)"


def series_matrix() -> pd.DataFrame:
    rows = {}
    with gzip.open(RAW / "GSE225349_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                rows["gsm"] = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = vals[0].split(":")[0].strip()
                rows[key] = [v.split(":", 1)[1].strip() if ":" in v else np.nan for v in vals]
    return pd.DataFrame(rows)


def main():
    (HERE / "output").mkdir(parents=True, exist_ok=True)
    A = pd.read_csv(RAW / "GSE225349_non-normalized_counts_covid_only.txt.gz", sep="\t", dtype={"gene": str, "ID REF": str})
    ann = pd.read_csv(RAW / "GSE225349_somalogic_protein_ids_validated.csv.gz", dtype=str)
    # aptamer id: the matrix writes SeqId 10000-28 as the number 10000.28, so a suffix ending in 0 loses it (10063-10
    # becomes 10063.1). Matched on the numeric form, which is unique in the annotation; the two files are also in the
    # same order (checked: every row agrees by position).
    num = ann["SeqId"].map(lambda s: round(float(s.replace("-", ".", 1)), 6))
    if not num.is_unique:
        raise ValueError("SeqId numeric forms are not unique")
    A["SeqId"] = A["gene"].astype(float).round(6).map(dict(zip(num, ann["SeqId"])))
    agree = int(A["SeqId"].notna().sum())
    by_position = int(sum(a == b for a, b in zip(A["SeqId"], ann["SeqId"])))
    sample_cols = [c for c in A.columns if c.startswith("Subj_")]
    wide = A.set_index("SeqId")[sample_cols].T  # samples x aptamers
    adapted = to_feature_space("somascan", wide, annotation=ann, aptamer_col="SeqId", uniprot_col="UniProt")

    G = series_matrix().rename(columns={"subject id": "subj", "blood sample id": "bsid"})
    G = G.set_index("bsid")
    CHK = pd.read_csv(CODE / "data/clinical/GSE215865/GSE215865_checklist_labels.csv").set_index("subject_id")
    meta = pd.DataFrame(index=adapted.X.index)
    meta.index.name = "sample_id"
    g = G.reindex(meta.index)
    meta["subject_id"] = "GSE215865:" + g["subj"]
    meta["lc_status"] = label_of(g["subj"], g["covid-19 status"], CHK["rebuilt"]).values
    meta["case_definition"] = "study_specific"
    meta["sample_matrix"] = "serum"
    meta["age"] = pd.to_numeric(g["age"], errors="coerce")
    meta["sex"] = g["Sex"].str.lower()
    meta["hospitalised"] = True
    meta["platform"] = PLATFORM
    meta["sample_phase"] = np.where(meta["lc_status"] == "never_infected", "not_applicable", "acute")
    meta["timepoint"] = pd.to_numeric(g["days since first sample"], errors="coerce")
    meta["collection_country"] = "United States"
    meta["overlaps_with_study"] = "GSE215865 (same participants and blood samples; Synapse syn35876735)"
    meta["notes"] = ("timepoint is days since the first biobank sample; acute hospital draws; lc_status from the post COVID checklist "
                     "(Synapse syn35876735, Thompson 2023 definition; Amendment 7); covid-19 status " + g["covid-19 status"].fillna("unknown"))
    write_cohort(STUDY, adapted, meta, root=CODE / "data/proteomics")
    S = [("aptamers in the matrix", len(A)), ("aptamers matched to the annotation (numeric SeqId)", agree), ("of which in the same position as the annotation", by_position),
         ("samples", len(meta)), ("subjects", meta["subject_id"].nunique()), ("UniProt features", adapted.X.shape[1]),
         ("transform", adapted.notes.get("transform")), ("aptamers without UniProt", adapted.notes.get("aptamers_without_uniprot"))]
    S += [(f"subjects with lc_status {k}", int(v)) for k, v in meta.drop_duplicates("subject_id")["lc_status"].value_counts().items()]
    pd.DataFrame(S, columns=["quantity", "value"]).to_csv(HERE / "output/build_summary.csv", index=False)
    print(pd.DataFrame(S, columns=["quantity", "value"]).to_string(index=False))
    meta["lc_status_symptoms_only"] = label_of(g["subj"], g["covid-19 status"], CHK["rebuilt_symptoms_only"]).values
    pd.set_option("display.width", 250)
    print(pd.crosstab(meta.drop_duplicates("subject_id")["lc_status"], meta.drop_duplicates("subject_id")["lc_status_symptoms_only"], margins=True))
    overlap(meta, "lc_status")


def label_of(subj: pd.Series, covid: pd.Series, rebuilt: pd.Series) -> pd.Series:
    r = subj.map(rebuilt).fillna("unknown")
    return pd.Series(np.select([covid != "Positive", r == "PASC", r == "no_PASC"], ["never_infected", "LC", "COVID_recovered"], "unknown"), index=subj.index)


def overlap(soma: pd.DataFrame, _=None):
    """Participants with RNA-seq, SomaScan and clinical data on the same first acute sample, by LC status, under the
    main (checklist) and the sensitivity (symptoms only) labels of Round2/0_Labels (Amendment 7)."""
    S = pd.read_csv(CODE / "Round2/0_Labels/output/sample_labels.csv")
    r = S[(S["study_id"] == "GSE215865") & (S["sample_phase_round2"] == "acute")].copy()
    r["bsid"] = r["sample_id"].str.replace(r"_Plate_\d+$", "", regex=True)
    T = pd.read_csv(CODE / "data/clinical/GSE215865/clinical_data_table.csv", low_memory=False, usecols=["Subject_ID", "Blood_Sample_ID", "Days_Since_First_Sample"])
    day = T.dropna(subset=["Blood_Sample_ID"]).set_index("Blood_Sample_ID")["Days_Since_First_Sample"]
    r["day"] = r["bsid"].map(day)
    first_rna = r.sort_values(["day", "bsid"]).drop_duplicates("subject_id")
    s = soma.reset_index()
    s = s[s["sample_phase"] == "acute"]
    s["day"] = s["sample_id"].map(day)
    first_soma = s.sort_values(["day", "sample_id"]).drop_duplicates("subject_id").set_index("subject_id")["sample_id"]
    soma_ids = set(s["sample_id"])
    clin_subj = set("GSE215865:" + T["Subject_ID"].astype(str))
    f = first_rna.set_index("subject_id")
    out = pd.DataFrame({"main": f["lc_status"], "symptoms_only": f["lc_status_symptoms_only"], "first_rna_sample": f["bsid"]})
    out["clinical"] = out.index.isin(clin_subj)
    out["first_rna_sample_has_somascan"] = out["first_rna_sample"].isin(soma_ids)
    out["first_rna_is_first_somascan"] = out["first_rna_sample"] == first_soma.reindex(out.index)
    out["any_somascan"] = out.index.isin(set(s["subject_id"]))
    rows = []
    for version in ("main", "symptoms_only"):
        for lb in ("LC", "COVID_recovered"):
            o = out[out[version] == lb]
            rows.append(dict(labels=version, label=lb, rna_first_acute_sample=len(o), with_clinical=int(o["clinical"].sum()),
                             with_somascan_any_sample=int(o["any_somascan"].sum()),
                             rna_and_somascan_on_the_same_first_sample=int(o["first_rna_sample_has_somascan"].sum()),
                             and_it_is_also_the_first_somascan_sample=int(o["first_rna_is_first_somascan"].sum()),
                             rna_somascan_clinical_same_first_sample=int((o["first_rna_sample_has_somascan"] & o["clinical"]).sum())))
        som = soma[soma["sample_phase"] == "acute"].drop_duplicates("subject_id")
        col = "lc_status" if version == "main" else "lc_status_symptoms_only"
        for lb in ("LC", "COVID_recovered"):
            rows.append(dict(labels=version, label=f"GSE225349 subjects {lb} (any acute SomaScan sample)", rna_first_acute_sample=np.nan,
                             with_somascan_any_sample=int((som[col] == lb).sum())))
    O = pd.DataFrame(rows)
    O.to_csv(HERE / "output/overlap.csv", index=False)
    pd.set_option("display.width", 250)
    print(O.T.to_string(header=False))


if __name__ == "__main__":
    main()
