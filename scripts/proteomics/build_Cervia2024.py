"""Build the Cervia2024 cohorts (Cervia-Hasler C, ..., Boyman O 2024 Science 383:eadg7942,
"Persistent complement dysregulation with signs of thromboinflammation in active Long Covid",
doi:10.1126/science.adg7942, PMID 38236961; Mendeley Data 10.17632/dvf6yvrg4x.1, CC BY 4.0).

Two layers measured on the same people (Zurich cohort, serum), written as two cohorts:

  Cervia2024_SomaScan  (somascan)       39 healthy controls (first visit) and 113 COVID-19 patients
                                        at acute infection (Sampling_month ONE) and 6-month
                                        follow-up (SIX); 265 samples
  Cervia2024_MS        (ms_proteomics)  DIA LC-MS/MS of 6-month follow-up sera (119 columns)

Raw files (data/proteomics/Cervia2024/raw):
  ..._Proteomics_Clinical_Data_6M_timepoint.xlsx           sheet Data: clinical + log10 RFU (6 months)
  ..._Proteomics_Clinical_Data_Acute_plus_healthy-*.xlsx  sheet Data: clinical + log10 RFU (acute, healthy)
  ..._Proteomics_Clinical_Data_Labels-*.xlsx              per subject PACS_6M_woDys, PACS_12M
  ..._Proteomics_SomaScan_Aptamer_Info-*.xlsx             AptName, TargetFullName, Target, Entrez, UniProt
  ..._Mass_Spectrometry_6M_timepoint-*.xlsx               Protein.Group x samples (DIA-NN maxLFQ)

Reshaping before the adapters (the adapters are unchanged):
  * SomaScan values are deposited as log10 RFU (Info sheet; paper Methods). They are returned to
    RFU (10**x) so that the somascan adapter applies its own log2 transform; no other change.
  * Only human protein aptamers are kept (``keep``): the annotation has no organism or type column,
    so "No protein" controls, HCE hybridisation controls, Fc_MOUSE, and targets named with a
    non-human UniProt mnemonic suffix (_MOUSE, _AEQVI, _HV2BE, _ECOLI, ...) or HIV are dropped.
  * The MS table is already one row per protein group and one column per sample; it is passed as is.

Sample ids are the authors' SampleId with "COV-ZH " written "COV_ZH_" (the MS column style), so the
same serum sample has the same sample_id in both layers. subject_id is the authors' SubjectID.
No normalisation, imputation, or batch correction here (values are the authors' processed values).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, swissprot_accessions, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Cervia2024"
PREFIX = "Cervia_Hasler_et_al_2024_"
F_SIX = "Proteomics_Clinical_Data_6M_timepoint.xlsx"
F_ACUTE = "Proteomics_Clinical_Data_Acute_plus_healthy-81EAkK.xlsx"
F_LABELS = "Proteomics_Clinical_Data_Labels-TAklSO.xlsx"
F_APT = "Proteomics_SomaScan_Aptamer_Info-KnqtVJ.xlsx"
F_MS = "Mass_Spectrometry_6M_timepoint-VOZ1Hu.xlsx"

SOMA_PLATFORM = "SomaScan 7k (v4; 7596 SOMAmers deposited)"
MS_PLATFORM = "LC-MS/MS DIA (Orbitrap Eclipse; DIA-NN library-free, maxLFQ)"
COMORBIDITIES = ["Hypertonia", "Diabetes", "Heart", "Lung", "Malignancy", "Immuno", "Kidney", "Cerebro", "Asthma", "Allergic_disease", "Autoimmune_diseases"]
NONHUMAN_SUFFIX = re.compile(r"_(?:MOUSE|RAT|BOVIN|AEQVI|HV2BE|HV1|THETH|GEOSE|PHOPY|NOH4|RHIML|HELSU|SAMNI|VESMG|ECOLI|XENLA|CYACA)\b")


def sid(sample_id: str) -> str:
    """'COV-ZH 1.2' -> 'COV_ZH_1.2' (the MS column naming)."""
    return str(sample_id).strip().replace("COV-ZH ", "COV_ZH_")


def human_aptamers(apt: pd.DataFrame) -> pd.Series:
    """Boolean keep mask indexed by AptName: human protein aptamers only."""
    name = apt["TargetFullName"].fillna("").astype(str)
    target = apt["Target"].fillna("").astype(str)
    no_protein = name.str.lower().eq("no protein")
    hyb = target.str.startswith("HCE")
    fc_mouse = target.eq("Fc_MOUSE")
    nonhuman = name.str.contains(NONHUMAN_SUFFIX) | target.str.contains(NONHUMAN_SUFFIX) | name.str.contains(r"\bHIV\b|Xenopus") | target.str.contains(r"\bHIV\b")
    keep = ~(no_protein | hyb | fc_mouse | nonhuman) & apt["UniProt"].notna()
    return pd.Series(keep.values, index=apt["AptName"].values)


def severity_band(r) -> str:
    if r.get("ICU") == "YES":
        return "icu_or_ventilated"
    if r.get("Oxygen") == "YES":
        return "hospitalised_oxygen"
    if r.get("Patient_Care_followup") == "Hospitalized":
        return "hospitalised_no_oxygen"
    if r.get("Patient_Care_followup") in ("Outpatient", "Mild"):
        return "ambulatory"
    return "unknown"


def load(raw: Path):
    six = pd.read_excel(raw / (PREFIX + F_SIX), sheet_name="Data")
    ac = pd.read_excel(raw / (PREFIX + F_ACUTE), sheet_name="Data")
    lab = pd.read_excel(raw / (PREFIX + F_LABELS), sheet_name=0).set_index("SubjectID")
    apt = pd.read_excel(raw / (PREFIX + F_APT))
    ms = pd.read_excel(raw / (PREFIX + F_MS))
    return six, ac, lab, apt, ms


def subject_table(six: pd.DataFrame, ac: pd.DataFrame, lab: pd.DataFrame) -> pd.DataFrame:
    """One row per authors' SubjectID with the subject level clinical fields."""
    a = ac.set_index("SubjectID")
    subj = a[["COVID", "Sex", "BMI", "Oxygen", "ICU", "Post_Vaccine"] + COMORBIDITIES].copy()
    subj["comorbidity_count"] = (a[COMORBIDITIES] == "YES").sum(axis=1)
    fu = six.set_index("SubjectID")[["Patient_Care", "COVID19_Severity", "COVID19_Severity_Grade"]]
    fu.columns = ["Patient_Care_followup", "severity_max", "severity_grade_max"]
    subj = subj.join(fu)
    subj = subj.join(lab[["MaxClin_Follow_up", "PACS_6M_woDys", "PACS_12M"]])
    # the healthy samples of four people who later had COVID-19 (and a separate SubjectID '<n>_2')
    later = {s.replace("_2_ID", "_ID"): s for s in subj.index if s.endswith("_2_ID")}
    subj["later_covid_subject_id"] = [later.get(s) if subj.loc[s, "COVID"] == "Healthy" else None for s in subj.index]
    return subj


def sample_row(r: pd.Series, subj: pd.DataFrame, platform: str) -> dict:
    s = r["SubjectID"]
    m = subj.loc[s]
    healthy = m["COVID"] == "Healthy"
    month = r["Sampling_month"]
    row = {
        "subject_id": s,
        "sample_matrix": "serum",
        "platform": platform,
        "collection_country": "Switzerland",
        "age": r["Age"],
        "sex": {"male": "male", "female": "female"}.get(str(r["Sex"]).strip().lower(), "unknown"),
        "bmi": m["BMI"],
        "comorbidity_count": m["comorbidity_count"],
        "original_sample_id": r["SampleId"],
        "sampling_month": month,
        "vaccinated_before_sampling": r["Post_Vaccine"] == "YES",
        "covid19_severity_max": None if healthy else m["severity_max"],
        "covid19_severity_grade_max": None if healthy else m["severity_grade_max"],
        "pacs_6m_wodys": None if healthy else m["PACS_6M_woDys"],
        "pacs_12m": None if healthy else m["PACS_12M"],
        "max_follow_up_attended": m["MaxClin_Follow_up"],
    }
    notes = ["BMI and comorbidities recorded at the first study visit"]
    if healthy:
        row.update(
            lc_status="never_infected",
            case_definition="not_applicable",
            control_type="never_infected_contemporary",
            sample_phase="not_applicable",
            timepoint="healthy_first_visit",
            hospitalised=False,
        )
        notes.append("healthy control recruited Apr 2020 to Apr 2021, no history of SARS-CoV-2 infection, negative anti-S1 IgA and IgG serology")
        if pd.notna(m["later_covid_subject_id"]):
            row["later_covid_subject_id"] = m["later_covid_subject_id"]
            notes.append(
                f"same sample number, sex and age as COVID-19 subject {m['later_covid_subject_id']} (probably the same person sampled before infection); "
                "kept as the authors' healthy control"
            )
        return {**row, "notes": "; ".join(notes)}

    lc = m["PACS_6M_woDys"]
    status = "LC" if lc == 1 else "COVID_recovered" if lc == 0 else "unknown"
    days = r["Days"]
    asympt = m["severity_grade_max"] == "Asymptomatic"
    row.update(
        lc_status=status,
        case_definition="study_specific" if status == "LC" else "not_applicable",
        acute_severity=severity_band(m),
        hospitalised=bool(m["Patient_Care_followup"] == "Hospitalized" or m["ICU"] == "YES" or m["Oxygen"] == "YES"),
        vaccinated_before_infection=False if subj.loc[s, "Post_Vaccine"] == "NO" else None,
        days_since_infection=days if pd.notna(days) and days >= 0 else None,
        time_anchor="positive_test" if asympt else "symptom_onset",
    )
    if month == "ONE":
        row.update(sample_phase="acute", control_type="not_applicable", timepoint="acute")
        notes.append("acute sample; lc_status is the 6-month outcome (PACS_6M_woDys)")
    else:
        row.update(sample_phase="post_acute", control_type="not_applicable" if status == "LC" else "COVID_recovered", timepoint="6M")
    if status == "COVID_recovered":
        notes.append("PACS_6M_woDys 0: no symptoms at 6 months, or only isolated smell/taste change at 6 months (8 such patients, not identifiable in the deposited files), includes 'recovered Long Covid' (symptoms >4 weeks, gone by 6 months)")
    if status == "LC":
        notes.append("6-month Long Covid: >=1 persisting COVID-19 related symptom at 6-month follow-up not explained otherwise, isolated smell/taste change excluded")
    if asympt:
        notes.append("asymptomatic: days counted from positive PCR")
    return {**row, "notes": "; ".join(notes)}


def build_somascan(six, ac, apt, subj) -> dict:
    keep = human_aptamers(apt)
    seq = [c for c in apt["AptName"] if c in six.columns]
    D = pd.concat([ac, six], ignore_index=True).copy()
    D["sample_id"] = D["SampleId"].map(sid)
    if D["sample_id"].duplicated().any():
        raise ValueError("duplicated SomaScan sample ids")
    log10 = D.set_index("sample_id")[seq].apply(pd.to_numeric, errors="coerce")
    rfu = np.power(10.0, log10)  # deposited log10 RFU -> RFU; the adapter takes log2
    adapted = to_feature_space("somascan", rfu, annotation=apt, aptamer_col="AptName", uniprot_col="UniProt", keep=keep)
    rows = {r["sample_id"]: sample_row(r, subj, SOMA_PLATFORM) for _, r in D.iterrows()}
    meta = pd.DataFrame.from_dict(rows, orient="index").loc[adapted.X.index]
    meta.index.name = "sample_id"
    extra = {
        "n_aptamers": len(seq),
        "n_aptamers_kept_human": int(keep.reindex(seq).fillna(False).sum()),
        "transform": adapted.notes.get("transform"),
        "aptamers_without_uniprot": adapted.notes.get("aptamers_without_uniprot"),
    }
    return write_cohort(f"{STUDY}_SomaScan", adapted, meta, extra)


def build_ms(ms, six, subj) -> dict:
    cols = [c for c in ms.columns if c != "Protein.Group"]
    adapted = to_feature_space("ms_proteomics", ms, id_col="Protein.Group", sample_cols=cols, preferred_accessions=swissprot_accessions())
    six_by_sid = six.copy().assign(sample_id=six["SampleId"].map(sid)).set_index("sample_id")
    rows = {}
    for c in cols:
        if c in six_by_sid.index:
            rows[c] = sample_row(six_by_sid.loc[c], subj, MS_PLATFORM)
            rows[c]["in_somascan"] = True
            continue
        # MS columns without a SomaScan 6-month sample: no deposited visit or clinical data
        num = c.replace("COV_ZH_", "").split(".")[0]
        cand = [s for s in (f"COV_ZH_{num}_2_ID", f"COV_ZH_{num}_ID") if s in subj.index and subj.loc[s, "COVID"] == "COVID"]
        row = {
            "subject_id": cand[0] if cand else f"COV_ZH_{num}_ID",
            "lc_status": "unknown",
            "case_definition": "not_applicable",
            "sample_matrix": "serum",
            "platform": MS_PLATFORM,
            "collection_country": "Switzerland",
            "control_type": "not_applicable",
            "sample_phase": "unknown",
            "timepoint": "unknown",
            "original_sample_id": c,
            "in_somascan": False,
        }
        if cand:
            m = subj.loc[cand[0]]
            row.update(
                sex={"male": "male", "female": "female"}.get(str(m["Sex"]).lower(), "unknown"),
                pacs_6m_wodys=m["PACS_6M_woDys"],
                pacs_12m=m["PACS_12M"],
                covid19_severity_max=m["severity_max"],
                notes=(
                    f"extra MS sample of subject {cand[0]} (subject inferred from the sample number); visit and sampling day not deposited, "
                    f"so lc_status left unknown (subject PACS_6M_woDys {m['PACS_6M_woDys']})"
                ),
            )
        else:
            row["notes"] = "MS sample whose subject is absent from all deposited clinical and SomaScan tables; no label"
        rows[c] = row
    meta = pd.DataFrame.from_dict(rows, orient="index").loc[adapted.X.index]
    meta.index.name = "sample_id"
    extra = {
        "n_protein_groups": len(ms),
        "n_extra_ms_samples": int((~meta["in_somascan"].astype(bool)).sum()),
        "transform": adapted.notes.get("transform"),
        "groups_contaminant_or_decoy": adapted.notes.get("groups_contaminant_or_decoy"),
        "groups_without_uniprot": adapted.notes.get("groups_without_uniprot"),
    }
    return write_cohort(f"{STUDY}_MS", adapted, meta, extra)


def build() -> list:
    raw = raw_dir(STUDY)
    six, ac, lab, apt, ms = load(raw)
    subj = subject_table(six, ac, lab)
    return [build_somascan(six, ac, apt, subj), build_ms(ms, six, subj)]


if __name__ == "__main__":
    for s in build():
        print("----")
        for k, v in s.items():
            print(f"{k}: {v}")
