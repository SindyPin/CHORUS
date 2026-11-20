"""GSE215865 clinical predictors from the Mount Sinai COVID-19 Biobank table (Synapse syn35876735; Thompson et al.,
Nat Med 2023), for early detection (Amendment 5, B and D).

The table has one row per subject per day (Days_Since_First_Sample, 0 = the first biobank blood sample). Only rows on
or before the day of each participant's first blood sample used here (the earliest RNA-seq sample, the Round2 sample;
its Blood_Sample_ID is the sample id without the library plate suffix) are read:
  time invariant columns    the subject's value (recorded at consent or enrolment)
  medications (MED_*)       given on any day up to the sample day (true), recorded and never given (false)
  laboratory (LAB_*_Median) the last daily median on or before the sample day
  SOFA_Score                the maximum on or before the sample day
  Severity                  the worst on or before the sample day (Moderate 1, Severe 2, Severe with end organ damage 3)
  VENTILATION_TYPE          the highest support on or before the sample day, mapped to the WHO Clinical Progression
                            Scale (who_severity): room air 4, supplemental oxygen 5, non invasive ventilation 6,
                            invasive ventilation 7 (lower bound of 7 to 9)
  Weight_KG                 the last value on or before the sample day
  hospitalised              in the COVID-19 hospital encounter on the sample day
  vaccinated_before_infection  first vaccine dose on or before the sample day (all samples are from 2020)
  days_since_infection      days since reported symptom onset on the sample day (timing: sensitivity analysis only)
Excluded: the post COVID checklist (defines the label), anything about days after the sample (encounter level flags,
days until discharge, death), SARS-CoV-2 test results (infection status), identifiers and dates.
"""

import re

import numpy as np
import pandas as pd

TABLE = "data/clinical/GSE215865/clinical_data_table.csv"
COMMON = {"Age": "age", "Sex": "sex", "Race_From_Consent": "race", "Ethnicity_From_Consent": "ethnicity", "BMI": "bmi",
          "SMOKING_STATUS": "smoking"}
SEVERITY = {"Moderate COVID-19": 1, "Severe COVID-19": 2, "Severe COVID-19 with EOD": 3}
VENT_CPS = {"NONE_ROOM_AIR": 4, "SUPPLEMENTAL_O2": 5, "NON_INVASIVE_VENTILATION": 6, "INVASIVE_VENTILATION": 7}
MISSING = {"Unknown/Prefer not to say", "NOT ASKED", "Unreviewed", "Unclear"}


def columns(T: pd.DataFrame) -> dict:
    """Decision for every column of the table: (status, kind, type, reason, used as)."""
    D = {}
    for c in T.columns:
        if c.startswith("Post_COVID19_"):
            D[c] = ("excluded", "symptoms", "", "post COVID checklist item: defines the label", "")
        elif c in COMMON:
            D[c] = ("eligible", "demographic", "", f"time invariant; used as {COMMON[c]}", COMMON[c])
        elif re.fullmatch(r"LAB_.+_Median", c):
            D[c] = ("eligible", "acute illness", "numeric", "laboratory value: last daily median on or before the sample day", c)
        elif c.startswith("LAB_"):
            D[c] = ("excluded", "acute illness", "", "another summary of a laboratory value kept as its daily median", "")
        elif c.startswith("MED_"):
            D[c] = ("eligible", "acute illness", "binary", "acute treatment given on or before the sample day (not treatment for LC)", c)
        elif c == "COMORBID_CANCER_DIAGNOSIS_DESCRIPTION":
            D[c] = ("excluded", "comorbidity", "", "free text (cancer is in COMORBID_CANCER_FLAG)", "")
        elif c.startswith("COMORBID_"):
            D[c] = ("eligible", "comorbidity", "binary", "comorbidity (time invariant)", c)
        else:
            D[c] = OTHER.get(c, ("excluded", "unclassified", "", "not classified", ""))
    return D


OTHER = {
    "Height_CM": ("eligible", "demographic", "numeric", "time invariant", "Height_CM"),
    "Weight_KG": ("eligible", "demographic", "numeric", "last value on or before the sample day", "Weight_KG"),
    "SOFA_Score": ("eligible", "acute illness", "numeric", "maximum on or before the sample day", "SOFA_Score"),
    "Severity": ("eligible", "acute illness", "numeric", "worst on or before the sample day (Moderate 1, Severe 2, Severe with EOD 3)", "Severity"),
    "VENTILATION_TYPE": ("excluded", "acute illness", "", "highest support on or before the sample day, used as who_severity", "who_severity"),
    "COVID19_Symptoms_At_Presentation_Reviewed": ("eligible", "acute illness", "categorical",
                                                  "COVID-19 symptoms at hospital presentation (acute; not the post COVID checklist); Unreviewed and Unclear as missing",
                                                  "COVID19_Symptoms_At_Presentation_Reviewed"),
    "COVID19_Within_Encounter": ("excluded", "acute illness", "", "on the sample day: used as hospitalised", "hospitalised"),
    "COVID19_Days_Since_Symptoms_Reported_Onset": ("timing", "timing", "", "on the sample day: used as days_since_infection (sensitivity analysis only)", "days_since_infection"),
    "COVID19_Vaccine_Received": ("excluded", "vaccination", "", "with the dose dates: used as vaccinated_before_infection", "vaccinated_before_infection"),
    "COVID19_Vaccine_Dose_1_Date": ("excluded", "vaccination", "", "used for vaccinated_before_infection", ""),
    "COVID19_Vaccine_Dose_2_Date": ("excluded", "vaccination", "", "used for vaccinated_before_infection", ""),
    "COVID19_Vaccine_Doses": ("excluded", "vaccination", "", "doses at any time, including after the sample", ""),
    "COVID19_Vaccine_Type": ("excluded", "vaccination", "", "product of a vaccine given after the sample", ""),
    "Encounter_ICU_Flag": ("excluded", "acute illness", "", "ICU at any point of the encounter: includes days after the sample", ""),
    "COVID19_Days_Until_Encounter_End": ("excluded", "timing", "", "days until discharge: after the sample", ""),
    "DECEASED_INDICATOR": ("excluded", "other", "", "death: after the sample", ""),
    "COVID19_Deceased_Within_Encounter": ("excluded", "other", "", "death in the encounter: after the sample", ""),
    "Severity_Ignoring_COVID19_Status": ("excluded", "acute illness", "", "duplicate of Severity", ""),
    "Severity_Ambiguous_EOD": ("excluded", "acute illness", "", "data quality flag of Severity", ""),
    "COVID19_Days_Since_Encounter_Start": ("excluded", "timing", "", "days since admission at the sample: set by sampling; the timing sensitivity uses symptom onset", ""),
    "Days_Since_First_Sample_Encounter": ("excluded", "timing", "", "timing of the first biobank sample in its encounter: set by sampling", ""),
    "Encounter_Date_Of_First_Sample": ("excluded", "timing", "", "date", ""),
    "COVID19_Symptoms_Reported_Onset_Date": ("excluded", "timing", "", "date (onset days used in the sensitivity analysis)", ""),
    "COVID19_Symptoms_Reported_Onset_Is_Exact": ("excluded", "timing", "", "precision of the onset date", ""),
    "Subject_ID": ("excluded", "design", "", "identifier", ""),
    "Blood_Sample_ID": ("excluded", "design", "", "identifier", ""),
    "Days_Since_First_Sample": ("excluded", "timing", "", "the day index of the table (used to select rows on or before the sample day)", ""),
}
for _c in ("Patient_Classification_At_First_Sample", "COVID19_Order_Result", "COVID19_Any_Antibody_Detected", "COVID19_Days_Since_Most_Recent_Antibody",
           "COVID19_Days_Since_Most_Recent_Order", "COVID19_Ever_Positive", "COVID19_Most_Recent_Antibody_Date", "COVID19_Most_Recent_Antibody_Result",
           "COVID19_Most_Recent_Order_Date", "COVID19_Most_Recent_Order_Result", "COVID19_No_Test_Result", "COVID19_PCR_First_Positive_Date",
           "COVID19_PCR_First_Test_Date", "SARSCoV2_Recovery_Status"):
    OTHER[_c] = ("excluded", "infection status", "", "SARS-CoV-2 test result or date: infection status (the cohort's eligibility), partly after the sample", "")


def boolean(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.upper().map({"TRUE": 1.0, "FALSE": 0.0})


def predictors(code, chosen: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """One row per Round2 subject (subject_id as in Round2) with the variables above; and the column decisions."""
    T = pd.read_csv(code / TABLE, low_memory=False)
    dec = columns(T)
    day_of = T.dropna(subset=["Blood_Sample_ID"]).set_index("Blood_Sample_ID")["Days_Since_First_Sample"]
    out = []
    for r in chosen.itertuples():
        subj = str(r.subject_id).replace("GSE215865:", "")
        bsid = re.sub(r"_Plate_\d+$", "", str(r.sample_id))
        cut = day_of.get(bsid, np.nan)
        g = T[T["Subject_ID"] == subj]
        before = g[g["Days_Since_First_Sample"] <= cut] if np.isfinite(cut) else g.iloc[0:0]
        row = {"subject_id": r.subject_id, "sample_day": cut}
        inv = g.iloc[0] if len(g) else pd.Series(dtype=object)
        for c, name in COMMON.items():
            v = inv.get(c, np.nan)
            row[name] = np.nan if (pd.isna(v) or str(v) in MISSING) else v
        for c, (status, kind, typ, why, used) in dec.items():
            if status not in ("eligible",) or c in COMMON:
                continue
            if c.startswith("COMORBID_"):
                row[c] = boolean(pd.Series([inv.get(c, np.nan)])).iloc[0]
            elif c == "Height_CM":
                row[c] = inv.get(c, np.nan)
            elif c.startswith("MED_"):
                b = boolean(before[c]).dropna()
                row[c] = b.max() if len(b) else np.nan
            elif c.startswith("LAB_") or c == "Weight_KG":
                x = before[["Days_Since_First_Sample", c]].dropna()
                row[c] = x.sort_values("Days_Since_First_Sample")[c].iloc[-1] if len(x) else np.nan
            elif c == "SOFA_Score":
                x = before[c].dropna()
                row[c] = x.max() if len(x) else np.nan
            elif c == "Severity":
                x = before[c].map(SEVERITY).dropna()
                row[c] = x.max() if len(x) else np.nan
            elif c == "COVID19_Symptoms_At_Presentation_Reviewed":
                v = inv.get(c, np.nan)
                row[c] = np.nan if (pd.isna(v) or str(v) in MISSING) else v
        x = before["VENTILATION_TYPE"].map(VENT_CPS).dropna()
        row["who_severity"] = x.max() if len(x) else np.nan
        on_day = g[g["Days_Since_First_Sample"] == cut]
        w = boolean(on_day["COVID19_Within_Encounter"]).dropna()
        row["hospitalised"] = w.iloc[0] if len(w) else np.nan
        d1 = inv.get("COVID19_Vaccine_Dose_1_Date", np.nan)
        rec = str(inv.get("COVID19_Vaccine_Received", "")).upper()
        row["vaccinated_before_infection"] = (1.0 if (pd.notna(d1) and d1 <= cut) else 0.0) if (pd.notna(d1) or rec == "FALSE") else np.nan
        on = on_day["COVID19_Days_Since_Symptoms_Reported_Onset"].dropna()
        row["days_since_infection"] = on.iloc[0] if len(on) else np.nan
        out.append(row)
    return pd.DataFrame(out), dec
