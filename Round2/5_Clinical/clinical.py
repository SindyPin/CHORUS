"""Shared code of Round2/5_Clinical (PROTOCOL.md, Amendment 5, B and C).

Clinical variables as predictors of Long COVID, for 3_Detection (post acute samples) and 2_EarlyDetection (acute
samples, outcome known later). Labels and sample phases come from Round2/0_Labels (sample_labels.csv); models,
folds, metrics and the resumable runner from Round2/0_Controls/common.py.

Units: a unit is one group of people. The two layers of Su2022 (Olink and metabolomics) are the same people, and so
are the two layers of Cervia2024 (SomaScan and MS): each is one unit, so that every person is counted once.
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0_Controls"))
import common as K  # noqa: E402

CODE = K.CODE
LABELS = CODE / "Round2/0_Labels/output/sample_labels.csv"
CACHE = HERE / "cache"
TASKS = {"3_Detection": "post_acute", "2_EarlyDetection": "acute"}
MIN_RECORDED = 0.5  # Amendment 5, B4
UNITS = {"Su2022": ["Su2022", "Su2022_metabolomics"], "Cervia2024": ["Cervia2024_SomaScan", "Cervia2024_MS"]}
UNIT_OF = {s: u for u, ss in UNITS.items() for s in ss}
FLAGGED = {"GSE270045", "GSE228320", "CORSIP2026", "Buggert2025"}  # PROTOCOL.md and Amendment 1, item 1
SPLITS = {"Buggert2025": ("collection_country", {"United Kingdom": "Buggert2025_UK", "Sweden": "Buggert2025_SE"})}
POSSIBLE_OVERLAP = {frozenset({"Chowdhury2024", "Chowdhury2026"}): "same laboratory and biobank; overlap not resolvable (Round1 metadata)"}

# ---------------------------------------------------------------------------------------------------- variables
# Every column recorded about a participant gets one decision. status: eligible (main analysis), early (eligible only
# for early detection), timing (sensitivity analysis only), excluded. kind: what the variable describes.
# type: numeric, binary (0/1) or categorical. unit: the unit of the odds ratio in 01_univariate.
WHO_CPS = {  # harmonised acute severity (Round1 schema) -> WHO Clinical Progression Scale (Marshall 2020, 0 to 10)
    "ambulatory": (2, "ambulatory mild disease (scale 1 to 3): 2, symptomatic and independent"),
    "hospitalised_no_oxygen": (4, "hospitalised, no oxygen therapy: 4"),
    "hospitalised_oxygen": (5, "hospitalised, oxygen by mask or nasal prongs: 5 (high flow or non invasive ventilation, 6, not separable)"),
    "icu_or_ventilated": (7, "intensive care or ventilated: 7, the lowest ventilated level (ICU without ventilation, 6, not separable)"),
}
WHO8_CPS = {  # GSE224615 Supplementary Table 3, WHO 8 point ordinal scale (R&D Blueprint 2020) -> Clinical Progression Scale
    "1": (2, "WHO ordinal 1, ambulatory without limitation: 2"),
    "2": (3, "WHO ordinal 2, ambulatory with limitation: 3"),
    "3": (4, "WHO ordinal 3, hospitalised without oxygen: 4"),
    "4-5": (5, "WHO ordinal 4 or 5, oxygen by mask or high flow: 5 (lower bound)"),
}

V = {}


def var(name, status, kind, vtype, reason, unit=""):
    V[name] = dict(status=status, kind=kind, type=vtype, reason=reason, unit=unit)


# eligible (Amendment 5, B1)
var("age", "eligible", "demographic", "numeric", "demographic", "per 10 years")
var("sex", "eligible", "demographic", "binary", "demographic (female = 1)", "female vs male")
var("sex_code", "eligible", "demographic", "categorical", "demographic; Lopez2023 SEX coded 1 or 2, coding undocumented")
var("bmi", "eligible", "demographic", "numeric", "demographic", "per 5 kg/m2")
var("race", "eligible", "demographic", "categorical", "demographic")
var("ethnicity", "eligible", "demographic", "categorical", "demographic")
var("who_severity", "eligible", "acute illness", "numeric", "acute severity mapped to the WHO Clinical Progression Scale (who_mapping.csv)", "per point")
var("hospitalised", "eligible", "acute illness", "binary", "acute hospitalisation", "yes vs no")
var("icu", "eligible", "acute illness", "binary", "intensive care or ventilation during the acute illness (from acute_severity)", "yes vs no")
var("disease_status_subgroup", "excluded", "acute illness", "", "disease status subgroup (Bansal2025): may define the label (Amendment 7, item 4)")
var("dexamethasone", "eligible", "acute illness", "categorical", "treatment of the acute illness, not of LC")
var("comorbidity_count", "eligible", "comorbidity", "numeric", "comorbidity", "per condition")
var("comorbidity", "eligible", "comorbidity", "binary", "comorbidity present (Bansal2025)", "present vs absent")
var("n_preexisting_conditions", "eligible", "comorbidity", "numeric", "comorbidity: number of conditions listed in preexisting_conditions_cc (GSE226260_WB)", "per condition")
var("smoking", "eligible", "other", "categorical", "smoking; Lopez2023 coded 1 or 2, coding undocumented")
var("t2dm", "eligible", "comorbidity", "categorical", "comorbidity; Lopez2023 coded 1 or 2, coding undocumented")
var("hypertension", "eligible", "comorbidity", "categorical", "comorbidity; Lopez2023 coded 1 or 2, coding undocumented")
var("obesity", "eligible", "comorbidity", "categorical", "comorbidity; Lopez2023 coded 1 or 2, coding undocumented")
var("vaccinated_before_infection", "eligible", "vaccination", "binary", "vaccination before infection", "yes vs no")
var("infection_wave", "eligible", "infection", "categorical", "infection wave")
var("variant_era", "eligible", "infection", "categorical", "infection variant era")
var("pregnancy", "eligible", "other", "categorical", "participant characteristic")
# state at the draw: known at the acute stage for early detection; at a post acute draw it describes the same time as the label
var("wos_at_draw", "early", "acute illness", "categorical", "WHO ordinal scale at the draw: acute state for early detection; at a post acute draw it is measured with the label, so excluded for detection")
var("encounter_location", "early", "acute illness", "categorical", "place of the draw (home, clinic, hospital, ICU): as wos_at_draw")
# timing (Amendment 5, B2): sensitivity analysis only
var("days_since_infection", "timing", "timing", "numeric", "time since infection: often set by recruitment design; sensitivity analysis only", "per 30 days")
var("days_since_infection_at_sample", "timing", "timing", "numeric", "time since infection (GSE226260_WB TPI): sensitivity analysis only", "per 30 days")
var("sampling_month", "excluded", "timing", "categorical", "follow up visit of the draw (one or six months): set by the study design")
var("days_since_symptom_resolution", "excluded", "timing", "numeric", "timing after the acute illness: set by the design")
var("timepoint", "excluded", "timing", "categorical", "study visit of the draw: set by the design")
var("time_anchor", "excluded", "timing", "categorical", "reference point of the timing columns: not about the participant")
var("timepoint.1", "excluded", "timing", "categorical", "duplicate of timepoint (study visit): set by the design")
var("lc_definition_strict_who", "excluded", "symptoms", "", "LC definition flag: defines the label")
# excluded: define the label or come after it (Amendment 5, B2)
for v_ in ("n_pasc_symptoms", "pasc_symptoms", "pasc_groups", "pasc_count_group", "symptom_domains", "symptom_count", "symptom_list",
           "composite_binary", "composite_binary_rec", "pacs_6m_wodys", "pacs_12m", "post_covid_symptom_severity_score", "dyspnoea_item"):
    var(v_, "excluded", "symptoms", "", "symptoms, symptom counts or scores: define the label")
var("lopez_symptoms", "excluded", "symptoms", "", "19 symptom columns (Lopez2023): define the label; coding undocumented")
for v_ in ("borg_cr10", "dyspnoea12", "nijmegen"):
    var(v_, "excluded", "score", "", "questionnaire score used with the LC definition")
var("dlco_pct_predicted", "excluded", "score", "", "outcome measure used as the case definition (GSE228320)")
var("lc_subgroup", "excluded", "symptoms", "", "LC subgroup: defines the label")
var("pasc_group_binary_xgboost", "excluded", "symptoms", "", "LC group: defines the label")
var("group_hr", "excluded", "symptoms", "", "study group (PASC, recovered, ICU): defines the label")
var("intervention_arm", "excluded", "treatment", "", "treatment given for LC (exercise trial)")
var("vaccinated_before_sampling", "excluded", "vaccination", "", "includes vaccination after infection: may follow the label")
var("vaccine_product", "excluded", "vaccination", "", "timing relative to infection not stated: may be after infection")
var("max_follow_up_attended", "excluded", "other", "", "follow up attendance: comes after the label")
# excluded: same information as another variable (one severity scale across cohorts)
var("acute_severity", "excluded", "acute illness", "", "mapped to who_severity, hospitalised and icu")
var("who_category_acute", "excluded", "acute illness", "", "mapped to who_severity (GSE224615)")
for v_ in ("peak_acute_wos", "covid19_severity_max", "covid19_severity_grade_max"):
    var(v_, "excluded", "acute illness", "", "same information as who_severity (one severity scale across cohorts)")
var("bmi_table_s2", "excluded", "demographic", "", "duplicate of bmi")
var("preexisting_conditions_cc", "excluded", "comorbidity", "", "free text list: used as n_preexisting_conditions")
var("inferred_sex", "excluded", "demographic", "", "sex inferred from expression, not recorded")
var("inferred_sex_method", "excluded", "other", "", "method of inferred_sex")
var("sex_code_deposited", "excluded", "demographic", "", "raw code behind sex_code")
var("weight", "excluded", "demographic", "", "flagged cohort only (CORSIP2026)")
var("height", "excluded", "demographic", "", "flagged cohort only (CORSIP2026)")
var("hospitalised_ever_18m", "excluded", "acute illness", "", "flagged cohort only (CORSIP2026); 18 month window")
var("vaccinated", "excluded", "vaccination", "", "flagged cohort only (CORSIP2026); timing not stated")
# excluded: about the study, the sample or the assay, not the participant
for v_ in ("control_type", "collection_country", "cohort_arm", "study_group", "case_definition", "lc_status", "lc_status_original",
           "sample_matrix", "omics_layer", "platform", "sample_phase", "notes", "label_source", "control_definition_status",
           "linked_healthy_subject_id", "overlaps_with_study", "in_somascan", "group_code", "group_prefix", "group_deposited",
           "table_s1_label", "collection_date"):
    var(v_, "excluded", "design", "", "study design or label column, not a participant characteristic")
# GSE215865: the Mount Sinai COVID-19 Biobank table (Synapse; synapse_gse215865.py). Its eligible columns keep their
# names; numeric ones get odds ratios per standard deviation in 01_univariate.
import synapse_gse215865 as SYN  # noqa: E402

_T = pd.read_csv(CODE / SYN.TABLE, nrows=0) if (CODE / SYN.TABLE).exists() else pd.DataFrame()
SYNAPSE_DECISIONS = SYN.columns(_T) if len(_T.columns) else {}
for _c, (_st, _kind, _typ, _why, _used) in SYNAPSE_DECISIONS.items():
    if _st == "eligible" and _used == _c:
        var(_c, "eligible", _kind, _typ, f"Synapse table: {_why}", "per SD" if _typ == "numeric" else ("yes vs no" if _typ == "binary" else ""))
TECHNICAL_HINTS = ("plate", "olink", "panel", "explore", "raw", "ms_", "precursors", "proteins_identified", "sdrf", "value_type",
                   "diann", "mq_version", "fasta", "file", "batch", "salmon", "biosample", "gsm", "original_sample_id", "repeated",
                   "n_values_qc", "run")

# acute illness and baseline variables known at the acute stage (Amendment 5, B3)
ACUTE_KINDS = {"demographic", "acute illness", "comorbidity", "vaccination", "infection", "other"}

# linked open tables (written by Round2/0_Labels/build_labels.py): cohort -> {variable: column}
LINKED = {
    "GSE224615": {"age": "age", "race": "race", "hospitalised": "hospitalized", "who_category_acute": "who_category_acute"},
    "GSE226260_WB": {"age": "Age", "sex": "Sex at birth", "race": "Race", "ethnicity": "Ethnicity",
                     "preexisting_conditions_cc": "preexisting_conditions_cc", "days_since_infection_at_sample": "TPI",
                     "symptom_count": "Symptom_count", "symptom_list": "symptom_list"},
    "GSE226260_PBMC": {"age": "Age (12/31/22)", "sex": "Gender", "race": "Race", "ethnicity": "Ethnicity", "vaccine_product": "Vaccine product"},
    "Lopez2023": {"age": "AGE", "sex_code": "SEX", "smoking": "SMOKING", "t2dm": "T2DM", "hypertension": "ARTERIAL HYPERTENSION",
                  "obesity": "OBESITY", "lopez_symptoms": "FEVER"},
    "Fineschi2026": {"infection_wave": "Wave", "post_covid_symptom_severity_score": "Score", "dyspnoea_item": "Dyspnoea"},
}
LINKED_FILE = {"GSE226260_WB": "data/clinical/GSE226260_WB/GSE226260_WB_clinical_linked.csv",
               "GSE226260_PBMC": "data/clinical/GSE226260_PBMC/GSE226260_PBMC_clinical_linked.csv"}
MISSING = {"", "nan", "none", "unknown", "missing", "na", "not_applicable", "none listed in epic"}


# Amendment 7, item 3 (early detection): a variable is in the main analysis only if it is known at or before the
# acute sample used. Peak severity, ICU over the whole stay, and anything whose timing cannot be established are used
# only in the sensitivity analysis "known by the end of the acute illness".
AT_SAMPLE, END_OF_ACUTE = "known at or before the acute sample", "known by the end of the acute illness"
BASELINE = {"age", "sex", "sex_code", "race", "ethnicity", "bmi", "smoking", "comorbidity_count", "comorbidity",
            "n_preexisting_conditions", "t2dm", "hypertension", "obesity", "infection_wave", "variant_era",
            "vaccinated_before_infection"}
AT_DRAW = {"wos_at_draw": "WHO ordinal scale on the day of the draw", "encounter_location": "place of the draw"}
TIMING = {  # (unit, variable): (timing, reason)
    ("Su2022", "who_severity"): (END_OF_ACUTE, "from peak_acute_wos: peak severity of the acute illness"),
    ("Su2022", "icu"): (END_OF_ACUTE, "from peak_acute_wos: ICU or ventilation at any time in the acute illness"),
    ("Su2022", "hospitalised"): (END_OF_ACUTE, "from peak_acute_wos: hospitalised at any time in the acute illness"),
    ("Cervia2024", "who_severity"): (END_OF_ACUTE, "from covid19_severity_max: peak severity of the acute illness"),
    ("Cervia2024", "icu"): (END_OF_ACUTE, "from covid19_severity_max: ICU or ventilation at any time in the acute illness"),
    ("Cervia2024", "hospitalised"): (END_OF_ACUTE, "hospitalised at any time in the acute illness; date relative to the sample not given"),
    ("GSE215865", "who_severity"): (AT_SAMPLE, "Synapse: highest ventilation support on or before the sample day"),
    ("GSE215865", "hospitalised"): (AT_SAMPLE, "Synapse: in the COVID-19 hospital encounter on the sample day"),
}


def timing_of(name: str, unit: str) -> tuple[str, str]:
    """Early detection timing of a variable in a unit (Amendment 7, item 3)."""
    if (unit, name) in TIMING:
        return TIMING[(unit, name)]
    if name in BASELINE:
        return AT_SAMPLE, ("baseline characteristic recorded before or at enrolment" if name not in ("comorbidity_count", "comorbidity", "n_preexisting_conditions", "t2dm", "hypertension", "obesity")
                           else "pre-existing condition")
    if name in AT_DRAW:
        return AT_SAMPLE, AT_DRAW[name]
    if unit == "GSE215865" and name in SYNAPSE_DECISIONS:
        return AT_SAMPLE, "Synapse: time invariant, or rows on or before the sample day only (synapse_gse215865.py)"
    return END_OF_ACUTE, "timing relative to the sample cannot be established"


def technical(col: str) -> bool:
    """Assay or file columns (by name); a column with a decision of its own is never technical ("raw" is in "wos_at_draw")."""
    return col not in V and any(h in col.lower() for h in TECHNICAL_HINTS)


def eligible(name: str, task: str) -> bool:
    d = V[name]
    if task == "2_EarlyDetection":
        return d["status"] in ("eligible", "early") and d["kind"] in ACUTE_KINDS
    return d["status"] == "eligible"


# ------------------------------------------------------------------------------------------------- participants

def labels() -> pd.DataFrame:
    S = pd.read_csv(LABELS)
    S["sample_id"], S["subject_id"] = S["sample_id"].astype(str), S["subject_id"].astype(str)
    return S


def meta_of(study: str) -> pd.DataFrame:
    """Sample metadata of a study (split cohorts selected by country)."""
    for parent, (col, kids) in SPLITS.items():
        for val, kid in kids.items():
            if study == kid:
                m = pd.read_csv(CODE / f"data/meta/{parent}_meta.csv", low_memory=False)
                return m[m[col].astype(str) == val].reset_index(drop=True)
    return pd.read_csv(CODE / f"data/meta/{study}_meta.csv", low_memory=False)


def task_samples(study: str, task: str, S: pd.DataFrame) -> pd.DataFrame:
    """Samples of a study in the task's phase with an LC or recovered label, with y (LC = 1)."""
    parent = next((p for p, (_, kids) in SPLITS.items() if study in kids.values()), study)
    m = meta_of(study)
    m["sample_id"], m["subject_id"] = m["sample_id"].astype(str), m["subject_id"].astype(str)
    s = S[(S["study_id"] == parent) & (S["sample_phase_round2"] == TASKS[task]) & S["lc_status"].isin(["LC", "COVID_recovered"])]
    s = s[s["sample_id"].isin(set(m["sample_id"]))]
    out = m.merge(s[["sample_id", "lc_status"]].rename(columns={"lc_status": "lab"}), on="sample_id")
    out["y"] = (out["lab"] == "LC").astype(int)
    return out


def studies_for(task: str, S: pd.DataFrame) -> list:
    """Unflagged studies with LC or recovered samples in the task's phase (split cohorts by their parts)."""
    out = []
    for st in sorted(S["study_id"].unique()):
        if not (CODE / f"data/meta/{st}_meta.csv").exists():
            continue
        for k in (list(SPLITS[st][1].values()) if st in SPLITS else [st]):
            if k not in FLAGGED and len(task_samples(k, task, S)):
                out.append(k)
    return out


def norm(x) -> str | float:
    s = str(x).strip()
    return np.nan if s.lower() in MISSING else s


def raw_table(unit_studies: list, task: str, S: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per participant of a unit: y, the chosen sample of every layer (earliest in the phase), and every
    participant column of the metadata and linked tables (raw values; decisions applied later)."""
    rows, chosen = [], []
    for st in unit_studies:
        m = task_samples(st, task, S)
        if not len(m):
            continue
        first = K.one_per_participant(m)
        chosen.append(first[["subject_id", "sample_id", "y"]].assign(study=st))
        cols = [c for c in m.columns if c not in ("sample_id", "subject_id", "y", "lab", "study_id") and not technical(c)]
        # participant level: the chosen sample's value, else the participant's first recorded value
        for sid, g in m.groupby("subject_id"):
            g = pd.concat([first[first["subject_id"] == sid], g])
            r = {"subject_id": sid, "study": st, "y": int(g["y"].iloc[0])}
            for c in cols:
                vals = g[c].map(norm).dropna()
                r[c] = vals.iloc[0] if len(vals) else np.nan
            rows.append(r)
        if st == "GSE215865" and SYNAPSE_DECISIONS:
            P, _ = SYN.predictors(CODE, first)
            P = P.set_index("subject_id")
            for r in rows:
                if r["study"] == st and r["subject_id"] in P.index:
                    for c, v in P.loc[r["subject_id"]].items():
                        if c != "sample_day" and pd.notna(v):
                            r[c] = v
        if st in LINKED:
            f = CODE / LINKED_FILE.get(st, f"data/clinical/{st}/{st}_clinical_linked.csv")
            T = pd.read_csv(f, low_memory=False)
            T["subject_id"] = T["subject_id"].astype(str)
            for v, col in LINKED[st].items():
                vals = T[["subject_id", col]].assign(**{col: T[col].map(norm)}).dropna().drop_duplicates("subject_id").set_index("subject_id")[col]
                for r in rows:
                    if r["study"] == st and r["subject_id"] in vals.index:
                        r[v] = vals[r["subject_id"]]
    R = pd.DataFrame(rows)
    if R.empty:
        return R, pd.DataFrame()
    # a unit's layers are the same people: one row per person (values from the first layer that has them)
    R = R.groupby("subject_id", sort=True).agg(lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan).reset_index()
    return derive(R), pd.concat(chosen, ignore_index=True)


def derive(R: pd.DataFrame) -> pd.DataFrame:
    """Derived variables: who_severity, icu, hospitalised as 0/1, n_preexisting_conditions."""
    R = R.copy()
    if "acute_severity" in R:
        w = R["acute_severity"].map(lambda s: WHO_CPS.get(s, (np.nan,))[0])
        R["who_severity"] = R["who_severity"].fillna(w) if "who_severity" in R else w
        R["icu"] = R["acute_severity"].map(lambda s: np.nan if s not in WHO_CPS else float(s == "icu_or_ventilated"))
    if "who_category_acute" in R:
        w = R["who_category_acute"].map(lambda s: WHO8_CPS.get(str(s), (np.nan,))[0])
        R["who_severity"] = R["who_severity"].fillna(w) if "who_severity" in R else w
    if "preexisting_conditions_cc" in R:
        R["n_preexisting_conditions"] = R["preexisting_conditions_cc"].map(lambda s: np.nan if pd.isna(s) else len([x for x in str(s).split(",") if x.strip()]))
    return R


def to_binary(name: str, x: pd.Series) -> pd.Series:
    s = x.astype(str).str.strip().str.lower()
    if name == "sex":
        return s.map({"female": 1.0, "male": 0.0})
    if name == "comorbidity":
        return s.map({"present": 1.0, "absent": 0.0})
    return s.map({"true": 1.0, "false": 0.0, "yes": 1.0, "no": 0.0, "1": 1.0, "0": 0.0, "1.0": 1.0, "0.0": 0.0})


def typed(name: str, x: pd.Series) -> pd.Series:
    t = V[name]["type"]
    if t == "numeric":
        return pd.to_numeric(x, errors="coerce")
    if t == "binary":
        return to_binary(name, x)
    return x.map(lambda s: np.nan if pd.isna(s) else str(s))


def decisions(R: pd.DataFrame, unit: str, task: str) -> list:
    """One row per column: status, reason, recorded fraction, kept in this unit and task."""
    out = []
    n = len(R)
    for c in [c for c in R.columns if c not in ("subject_id", "study", "y")]:
        d = V.get(c)
        if d is None:
            out.append(dict(task=task, unit=unit, variable=c, status="excluded", kind="unclassified", type="",
                            reason="not a clinical variable of the protocol list (left out; see README)", n=n,
                            n_recorded=int(R[c].notna().sum()), missing_rate=round(1 - R[c].notna().mean(), 3), kept=False,
                            timing="", timing_reason="", kept_end_of_acute_sensitivity=False))
            continue
        x = typed(c, R[c]) if d["type"] else R[c]
        rec = x.notna()
        why = d["reason"]
        kept = False
        if d["status"] == "timing":
            kept_t = rec.mean() >= MIN_RECORDED and x[rec].nunique() > 1
            why = f"{why}; {'used in the timing sensitivity analysis' if kept_t else 'not used: ' + ('recorded for under 50 percent' if rec.mean() < MIN_RECORDED else 'constant')}"
        elif not eligible(c, task):
            why = why if d["status"] != "early" else f"{why} (detection: excluded)"
            if task == "2_EarlyDetection" and d["kind"] not in ACUTE_KINDS and d["status"] == "eligible":
                why = "not known at the acute stage"
        elif rec.mean() < MIN_RECORDED:
            why = f"eligible; dropped: recorded for {rec.mean():.0%} of participants (under 50 percent)"
        elif x[rec].nunique() < 2:
            why = "eligible; dropped: constant in this unit"
        else:
            kept = True
            why = "eligible; kept"
        tm, tr = ("not applicable (post acute task)", "") if task != "2_EarlyDetection" else timing_of(c, unit)
        kept_end = kept  # sensitivity "known by the end of the acute illness": every variable that passes the rules
        if kept and task == "2_EarlyDetection" and tm == END_OF_ACUTE:
            kept = False
            why = f"eligible; main analysis: excluded ({tr}); used only in the sensitivity analysis '{END_OF_ACUTE}' (Amendment 7, item 3)"
        out.append(dict(task=task, unit=unit, variable=c, status=d["status"], kind=d["kind"], type=d["type"], reason=why, n=n,
                        n_recorded=int(rec.sum()), missing_rate=round(1 - rec.mean(), 3), kept=kept,
                        timing=tm if d["status"] in ("eligible", "early") else "", timing_reason=tr if d["status"] in ("eligible", "early") else "",
                        kept_end_of_acute_sensitivity=kept_end if task == "2_EarlyDetection" else False))
    return out


def kept_vars(D: pd.DataFrame, unit: str, task: str, end_of_acute: bool = False) -> list:
    """Kept variables of a unit; end_of_acute: the early detection sensitivity set (Amendment 7, item 3)."""
    d = D[(D["unit"] == unit) & (D["task"] == task) & D["kept_end_of_acute_sensitivity" if end_of_acute else "kept"]]
    return list(d["variable"])


def timing_var(D: pd.DataFrame, unit: str, task: str) -> str | None:
    d = D[(D["unit"] == unit) & (D["task"] == task) & (D["status"] == "timing") & D["reason"].str.contains("used in the timing")]
    return d["variable"].iloc[0] if len(d) else None


# ---------------------------------------------------------------------------------------------- encoding

def arrays(R: pd.DataFrame, variables: list) -> tuple[np.ndarray, np.ndarray, list, list]:
    """(numeric block, categorical block as codes, numeric names, categorical names). Category codes are fixed per
    unit from the sorted category names (no labels used); one hot encoding happens inside the training fold."""
    num, cat, nn, cn = [], [], [], []
    for v in variables:
        x = typed(v, R[v]) if v in R else pd.Series(np.nan, index=R.index)
        if V[v]["type"] in ("numeric", "binary"):
            num.append(x.to_numpy(dtype=float))
            nn.append(v)
        else:
            cats = sorted(x.dropna().unique())
            cat.append(x.map({c: i for i, c in enumerate(cats)}).to_numpy(dtype=float))
            cn.append(v + "=" + "|".join(cats))
    n = len(R)
    N = np.column_stack(num) if num else np.zeros((n, 0))
    C = np.column_stack(cat) if cat else np.zeros((n, 0))
    return N, C, nn, cn


def encode(Ntr, Ctr, Nte, Cte):
    """Inside the training fold (Amendment 5, B4): numeric gaps filled with the training median, categorical gaps with
    the most frequent training category, a missing indicator for every variable with gaps in the training fold, one hot
    encoding on the training categories, columns constant in the training fold dropped, standardisation with training
    statistics. No feature selection (B5). B_scaled z scores the numeric columns within each cohort before this
    (run_fold)."""
    tr_cols, te_cols = [], []
    for j in range(Ntr.shape[1]):
        x, z = Ntr[:, j], Nte[:, j]
        miss = np.isnan(x)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            med = np.nanmedian(x) if (~miss).any() else 0.0
        tr_cols.append(np.where(miss, med, x))
        te_cols.append(np.where(np.isnan(z), med, z))
        if miss.any():
            tr_cols.append(miss.astype(float))
            te_cols.append(np.isnan(z).astype(float))
    for j in range(Ctr.shape[1]):
        x, z = Ctr[:, j], Cte[:, j]
        miss = np.isnan(x)
        if (~miss).sum() == 0:
            continue
        vals, counts = np.unique(x[~miss], return_counts=True)
        mode = vals[np.argmax(counts)]
        xf, zf = np.where(miss, mode, x), np.where(np.isnan(z), mode, z)
        for v in vals:
            tr_cols.append((xf == v).astype(float))
            te_cols.append((zf == v).astype(float))
        if miss.any():
            tr_cols.append(miss.astype(float))
            te_cols.append(np.isnan(z).astype(float))
    if not tr_cols:
        raise ValueError("no clinical column left in the training fold")
    A, B = np.column_stack(tr_cols), np.column_stack(te_cols)
    keep = A.std(0) > 0
    A, B = A[:, keep], B[:, keep]
    if not A.shape[1]:
        raise ValueError("every clinical column is constant in the training fold")
    mu, sd = A.mean(0), A.std(0)
    return (A - mu) / sd, (B - mu) / sd


def run_fold(job, d, tr, te):
    """Arms of 5_Clinical (called by common.run_fold through job["module"])."""
    task, name, arm = job["task"], job["model"], job["arm"]
    y = d["y"]
    N, C = d["N"], d["C"]
    if arm in ("clinical", "A_all", "A_shared"):
        A, B = encode(N[tr], C[tr], N[te], C[te])
        return K.fit_predict(name, task, A, y[tr], B)
    if arm in ("B_naive", "B_scaled"):
        Np, Cp, yp = d["Np"], d["Cp"], d["yp"]
        Ntr, Ctr, Nte = np.vstack([N[tr], Np]), np.vstack([C[tr], Cp]), N[te]
        if arm == "B_scaled" and N.shape[1]:
            Ntr = np.vstack([K.zscore_within(N[tr]), K.zscore_within(Np)])
            Nte = K.zscore_within(N[te], ref=N[tr])
        A, B = encode(Ntr, Ctr, Nte, C[te])
        return K.fit_predict(name, task, A, np.concatenate([y[tr], yp]), B)
    if arm == "omics":
        A, B, _ = K.prepare(d["X"][tr], y[tr], d["X"][te], task, job["k"])
        return K.fit_predict(name, task, A, y[tr], B)
    if arm == "joined":
        A1, B1, _ = K.prepare(d["X"][tr], y[tr], d["X"][te], task, job["k"])
        A2, B2 = encode(N[tr], C[tr], N[te], C[te])
        return K.fit_predict(name, task, np.hstack([A1, A2]), y[tr], np.hstack([B1, B2]))
    raise KeyError(arm)


def models() -> list:
    """Every protocol model except the expression models (they embed RNA and have no clinical or joined form, and
    every comparison here pairs the same model across arms); Mitra stays listed so that its removal is recorded."""
    return K.CLASSICAL + K.ENSEMBLE + K.FOUNDATION
