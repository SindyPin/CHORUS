"""Build the Zoodsma2022 Olink cohorts (Zoodsma M, de Nooijer AH, Grondman I, ..., Netea MG, Li Y 2022
Front Immunol, "Targeted proteomics identifies circulating biomarkers associated with active COVID-19
and post-COVID-19", doi:10.3389/fimmu.2022.1027122, PMID 36405747, PMC9670186). Data: GitHub
CiiM-Bioinformatics-group/Postcovid-targeted-proteomics at commit d45e3607 (MIT licence).

Three cohorts, three separate Olink NPX Manager projects, so three CHORUS cohorts:

  Zoodsma2022_Hannover  cohort 3 (MHH): post-COVID-19 convalescents and healthy blood donors;
                        Inflammation, Cardiometabolic, Cardiovascular II, Neurology (one run, 3 plates)
  Zoodsma2022_Breda     cohort 1 (Amphia Hospital Breda): hospitalised COVID-19, longitudinal;
                        Inflammation only (one run, 7 plates C4-C10)
  Zoodsma2022_Nijmegen  cohort 2 (Radboudumc, Janssen et al. 2021 J Infect Dis): hospitalised
                        COVID-19, longitudinal; Inflammation, Cardiometabolic, Cardiovascular II in
                        three separate runs (plate 1, plate 2, plates 3-7), each normalised alone

Values: the raw NPX exported by Olink NPX Manager ("raw data/Cohort *.xlsx", which the repository
describes as "raw data as obtained from Olink"), not the authors' bridged files (data/*/integration_*,
data/Radboud/merged_data) or data/data.RData. One matrix column per measured sample per run. Bridge
samples (samples of one cohort measured in the run of another; the authors used them to bridge-normalise
to MHH) are kept in the run where they were measured, flagged ``bridge_sample``, with the labels of
their home cohort. No bridge normalisation, no normalisation, no imputation, no cross-plate or cross-run
correction. Vendor QC only: values of a sample on a panel whose Olink "QC Warning" is not "Pass" are set
missing. Values below LOD are kept; LOD and the missing data frequency rows are not used. No assay
filtering (the authors' 20% / 80% detection filters are not applied).

Excluded wells: Nijmegen "control2-1" (a control well on Cardiovascular II plate 1), and every sample of
a Nijmegen patient absent from data/Radboud/patient_information_Aline.xlsx, the file that the authors'
load_data.R uses to keep "only the samples we are allowed to use" (16 patients, also when measured as
bridge samples in the Hannover or Breda run).

Labels: no participant has a long COVID / recovered label. Hannover convalescents (authors:
"post-COVID-19 symptom-free individuals"; Discussion: "largely symptom-free"; no per-person symptom
data deposited) are lc_status unknown, sample_phase post_acute. Hannover healthy donors are
never_infected (seronegative, recruited April to October 2020: never_infected_contemporary). All
hospitalised samples are acute with lc_status unknown (no later outcome).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Zoodsma2022"
REPO = "Postcovid-targeted-proteomics"
F_BREDA = "raw data/Cohort 1- Breda.xlsx"
F_NIJM = "raw data/Cohort 2- Nijmegen.xlsx"
F_HANN = "raw data/Cohort 3- Hannover.xlsx"
F_MHH_ANNOT = "data/MHH/patient_information_complete.xlsx"
F_BREDA_ANNOT = "data/Breda/patient_info_breda.xlsx"
F_BREDA_MANIFEST = "data/Breda/Olink-Sample Manifest Plate C4 - C10_foranalysis.xlsx"
F_BREDA_ARCHIVE = "data/Breda/2020 05 02 basis EDTA plasma buisarchief proteomics anoniem (1).xlsx"
F_RUMC_ALLOWED = "data/Radboud/patient_information_Aline.xlsx"

PANEL_NAME = {
    "INFLAMMATION": "Inflammation",
    "CARDIOMETABOLIC": "Cardiometabolic",
    "CARDIOVASCULAR II": "Cardiovascular II",
    "NEUROLOGY": "Neurology",
}
# the authors' PCA outliers (R/1 - QC.R), removed from their analyses; kept here and flagged
PCA_OUTLIERS = {"TK_41", "RUMC_0112 W1T1"}
OVERLAP = {
    "Hannover": "Bonifacius A et al. 2021 Immunity doi:10.1016/j.immuni.2021.01.008 (cohort partially described)",
    "Nijmegen": "Janssen NA et al. 2021 J Infect Dis doi:10.1093/infdis/jiab065 (same cohort and Olink data)",
    "Breda": "",
}
COUNTRY = {"Hannover": "Germany", "Nijmegen": "Netherlands", "Breda": "Netherlands"}


# ---------------------------------------------------------------------------------- raw NPX
def clean_id(s) -> str:
    """Whitespace collapsed; 'RUMC_0011W2T2' -> 'RUMC_0011 W2T2'; 'PC - HC04' -> 'PC-HC04'."""
    s = re.sub(r"\s+", " ", str(s)).strip()
    s = re.sub(r"^(RUMC_\d{4})(W\d)", r"\1 \2", s)
    s = re.sub(r"^PC\s*-\s*", "PC-", s)
    return s


def read_npx_sheet(path: Path, sheet: str, run: str) -> tuple:
    """One NPX Manager wide sheet -> (long table, run information)."""
    d = pd.read_excel(path, sheet_name=sheet, header=None)
    c0 = d[0].astype(str).str.strip()
    hdr = d.index[c0 == "Assay"][0]
    lod = d.index[c0 == "LOD"][0]
    head = {str(k).strip(): v for k, v in d.iloc[:hdr, :2].dropna(how="all").values.tolist()}
    panel_raw = str(head.get("Panel"))
    m = re.match(r"^(?:Olink Target 96 |Olink )?(.+?)\s*\(v\.(\d+)\)$", panel_raw)
    panel = PANEL_NAME[m.group(1).upper()]
    version = m.group(2)
    npx_manager = next((str(v) for k, v in head.items() if "NPX Manager" in str(v)), None) or next(
        (k for k in head if "NPX Manager" in k), ""
    )
    tail = {str(k).strip(): v for k, v in d.iloc[lod:, :2].dropna(how="all").values.tolist()}
    normalisation = head.get("Normalization Method", tail.get("Normalization"))
    assay, uni, oid = d.iloc[hdr].tolist(), d.iloc[hdr + 1].tolist(), d.iloc[hdr + 2].tolist()
    body = d.iloc[hdr + 3 : lod]
    body = body[body[0].notna()]
    cols = [i for i, o in enumerate(oid) if isinstance(o, str) and o.startswith("OID")]
    qi = assay.index("QC Warning")
    pi = assay.index("Plate ID") if "Plate ID" in assay else None
    ids = body[0].map(clean_id)
    # a sample id that occurs twice in one run gets a well suffix
    dup = ids.duplicated(keep=False)
    ids = ids.where(~dup, ids + "_dup" + (ids.groupby(ids).cumcount() + 1).astype(str))
    wide = body[cols].apply(pd.to_numeric, errors="coerce")
    wide.columns = [oid[i] for i in cols]
    wide.index = ids.values
    L = wide.stack(future_stack=True).rename("NPX").reset_index()
    L.columns = ["SampleID", "OlinkID", "NPX"]
    L["UniProt"] = L["OlinkID"].map({oid[i]: uni[i] for i in cols})
    L["Assay"] = L["OlinkID"].map({oid[i]: assay[i] for i in cols})
    qc = pd.Series(body[qi].astype(str).str.strip().values, index=ids.values)
    plate = pd.Series(body[pi].astype(str).values if pi is not None else [""] * len(body), index=ids.values)
    L["qc_warning"] = L["SampleID"].map(qc)
    L["vendor_qc"] = np.where(L["qc_warning"] == "Pass", "PASS", "FAIL")
    L["panel"] = panel
    L["run"] = run
    L["plate"] = L["SampleID"].map(plate)
    info = {
        "panel": panel,
        "panel_version": version,
        "npx_manager": npx_manager.replace("Olink NPX Manager ", ""),
        "normalisation": "" if normalisation is None or pd.isna(normalisation) else str(normalisation),
        "project": next((k for k, v in head.items() if "NPX Manager" in str(v) and k != "nan"), ""),
    }
    return L, info


def read_cohort_npx(path: Path, sheets: dict) -> tuple:
    """sheets: {sheet name: run label}. Returns (long table, per run x panel info table)."""
    parts, infos = [], []
    for sheet, run in sheets.items():
        L, info = read_npx_sheet(path, sheet, run)
        parts.append(L)
        infos.append({"run": run, "sheet": sheet, **info})
    return pd.concat(parts, ignore_index=True), pd.DataFrame(infos)


def run_samples(L: pd.DataFrame) -> pd.DataFrame:
    """One row per (sample, run): plates, panels, QC warnings."""
    g = L.groupby(["SampleID", "run"])
    out = pd.DataFrame(
        {
            "panels": g["panel"].apply(lambda s: ";".join(sorted(s.unique()))),
            "plates": g["plate"].apply(lambda s: ";".join(sorted(x for x in s.unique() if x))),
            "qc_warning_panels": g.apply(
                lambda t: ";".join(sorted(t.loc[t["vendor_qc"] != "PASS", "panel"].unique())), include_groups=False
            ),
            "n_values_qc_masked": g["vendor_qc"].apply(lambda s: int((s != "PASS").sum())),
        }
    ).reset_index()
    return out


# ---------------------------------------------------------------------------------- annotation
def sex_of(v) -> str:
    s = str(v).strip().lower()
    return {"male": "male", "female": "female", "man": "male", "vrouw": "female"}.get(s, "unknown")


class Labels:
    """Per-sample labels from the deposited annotation of each cohort."""

    def __init__(self, raw: Path):
        r = raw / REPO
        # Hannover (MHH)
        a = pd.read_excel(r / F_MHH_ANNOT)
        a["SampleID"] = a["SampleID"].map(clean_id)
        self.mhh = a.drop_duplicates("SampleID").set_index("SampleID")
        # Nijmegen (Radboud): sample annotation and the authors' list of usable patients
        n = pd.read_excel(r / F_NIJM, sheet_name="Annot")
        n["full_id"] = n["full_id"].map(clean_id)
        self.rumc = n.set_index("full_id")
        al = pd.read_excel(r / F_RUMC_ALLOWED)
        self.rumc_allowed = al.set_index("Record_ID")
        # Breda: Olink well code -> sample, patient info, admission dates from the EDTA archive
        b = pd.read_excel(r / F_BREDA_ANNOT, dtype={"patient_nr": str, "date": str})
        self.breda_info = b.set_index("Olink.ID")
        man = pd.read_excel(r / F_BREDA_MANIFEST, sheet_name=None)
        mm = pd.concat([t.iloc[:, [1, 2, 3]].set_axis(["well", "sid", "olink"], axis=1).assign(plate=k) for k, t in man.items()])
        mm["olink"] = mm["olink"].astype(str).str.strip()
        self.breda_manifest = mm.set_index("olink")
        e = pd.read_excel(r / F_BREDA_ARCHIVE, header=1, dtype={"studienummer KCHL lijst": str})
        e = e[e["studienummer KCHL lijst"].notna()].copy()
        e["nr"] = e["studienummer KCHL lijst"].str.strip().str.zfill(3)
        self.breda_adm = e.groupby("nr")["Datum opname"].first()
        self.breda_bmi = e.groupby("nr")["BMI"].first().map(lambda v: pd.to_numeric(str(v).replace(",", "."), errors="coerce"))

    def rumc_is_allowed(self, sid: str) -> bool:
        return sid.split(" ")[0] in self.rumc_allowed.index

    # each returns a dict of metadata for one sample id of that cohort
    def hannover(self, sid: str) -> dict:
        if sid not in self.mhh.index:
            return {"subject_id": sid, "lc_status": "unknown", "case_definition": "not_applicable", "control_type": "not_applicable",
                    "sample_phase": "unknown", "home_cohort": "Hannover", "notes": ["not in patient_information_complete.xlsx"]}
        m = self.mhh.loc[sid]
        cond = str(m["condition"])
        row = {
            "subject_id": sid,
            "home_cohort": "Hannover",
            "age": m["age"],
            "sex": sex_of(m["gender"]),
            "sampling_date": pd.to_datetime(m["date"]).date().isoformat() if pd.notna(m["date"]) else None,
            "authors_condition": cond,
            "days_since_symptom_resolution": m["time.convalescence"],
            "who_score_authors": m["WHO.score"],
            "disease_severity_score_authors": m["disease.severity.score"],
            "antibody_S1_ratio": m["antibody.S1"],
            "antibody_NCP_ratio": m["antibody.ncp"],
            "timepoint": "single_visit",
            "notes": [],
        }
        if cond == "postcovid":
            row.update(lc_status="unknown", case_definition="not_applicable", control_type="not_applicable",
                       sample_phase="post_acute", variant_era="pre_alpha")
            row["notes"].append("post-COVID-19 convalescent (PCR proven, mild); authors call the group symptom-free / largely symptom-free; "
                                "no per-person symptom data, LC status unknown")
        elif cond == "healthy":
            row.update(lc_status="never_infected", case_definition="not_applicable", control_type="never_infected_contemporary",
                       sample_phase="not_applicable")
            row["notes"].append("healthy blood donor, SARS-CoV-2 seronegative, recruited 2020")
        else:
            raise ValueError(f"unexpected MHH condition {cond} for {sid}")
        return row

    def nijmegen(self, sid: str) -> dict:
        pat, _, tp = sid.partition(" ")
        tp = re.sub(r"_dup\d+$", "", tp)
        row = {"subject_id": pat, "home_cohort": "Nijmegen", "timepoint": tp, "lc_status": "unknown",
               "case_definition": "not_applicable", "control_type": "not_applicable", "sample_phase": "acute",
               "hospitalised": True, "variant_era": "pre_alpha", "notes": []}
        base = re.sub(r"_dup\d+$", "", sid)
        if base in self.rumc.index:
            m = self.rumc.loc[base]
            if isinstance(m, pd.DataFrame):
                m = m.iloc[0]
            src = "Annot"
        else:
            m = self.rumc_allowed.loc[pat]
            src = "patient_information_Aline"
            row["notes"].append("sample not in the Nijmegen Annot sheet; patient-level data from patient_information_Aline.xlsx")
        cond = str(m["condition"])
        pcr = m.get("PCR_proven_COVID-19")
        row.update(
            age=m["age"], sex=sex_of(m["gender"]), bmi=m["BMI"], authors_condition=cond,
            acute_severity="icu_or_ventilated" if cond == "ICU" else "unknown",
            pcr_proven=("no" if str(pcr).strip().lower() == "no" else "") if src == "Annot" else "",
        )
        if str(pcr).strip().lower() == "no":
            row["notes"].append("COVID-19 not PCR proven (clinical features and CT per paper)")
        row["notes"].append(f"hospitalised COVID-19 ({cond}); week/timepoint label {tp}; days since admission not deposited")
        return row

    def breda(self, sid: str) -> dict:
        row = {"lc_status": "unknown", "case_definition": "not_applicable", "control_type": "not_applicable",
               "home_cohort": "Breda", "notes": []}
        if sid in self.breda_info.index:
            m = self.breda_info.loc[sid]
            nr = str(m["patient_nr"]).zfill(3)
            # two dates carry typos ("01-040-2020", "04-04-,2020"): day, first two digits of month, year
            parts = [p for p in re.sub(r"[^0-9-]", "", str(m["date"])).split("-") if p]
            date = pd.to_datetime(f"{parts[0]}-{parts[1][:2]}-{parts[2]}", format="%d-%m-%Y", errors="coerce") if len(parts) == 3 else pd.NaT
            adm = self.breda_adm.get(nr)
            days = (date - adm).days if pd.notna(date) and pd.notna(adm) else np.nan
            icu = str(m["ICU"]).strip().lower()
            row.update(
                subject_id=str(m["breda_nr"]).strip(), breda_sample=re.sub(r"\s+", " ", str(m["Sample.ID"])).strip(),
                timepoint=str(m["timepoint"]).strip(), sample_phase="acute", hospitalised=True, variant_era="pre_alpha",
                age=m["age"], sex=sex_of(m["gender"]), bmi=self.breda_bmi.get(nr),
                authors_condition={"yes": "ICU", "no": "non-ICU"}.get(icu, ""),
                acute_severity={"yes": "icu_or_ventilated"}.get(icu, "unknown"),
                sampling_date=date.date().isoformat() if pd.notna(date) else None,
                days_since_infection=days if pd.notna(days) else np.nan,
                time_anchor="admission" if pd.notna(days) else None,
            )
            row["notes"].append("hospitalised COVID-19" + (f" ({row['authors_condition']}, patient ever on ICU)" if icu == "yes" else
                                (" (non-ICU)" if icu == "no" else "; ICU status, age and sex not deposited")))
            if pd.isna(days):
                row["notes"].append("admission or sampling date not available")
            else:
                row["notes"].append("days_since_infection = sampling date minus admission date (EDTA archive)")
        else:
            man = self.breda_manifest.sid.get(sid)
            man = None if man is None or pd.isna(man) else re.sub(r"\s+", " ", str(man)).strip()
            row.update(subject_id=f"Breda_unannotated_{sid}" if man is None else man.split(" ")[0].replace("Breda-", "Breda_"),
                       breda_sample=man or "", sample_phase="unknown", timepoint="")
            row["notes"].append(f"well {sid} not in the authors' patient_info_breda.xlsx (manifest id: {man or 'none'}); not analysed by the authors")
        return row


# ---------------------------------------------------------------------------------- per cohort
def label_sample(lab: Labels, sid: str, run_cohort: str) -> dict:
    """Labels of one measured sample of the run of ``run_cohort``; bridge samples carry their home labels."""
    base = re.sub(r"_dup\d+$", "", sid)
    if base.startswith("RUMC_"):
        row = lab.nijmegen(sid)
    elif re.match(r"^(Co|TK)_?\d", base):
        row = lab.hannover(base)
    elif re.match(r"^[A-H]-\d+$", base):
        row = lab.breda(base)
    elif base.startswith("PC-"):
        row = {"subject_id": base, "lc_status": "unknown", "case_definition": "not_applicable", "control_type": "not_applicable",
               "sample_phase": "unknown", "home_cohort": "other_project", "timepoint": "",
               "notes": ["sample of another project measured on the Breda plates (PC-PAT/PC-HC); not part of this study, no annotation"]}
    else:
        raise ValueError(f"unexpected sample id {sid}")
    row["bridge_sample"] = row["home_cohort"] != run_cohort and row["home_cohort"] != "other_project"
    if row["bridge_sample"]:
        row["notes"].insert(0, f"bridge sample of the {row['home_cohort']} cohort measured in the {run_cohort} run "
                               "(authors' bridging to MHH); not bridge normalised here")
    row["authors_pca_outlier"] = base in PCA_OUTLIERS
    if row["authors_pca_outlier"]:
        row["notes"].append("authors' PCA outlier (removed from their analyses)")
    return row


def build_one(cohort: str, L: pd.DataFrame, runinfo: pd.DataFrame, lab: Labels, drop_ids=()) -> dict:
    n_dropped = len(set(drop_ids))
    L = L[~L["SampleID"].isin(set(drop_ids))].copy()
    multi_run = L["run"].nunique() > 1
    L["sample_id"] = (L["SampleID"].str.replace(" ", "_") + ("_" + L["run"] if multi_run else "")).values
    rs = run_samples(L)
    rs["sample_id"] = rs["SampleID"].str.replace(" ", "_") + ("_" + rs["run"] if multi_run else "")
    rs = rs.set_index("sample_id")

    adapted = to_feature_space("olink", L, sample_col="sample_id", assay_col="OlinkID", uniprot_col="UniProt", value_col="NPX",
                               qc_col="vendor_qc", qc_fail=("FAIL",))
    # samples whose values are all missing in the export (the adapter's pivot drops them)
    empty = sorted(set(rs.index) - set(adapted.X.index))
    order = [x for x in rs.sort_values(["run", "SampleID"]).index if x in adapted.X.index]
    adapted.X = adapted.X.loc[order]
    rs = rs.loc[order]

    rows = []
    for sid, r in rs.iterrows():
        row = label_sample(lab, r["SampleID"], cohort)
        info = runinfo[runinfo["run"] == r["run"]]
        panels = r["panels"].split(";")
        vers = {p: v for p, v in zip(info["panel"], info["panel_version"])}
        row["platform"] = "Olink Target 96 " + ", ".join(f"{p} (v.{vers[p]})" for p in panels) + "; NPX"
        row["sample_matrix"] = "plasma"
        row["collection_country"] = COUNTRY.get(row["home_cohort"], "")
        row["overlaps_with_study"] = OVERLAP.get(row["home_cohort"], "")
        row["original_sample_id"] = r["SampleID"]
        row["olink_run"] = r["run"]
        row["olink_plates"] = r["plates"]
        row["panels_measured"] = r["panels"]
        row["npx_manager"] = ";".join(sorted(set(info["npx_manager"])))
        row["normalisation_method"] = ";".join(sorted(set(x for x in info["normalisation"] if x))) or "not stated"
        row["qc_warning_panels"] = r["qc_warning_panels"]
        row["n_values_qc_masked"] = r["n_values_qc_masked"]
        if r["qc_warning_panels"]:
            row["notes"].append(f"Olink QC Warning on {r['qc_warning_panels']}: those values set missing")
        rows.append(row)
    meta = pd.DataFrame(rows, index=rs.index)
    meta.index.name = "sample_id"
    meta["notes"] = meta["notes"].map(lambda n: "; ".join(n))
    for c in ["time_anchor", "acute_severity", "variant_era", "hospitalised", "sex"]:
        if c in meta.columns:
            meta[c] = meta[c].where(meta[c].notna() & (meta[c].astype(str) != ""), None)

    extra = {
        "n_runs": int(L["run"].nunique()),
        "n_excluded_samples": n_dropped,
        "samples_all_missing_dropped": ";".join(empty),
        "n_bridge_samples": int(meta["bridge_sample"].sum()),
        "n_qc_masked_values": adapted.notes.get("qc_masked_values"),
        "assays_without_uniprot": adapted.notes.get("assays_without_uniprot"),
        "n_assays": int(L["OlinkID"].nunique()),
        "n_postcovid_unknown": int(((meta["lc_status"] == "unknown") & (meta["sample_phase"] == "post_acute")).sum()),
        "n_with_days": int(pd.to_numeric(meta.get("days_since_infection"), errors="coerce").notna().sum()) if "days_since_infection" in meta else 0,
        "n_timepoint": ",".join(f"{k}={v}" for k, v in meta["timepoint"].fillna("").value_counts().sort_index().items()),
        "values": "raw NPX (NPX Manager export), vendor QC Warning masked, not bridged",
    }
    return write_cohort(f"{STUDY}_{cohort}", adapted, meta, extra)


def build() -> list:
    raw = raw_dir(STUDY)
    r = raw / REPO
    lab = Labels(raw)

    # Nijmegen samples of patients the authors were not allowed to use, and a control well
    def excluded(ids) -> list:
        out = []
        for s in ids:
            if s.startswith("control"):
                out.append(s)
            elif s.startswith("RUMC_") and not lab.rumc_is_allowed(re.sub(r"_dup\d+$", "", s)):
                out.append(s)
        return out

    out = []
    L, info = read_cohort_npx(r / F_HANN, {"INF": "MHH", "CM": "MHH", "CVD II": "MHH", "NEU": "MHH"})
    out.append(build_one("Hannover", L, info, lab, excluded(L["SampleID"].unique())))

    L, info = read_cohort_npx(r / F_BREDA, {"INF": "Breda"})
    out.append(build_one("Breda", L, info, lab, excluded(L["SampleID"].unique())))

    sheets = {
        "INF plat 1": "run1", "CM plate 1": "run1", "CVD plate 1": "run1",
        "INF plate 2": "run2", "CM plate 2": "run2", "CVD plate 2": "run2",
        "INF plate 3-7": "run3", "CM plates 3-7": "run3", "CVD plates 3-7": "run3",
    }
    L, info = read_cohort_npx(r / F_NIJM, sheets)
    out.append(build_one("Nijmegen", L, info, lab, excluded(L["SampleID"].unique())))
    return out


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
        print()
