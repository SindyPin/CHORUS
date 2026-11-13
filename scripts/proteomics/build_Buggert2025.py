"""Build the Buggert2025 Olink cohort (Gao Y, ..., Buggert M 2025 Nat Immunol,
doi:10.1038/s41590-025-02135-5, PMID 40307449; Zenodo 10.5281/zenodo.14772494, CC-BY-4.0).

Raw files (data/proteomics/Buggert2025/raw):
  plate1_olink.csv, plate2_olink.csv, plate3_olink.csv   Olink Explore long NPX exports (";")
  meta_UKsamples.csv   primary cohort (University Hospital Llandough, UK): CA = long COVID, CO = convalescent
  meta_SEsamples.csv   secondary cohort (Karolinska University Hospital, Sweden): KLIMP = long COVID

One matrix column per run (sample_id = <subject_id>_p<plate file>). Subjects measured on two plate
files (the authors' bridge samples) keep both runs; no bridge normalisation is applied here.
Vendor QC only: values whose QC_Warning or Assay_Warning is not PASS are set missing. Values
below LOD are kept. No normalisation, imputation or cross-plate correction.
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

STUDY = "Buggert2025"
PLATES = {1: "plate1_olink.csv", 2: "plate2_olink.csv", 3: "plate3_olink.csv"}
PLATE_PLATFORM = {
    1: "Olink Explore 3072 (Explore 3.0.5)",
    2: "Olink Explore 3072 (Explore 3.6.0)",
    3: "Olink Explore 3072 (Explore 3.10.0)",
}


def norm_subject(s: str) -> str:
    """CA14 / CA014 -> CA014; CO08 -> CO008; KLIMP ids unchanged."""
    s = str(s).strip()
    m = re.match(r"^(CA|CO)0*(\d+)$", s)
    return f"{m.group(1)}{int(m.group(2)):03d}" if m else s


def read_plates(raw: Path) -> pd.DataFrame:
    parts = []
    for p, fn in PLATES.items():
        d = pd.read_csv(raw / fn, sep=";", low_memory=False)
        d = d[d["Sample_Type"] == "SAMPLE"].copy()  # drop Olink CONTROL wells
        d["plate_file"] = p
        d["subject_id"] = d["SampleID"].map(norm_subject)
        d["run_id"] = d["subject_id"] + f"_p{p}"
        parts.append(d)
    L = pd.concat(parts, ignore_index=True)
    # the export labels the NT-proBNP assay (OID20125) "NTproBNP" instead of an accession, so the
    # adapter drops it; it is not mapped to P16860 because that key already holds the BNP assay
    # (OID20049) and the two assays would be averaged
    # combined vendor flag: anything not PASS in QC_Warning (sample) or Assay_Warning (assay)
    ok = (L["QC_Warning"].astype(str) == "PASS") & (L["Assay_Warning"].astype(str) == "PASS")
    L["vendor_qc"] = np.where(ok, "PASS", "FAIL")
    return L


def run_table(L: pd.DataFrame) -> pd.DataFrame:
    g = L.groupby("run_id")
    runs = pd.DataFrame(
        {
            "subject_id": g["subject_id"].first(),
            "original_sample_id": g["SampleID"].first(),
            "plate_file": g["plate_file"].first(),
            "olink_project": g["PlateID"].first().str.split("_").str[0],
            "olink_plate_ids": g["PlateID"].apply(lambda s: ";".join(sorted(s.unique()))),
            "explore_version": g["ExploreVersion"].first().astype(str),
            "panel_lots": g["Panel_Lot_Nr"].apply(lambda s: ";".join(sorted(s.astype(str).unique()))),
            "n_values_qc_masked": g["vendor_qc"].apply(lambda s: int((s != "PASS").sum())),
        }
    )
    runs["platform"] = runs["plate_file"].map(PLATE_PLATFORM)
    n_runs = runs.groupby("subject_id").size()
    runs["repeated_across_plates"] = runs["subject_id"].map(n_runs).gt(1)
    return runs


def build() -> list:
    raw = raw_dir(STUDY)
    L = read_plates(raw)
    runs = run_table(L)

    adapted = to_feature_space(
        "olink",
        L,
        sample_col="run_id",
        assay_col="OlinkID",
        uniprot_col="UniProt",
        value_col="NPX",
        qc_col="vendor_qc",
        qc_fail=("FAIL",),
    )
    order = sorted(adapted.X.index, key=lambda r: (runs.loc[r, "plate_file"], r))
    adapted.X = adapted.X.loc[order]
    runs = runs.loc[order]

    uk = pd.read_csv(raw / "meta_UKsamples.csv", encoding="utf-8-sig")
    uk["SampleID"] = uk["SampleID"].map(norm_subject)
    uk = uk.set_index("SampleID")
    se = pd.read_csv(raw / "meta_SEsamples.csv", encoding="utf-8-sig")
    se["SampleID"] = se["SampleID"].astype(str).str.strip()
    se = se.set_index("SampleID")

    sex_map = {"f": "female", "m": "male", "female": "female", "male": "male"}
    rows = []
    for rid, r in runs.iterrows():
        s = r["subject_id"]
        row = {"subject_id": s, "sample_matrix": "plasma", "platform": r["platform"], "timepoint": "single_visit", "time_anchor": None}
        notes = []
        if s.startswith("CA") or s.startswith("CO"):
            m = uk.loc[s]
            is_lc = s.startswith("CA")
            row.update(
                lc_status="LC" if is_lc else "COVID_recovered",
                case_definition="NICE_NG188" if is_lc else "not_applicable",
                control_type="not_applicable" if is_lc else "COVID_recovered",
                sample_phase="post_acute",
                collection_country="United Kingdom",
                acute_severity="ambulatory",
                hospitalised=False,
                age=m["Age"],
                sex=sex_map.get(str(m["Sex"]).strip().lower(), "unknown"),
                bmi=m["BMI"],
                cohort_arm="primary_UK_Llandough",
                borg_cr10=np.nan,
                dyspnoea12=m["Dyspnoea 12 (/36)"],
                nijmegen=m["Nijmegen (/64)"],
            )
            notes.append("primary cohort (UK); mild non-hospitalised acute COVID-19 per paper Methods (Donors)")
        elif s.startswith("KLIMP"):
            row.update(
                sample_phase="post_acute",
                collection_country="Sweden",
                cohort_arm="secondary_SE_Karolinska",
                dyspnoea12=np.nan,
                nijmegen=np.nan,
            )
            if s in se.index:
                m = se.loc[s]
                row.update(
                    lc_status="LC",
                    case_definition="study_specific",
                    control_type="not_applicable",
                    age=m["Age"],
                    sex=sex_map.get(str(m["Sex"]).strip().lower(), "unknown"),
                    bmi=m["BMI"],
                    borg_cr10=m["BORG-CR10"],
                )
                notes.append("secondary cohort (Sweden) long COVID; case definition not stated in the paper")
            else:
                row.update(lc_status="unknown", case_definition="not_applicable", control_type="not_applicable", borg_cr10=np.nan)
                notes.append("KLIMP sample on plates 2 and 3 but absent from meta_SEsamples.csv; group not deposited")
        else:
            raise ValueError(f"unexpected subject id {s}")
        if r["repeated_across_plates"]:
            notes.append("subject measured on two plate files (authors' bridge sample); runs kept separate, not bridge normalised")
        row["notes"] = "; ".join(notes)
        rows.append(row)

    meta = pd.DataFrame(rows, index=runs.index)
    meta.index.name = "sample_id"
    for c in ["original_sample_id", "plate_file", "olink_project", "olink_plate_ids", "explore_version", "panel_lots", "repeated_across_plates", "n_values_qc_masked"]:
        meta[c] = runs[c]

    extra = {
        "n_runs_repeated": int(runs["repeated_across_plates"].sum()),
        "n_subjects_repeated": int(runs.loc[runs["repeated_across_plates"], "subject_id"].nunique()),
        "qc_masked_values": adapted.notes.get("qc_masked_values"),
        "assays_without_uniprot": adapted.notes.get("assays_without_uniprot"),
        "n_assays": int(L["OlinkID"].nunique()),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
