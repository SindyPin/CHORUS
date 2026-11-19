"""Build the Bansal2025 Olink cohort (Bansal A, Olechnowicz SWZ, et al. 2026 Commun Med, "Divergent inflammatory and
neurology-related protein levels in long COVID following primary and breakthrough SARS-CoV-2 infections",
doi:10.1038/s43856-026-01541-6, PMID 41975235, PMC13442878; Zenodo 10.5281/zenodo.15237087, version 15237088,
CC BY 4.0). COVID PROFILE cohort, Victoria, Australia.

Raw file (data/proteomics/Bansal2025/raw): baseline57.xlsx, the cross-sectional baseline set: one row per participant
and Olink Target 96 assay (Inflammation and Neurology panels, 182 assays, NPX) with the participant's metadata
(study_group, Long_COVID status, gender_at_birth, age, comorbidity, disease_status_subgroup, vaccination, antibody
levels, Time). The longitudinal workbook (Longitudinal_paired34.xlsx) is not used.

Labels: ``Long_COVID status`` per participant. Long COVID = LC (WHO 2021: symptoms persisting or new 3 months after
infection, lasting at least 2 months); Recovered = COVID_recovered (symptoms shorter than 90 days); Healthy =
never_infected (unvaccinated, SARS-CoV-2 naive, never_infected_contemporary). COVID-19 samples were drawn 6 to 9 months
after the first infection and before COVID-19 vaccination (``Time`` Pre-vax for every row).

Vendor QC only: values whose QC_Warning is not Pass are set missing. The ``Normalization`` column is not used.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Bansal2025"
PLATFORM = "Olink Target 96 Inflammation and Neurology (NPX)"
STATUS = {"Long COVID": "LC", "Recovered": "COVID_recovered", "Healthy": "never_infected"}
SEVERITY = {"Moderate/Mild": "ambulatory", "Asymptomatic": "ambulatory", "Severe": "unknown"}


def build() -> list:
    L = pd.read_excel(raw_dir(STUDY) / "baseline57.xlsx")
    L["sample"] = STUDY + "_" + L["Study Participant"].astype(str)
    per = L.groupby("sample")[["Long_COVID status", "study_group", "gender_at_birth", "age", "Time"]].nunique()
    if (per > 1).any().any():
        raise ValueError("participant with more than one metadata value")
    unknown = set(L["Long_COVID status"]) - set(STATUS)
    if unknown:
        raise ValueError(f"unexpected Long_COVID status values {sorted(unknown)}")
    L["vendor_qc"] = np.where(L["QC_Warning"].astype(str) == "Pass", "PASS", "FAIL")
    adapted = to_feature_space("olink", L, sample_col="sample", assay_col="OlinkID", uniprot_col="UniProt", value_col="NPX", qc_col="vendor_qc", qc_fail=("FAIL",))
    adapted.X = adapted.X.loc[sorted(adapted.X.index, key=lambda x: int(x.split("_")[1]))]
    s = L.drop_duplicates("sample").set_index("sample").loc[adapted.X.index]
    st = s["Long_COVID status"].map(STATUS)
    meta = pd.DataFrame({
        "subject_id": s.index,
        "lc_status": st,
        "case_definition": np.where(st == "LC", "WHO_2021", "not_applicable"),
        "control_type": st.map({"COVID_recovered": "COVID_recovered", "never_infected": "never_infected_contemporary"}).fillna("not_applicable"),
        "sample_phase": np.where(st == "never_infected", "not_applicable", "post_acute"),
        "sample_matrix": "plasma",
        "time_anchor": np.where(st == "never_infected", None, "positive_test"),
        "acute_severity": s["disease_status_subgroup"].map(SEVERITY),
        "age": s["age"],
        "sex": s["gender_at_birth"].str.lower(),
        "platform": PLATFORM,
        "timepoint": np.where(st == "never_infected", "enrolment", "6_to_9_months"),
        "collection_country": "Australia",
        "study_group": s["study_group"],
        "disease_status_subgroup": s["disease_status_subgroup"],
        "comorbidity": s["Active Comorbidity Present/Absent "],
        "vaccinated_before_sampling": False,
        "olink_plate": s["PlateID"].astype(str),
        "notes": st.map({
            "LC": "Long COVID (WHO definition); plasma 6 to 9 months after the first infection, before vaccination",
            "COVID_recovered": "symptoms shorter than 90 days; plasma 6 to 9 months after the first infection, before vaccination",
            "never_infected": "unvaccinated SARS-CoV-2 naive healthy participant, plasma at enrolment",
        }),
    }, index=adapted.X.index)
    meta.index.name = "sample_id"
    plates = L.groupby("sample")["PlateID"].apply(lambda p: ";".join(sorted(p.astype(str).unique())))
    meta["olink_plate"] = plates.reindex(meta.index).values
    extra = {
        "n_assays": int(L["OlinkID"].nunique()),
        "qc_masked_values": adapted.notes.get("qc_masked_values"),
        **{f"n_{k}": int((st == k).sum()) for k in ("LC", "COVID_recovered", "never_infected")},
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for r in build():
        for k, v in r.items():
            print(f"{k}: {v}")
