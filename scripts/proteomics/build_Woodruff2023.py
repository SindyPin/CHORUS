"""Build the Woodruff2023 Olink cohort (Woodruff MC, Bonham KS, et al. 2023 Nat Commun 14:4201, "Chronic inflammation,
neutrophil activity, and autoreactivity splits long COVID", doi:10.1038/s41467-023-40012-7, PMID 37452024,
PMC10349085; analysis code release Zenodo 10.5281/zenodo.8092298, Covid-v0.4.0, CC BY 4.0).

Raw files (data/proteomics/Woodruff2023/raw):
  Covid-v0.4.0.tar.gz   the release (md5 checked); only data/complete.csv is taken out of it
  complete.csv          long table, one row per sample and Olink assay (Olink Explore 3072, NPX) with the authors'
                        metadata: sample, group_hr, olinkID, uniprot, gene, NPX, group, sex, race, disease_phase,
                        initial_disease_severity, dpso (days post symptom onset)

Labels: ``group_hr`` per sample. PASC = LC (case_definition study_specific: Emory PASC clinics, documented infection or
the CDC clinical case definition, new or worsening symptoms more than 14 days after onset; paper Methods). CR =
COVID_recovered (donors with uncomplicated recoveries, sampled at similar intervals; the paper lists 26, the table has
27: sample 1433-COV-3140, which the deposit's data/README.md notes was listed in both PASC and CR in the Olink table and
as CR in the metadata). ICU Pre-Dex and ICU Post-Dex = acute ICU samples before or after the start of dexamethasone,
lc_status unknown (no later outcome), sample_phase acute, so they are outside the LC contrast. SLE samples are dropped.

Values: NPX as deposited (no QC column in the table; nothing masked). Each sample is one subject (no id repeats).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Woodruff2023"
PLATFORM = "Olink Explore 3072 (NPX)"
STATUS = {"PASC": "LC", "CR": "COVID_recovered", "ICU Pre-Dex": "unknown", "ICU Post-Dex": "unknown"}
SEVERITY = {"out": "ambulatory", "hosp": "unknown", "icu": "icu_or_ventilated"}  # hosp: oxygen not recorded
SEX = {"f": "female", "m": "male"}


def build() -> list:
    L = pd.read_csv(raw_dir(STUDY) / "complete.csv", low_memory=False)
    n_sle = L.loc[L["group_hr"] == "SLE", "sample"].nunique()
    L = L[L["group_hr"] != "SLE"].copy()
    unknown = set(L["group_hr"]) - set(STATUS)
    if unknown:
        raise ValueError(f"unexpected group_hr values {sorted(unknown)}")
    adapted = to_feature_space("olink", L, sample_col="sample", assay_col="olinkID", uniprot_col="uniprot", value_col="NPX")
    adapted.X = adapted.X.loc[sorted(adapted.X.index)]
    s = L.drop_duplicates("sample").set_index("sample").loc[adapted.X.index]
    g = s["group_hr"]
    acute = g.str.startswith("ICU")
    meta = pd.DataFrame({
        "subject_id": STUDY + "_" + s.index.to_series(),
        "lc_status": g.map(STATUS),
        "case_definition": np.where(g == "PASC", "study_specific", "not_applicable"),
        "control_type": np.where(g == "CR", "COVID_recovered", "not_applicable"),
        "sample_phase": np.where(acute, "acute", "post_acute"),
        "sample_matrix": "plasma",
        "days_since_infection": s["dpso"].astype(float),
        "time_anchor": "symptom_onset",
        "acute_severity": s["initial_disease_severity"].map(SEVERITY).fillna("unknown"),
        "hospitalised": s["initial_disease_severity"].map({"out": False, "hosp": True, "icu": True}),
        "sex": s["sex"].map(SEX).fillna("unknown"),
        "platform": PLATFORM,
        "timepoint": "single_visit",
        "collection_country": "United States",
        "group_hr": g,
        "dexamethasone": g.map({"ICU Pre-Dex": "before", "ICU Post-Dex": "after"}),
        "race": s["race"],
        "notes": g.map({
            "PASC": "PASC clinic patient: documented infection or CDC clinical case definition, new or worsening symptoms more than 14 days after onset",
            "CR": "uncomplicated recovery from COVID-19, sampled at a similar interval after onset",
            "ICU Pre-Dex": "acute ICU sample before dexamethasone; later outcome unknown",
            "ICU Post-Dex": "acute ICU sample after dexamethasone; later outcome unknown",
        }),
    }, index=adapted.X.index)
    meta.index.name = "sample_id"
    extra = {
        "n_assays": int(L["olinkID"].nunique()),
        "n_SLE_dropped": int(n_sle),
        **{f"n_{k.replace(' ', '_')}": int((g == k).sum()) for k in STATUS},
        "n_LC_dpso_below_28": int(((g == "PASC") & (s["dpso"] < 28)).sum()),
        "n_CR_dpso_below_28": int(((g == "CR") & (s["dpso"] < 28)).sum()),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for r in build():
        for k, v in r.items():
            print(f"{k}: {v}")
