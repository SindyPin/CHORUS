"""Build the Fineschi2026 Olink cohort (Fineschi S, Klar J, Schuster J, Bergquist J, Dahl N 2026 Front Immunol
17:1775044, "Plasma proteomic profile reveals persistent immune activation in post-acute sequelae of SARS-CoV-2
infection", doi:10.3389/fimmu.2026.1775044, PMID 41808838, PMC12968220; SciLifeLab figshare
10.17044/scilifelab.31316779, CC BY 4.0).

Raw files (data/proteomics/Fineschi2026/raw):
  XL-4180_raw_data.csv     Olink Explore 384 Inflammation long NPX export (";"), two plates
  README.txt, MANIFEST.txt, XL-4180_Analysis_Report_2025-06-04.pdf

Labels: the deposit has no separate sample sheet. The group is the prefix of SampleID: PA-nn (patient) and KO-nn
(kontroll); the README gives 92 patients and 73 controls, which are the counts of PA and KO ids in the export. The
same prefix carries the group in the authors' PBMC RNA-seq deposit (10.17044/scilifelab.28832492, SamplesPC.txt),
where every PA id has Condition PAT and every KO id CTL.

Definitions (paper Methods, Participants): patients met the WHO post COVID-19 condition criteria; controls are
age and sex matched people infected in the same period who recovered without persistent symptoms. All were
infected between February 2020 and June 2022 without hospitalisation; blood was drawn in April-May 2023 or
November-December 2024 (mean 34 months after infection). Uppsala, Sweden.

Vendor QC only: values whose QC_Warning or Assay_Warning is not PASS are set missing. Values below LOD are kept.
Olink CONTROL wells are dropped. No normalisation or imputation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Fineschi2026"
PLATFORM = "Olink Explore 384 Inflammation (Explore 3.11.1)"
PBMC_SAMPLES = ROOT / "data" / "counts" / "Fineschi2025_check" / "SamplesPC.txt"  # label file of 10.17044/scilifelab.28832492


def build() -> list:
    raw = raw_dir(STUDY)
    L = pd.read_csv(raw / "XL-4180_raw_data.csv", sep=";", low_memory=False)
    L = L[L["Sample_Type"] == "SAMPLE"].copy()
    ok = (L["QC_Warning"].astype(str) == "PASS") & (L["Assay_Warning"].astype(str) == "PASS")
    L["vendor_qc"] = np.where(ok, "PASS", "FAIL")
    adapted = to_feature_space("olink", L, sample_col="SampleID", assay_col="OlinkID", uniprot_col="UniProt", value_col="NPX", qc_col="vendor_qc", qc_fail=("FAIL",))
    adapted.X = adapted.X.loc[sorted(adapted.X.index)]

    g = L.groupby("SampleID")
    pbmc = pd.read_csv(PBMC_SAMPLES, sep="\t").dropna(subset=["Sample"]) if PBMC_SAMPLES.exists() else None
    in_pbmc = set(pbmc["Sample"].str.replace(".", "-", regex=False)) if pbmc is not None else set()
    rows = []
    for s in adapted.X.index:
        prefix = s.split("-")[0]
        if prefix not in ("PA", "KO"):
            raise ValueError(f"unexpected sample id {s}")
        lc = prefix == "PA"
        rows.append({
            "subject_id": f"{STUDY}_{s}",
            "lc_status": "LC" if lc else "COVID_recovered",
            "case_definition": "WHO_2021" if lc else "not_applicable",
            "control_type": "not_applicable" if lc else "COVID_recovered",
            "sample_phase": "post_acute",
            "sample_matrix": "plasma",
            "platform": PLATFORM,
            "timepoint": "single_visit",
            "time_anchor": "unknown",
            "acute_severity": "ambulatory",
            "hospitalised": False,
            "collection_country": "Sweden",
            "overlaps_with_study": "Fineschi2025_PBMC (scilifelab.28832492)" if s in in_pbmc else None,
            "olink_plate": ";".join(sorted(g.get_group(s)["PlateID"].unique())),
            "n_values_qc_masked": int((g.get_group(s)["vendor_qc"] != "PASS").sum()),
            "notes": ("group from the SampleID prefix (PA patient, KO control; README 92 / 73); "
                      + ("WHO post COVID-19 condition" if lc else "infected in the same period, recovered without persistent symptoms")
                      + "; not hospitalised; blood April-May 2023 or November-December 2024"),
        })
    meta = pd.DataFrame(rows, index=adapted.X.index)
    meta.index.name = "sample_id"
    extra = {
        "qc_masked_values": adapted.notes.get("qc_masked_values"),
        "n_assays": int(L["OlinkID"].nunique()),
        "n_PA": int((meta["lc_status"] == "LC").sum()),
        "n_KO": int((meta["lc_status"] == "COVID_recovered").sum()),
        "n_in_pbmc_rnaseq": int(meta["overlaps_with_study"].notna().sum()),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
