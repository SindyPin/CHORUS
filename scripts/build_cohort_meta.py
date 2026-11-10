"""Build ClinicalSchema metadata tables for the transcriptomic cohorts that have LC labels.

Inputs are the per GSM GEO sample tables (data/meta/<GSE>_geo_samples.csv: title, source and
one ch_<characteristic> column per GEO characteristic) and, for GSE215865, the TACO clinical
table. Each output (data/meta/<study_id>_meta.csv) has the schema columns plus ``gsm`` and
``timepoint``, uses the matrix column ids as sample_id, and is validated with
validate_metadata(strict=True) before it is written. Subject ids are prefixed with the series
accession so that unrelated studies cannot collide in deduplicate_subjects.

  python scripts/build_cohort_meta.py --taco ../TACO/Data/411_DCEG_Paths_average.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from chorus.schema import ClinicalSchema, validate_metadata

COLUMNS = ClinicalSchema().columns + ["gsm", "timepoint"]


def finish(df: pd.DataFrame, study_id: str) -> pd.DataFrame:
    df = df.copy()
    df["study_id"] = study_id
    df["omics_layer"] = "bulk_rna"
    df = df.reindex(columns=COLUMNS)
    validate_metadata(df, ClinicalSchema(), strict=True)
    return df


def gse215865(geo: pd.DataFrame, taco: pd.DataFrame) -> pd.DataFrame:
    """Thompson 2023 (Mount Sinai): labels from the TACO clinical table joined on subject id.

    Long_COVID 1 with a COVID-19 positive subject is LC, 0 with positive is COVID_recovered and
    0 with negative is never_infected (GEO ``covid-19 positive``). A subject labelled 1 but COVID
    negative is contradictory and excluded (lc_status unknown, with a note); subjects missing from
    the TACO table are unknown.
    """
    lc = taco.drop_duplicates("Subject_ID").set_index("Subject_ID")["Long_COVID"]
    positive = geo.groupby("ch_subject_id")["ch_covid-19_positive"].agg(lambda v: (v == "Positive").any())
    subj = geo["ch_subject_id"]
    flag, pos = subj.map(lc), subj.map(positive)
    status = pd.Series("unknown", index=geo.index)
    status[(flag == 1) & pos] = "LC"
    status[(flag == 0) & pos] = "COVID_recovered"
    status[(flag == 0) & ~pos] = "never_infected"
    notes = "timepoint is days since first sample (" + geo["ch_sampling_time_point_label"] + "), not days since infection; Long_COVID from TACO (Synapse syn35874390)"
    notes[flag.isna()] = "no Long_COVID value in the TACO table"
    notes[(flag == 1) & ~pos] = "excluded: Long_COVID 1 in TACO but COVID-19 negative in GEO; check against Synapse syn35874390"
    df = pd.DataFrame(
        {
            "sample_id": geo["ch_blood_sample_id"] + "_" + geo["ch_library_prep_plate"],
            "subject_id": "GSE215865:" + subj,
            "lc_status": status,
            "case_definition": "study_specific",
            "sample_matrix": "whole_blood",
            "sex": geo["ch_Sex"].str.lower(),
            "timepoint": geo["ch_days_since_first_sample"],
            "gsm": geo["gsm"],
            "notes": notes,
        }
    )
    return finish(df, "GSE215865")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default="data/meta")
    ap.add_argument("--taco", required=True, help="TACO Data/411_DCEG_Paths_average.csv")
    a = ap.parse_args()
    d = Path(a.meta_dir)
    built = {"GSE215865": gse215865(pd.read_csv(d / "GSE215865_geo_samples.csv"), pd.read_csv(a.taco, usecols=["Subject_ID", "Long_COVID"]))}
    for sid, df in built.items():
        df.to_csv(d / f"{sid}_meta.csv", index=False)
        print(f"{sid}: {len(df)} samples, {df['subject_id'].nunique()} subjects; lc_status {df['lc_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
