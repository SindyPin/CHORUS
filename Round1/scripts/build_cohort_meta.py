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


def gse226260_wb(geo: pd.DataFrame) -> pd.DataFrame:
    """Aid 2025 whole blood subset. Titles are <subject>.<days>.<n>[PASC|Recovered].

    The second field agrees with the GEO draw window in every sample (Acute 15 to 30 days,
    3_6M 79 to 177, 6_12M 181 to 363, >12M 368 to 1678), so it is read as days since infection;
    the paper counts days "after initial infection" without naming the anchor.
    """
    g = geo[geo["ch_tissue"] == "Whole blood"]
    parts = g["title"].str.split(".", expand=True)
    df = pd.DataFrame(
        {
            "sample_id": g["gsm"],
            "subject_id": "GSE226260:" + parts[0],
            "lc_status": g["ch_current_long_covid_status"].map({"PASC": "LC", "Recovered": "COVID_recovered"}),
            "case_definition": "study_specific",
            "sample_matrix": "whole_blood",
            "days_since_infection": pd.to_numeric(parts[1]),
            "time_anchor": "unknown",
            "timepoint": g["ch_time_of_blood_draw"],
            "gsm": g["gsm"],
        }
    )
    return finish(df, "GSE226260_WB")


def gse226260_pbmc(geo: pd.DataFrame) -> pd.DataFrame:
    """Aid 2025 PBMC subset; timepoint is the collection date (no infection date in GEO)."""
    g = geo[geo["ch_cell_type"] == "PBMC"]
    df = pd.DataFrame(
        {
            "sample_id": g["gsm"],
            "subject_id": "GSE226260:" + g["ch_subject_id"].astype(int).astype(str),
            "lc_status": g["ch_pasc_status"].map({"PASC": "LC", "NOPASC": "COVID_recovered"}),
            "case_definition": "study_specific",
            "sample_matrix": "PBMC",
            "age": pd.to_numeric(g["ch_age"], errors="coerce"),
            "sex": g["ch_gender"].str.lower(),
            "timepoint": g["ch_timepoint"],
            "gsm": g["gsm"],
        }
    )
    df.loc[g["title"].str.contains("-s049_"), "notes"] = "carried the malformed tail of GSE226260_combined.ExprData; check before use"
    return finish(df, "GSE226260_PBMC")


def gse224615(geo: pd.DataFrame) -> pd.DataFrame:
    """Yin 2024 (LIINC): all PBMC at 8 months; the paper applies the WHO 3 month definition."""
    df = pd.DataFrame(
        {
            "sample_id": geo["gsm"],
            "subject_id": "GSE224615:" + geo["ch_pid"].astype(str),
            "lc_status": geo["ch_lc_status"].map({"LC": "LC", "Non-LC": "COVID_recovered"}),
            "case_definition": "WHO_2021",
            "sample_matrix": "PBMC",
            "days_since_infection": 243,
            "time_anchor": "unknown",
            "sex": geo["ch_Sex"].str.lower(),
            "timepoint": geo["ch_time_(collection_point)"],
            "gsm": geo["gsm"],
            "notes": "days_since_infection is 8 months (243 days), the cohort visit, not a per sample date",
        }
    )
    return finish(df, "GSE224615")


def gse251849(geo: pd.DataFrame, matrix_columns: list) -> pd.DataFrame:
    """Greene 2024: the count matrix uses its own ids, so labels come from the id prefix.

    Group sizes agree with GEO (7 Control, 5 Convalescent, 11 LongCOVID), but which column is
    which GSM within a group is not recorded, so per sample sex and the GSM are left empty.
    """
    ids = pd.Series(matrix_columns)
    group = ids.str.extract(r"^(Control|Convalescent|LongCOVID)")[0]
    expected = {"Control": "Healthy", "Convalescent": "COVID convalescent"}
    for k, v in expected.items():
        assert (group == k).sum() == (geo["ch_group"] == v).sum(), k
    assert (group == "LongCOVID").sum() == geo["ch_group"].str.startswith("Long COVID").sum()
    df = pd.DataFrame(
        {
            "sample_id": ids,
            "subject_id": "GSE251849:" + ids,
            "lc_status": group.map({"Control": "never_infected", "Convalescent": "COVID_recovered", "LongCOVID": "LC"}),
            "case_definition": "study_specific",
            "sample_matrix": "PBMC",
            "notes": "matrix id; GSM and sex not mappable within group",
        }
    )
    return finish(df, "GSE251849")


def gse270045(geo: pd.DataFrame) -> pd.DataFrame:
    """Shahbaz 2025: the paper states the GEO 'Healthy Control' samples recovered from COVID-19."""
    status = np.where(geo["title"].str.startswith("Long Covid"), "LC", np.where(geo["title"].str.startswith("Healthy Control"), "COVID_recovered", "unknown"))
    df = pd.DataFrame(
        {
            "sample_id": geo["ch_sample_id"],
            "subject_id": "GSE270045:" + geo["ch_sample_id"],
            "lc_status": status,
            "case_definition": "study_specific",
            "sample_matrix": "whole_blood",
            "gsm": geo["gsm"],
        }
    )
    return finish(df, "GSE270045")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta-dir", default="data/meta")
    ap.add_argument("--counts-dir", default="data/counts")
    ap.add_argument("--taco", required=True, help="TACO Data/411_DCEG_Paths_average.csv")
    a = ap.parse_args()
    d = Path(a.meta_dir)
    geo = {g: pd.read_csv(d / f"{g}_geo_samples.csv") for g in ["GSE215865", "GSE226260", "GSE224615", "GSE251849", "GSE270045"]}
    gse251849_cols = pd.read_csv(Path(a.counts_dir) / "GSE251849" / "GSE251849_Counts.txt.gz", sep="\t", index_col=0, nrows=0).columns.tolist()
    built = {
        "GSE226260_WB": gse226260_wb(geo["GSE226260"]),
        "GSE226260_PBMC": gse226260_pbmc(geo["GSE226260"]),
        "GSE224615": gse224615(geo["GSE224615"]),
        "GSE251849": gse251849(geo["GSE251849"], gse251849_cols),
        "GSE270045": gse270045(geo["GSE270045"]),
        "GSE215865": gse215865(geo["GSE215865"], pd.read_csv(a.taco, usecols=["Subject_ID", "Long_COVID"])),
    }
    inferred = d / "inferred_sex.csv"  # written by scripts/infer_sex.py; kept when the tables are rebuilt
    inferred = pd.read_csv(inferred) if inferred.exists() else None
    for sid, df in built.items():
        if inferred is not None:
            m = inferred[inferred["study_id"] == sid].set_index("sample_id")
            df["inferred_sex"] = df["sample_id"].map(m["inferred_sex"])
            df["inferred_sex_method"] = df["sample_id"].map(m["provenance"])
        df.to_csv(d / f"{sid}_meta.csv", index=False)
        print(f"{sid}: {len(df)} samples, {df['subject_id'].nunique()} subjects; lc_status {df['lc_status'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
