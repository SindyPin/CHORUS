"""Harmonised clinical metadata schema for pooled Long COVID omics cohorts.

The schema is deliberately minimal. Every field has a controlled vocabulary or a unit so
that cohorts curated by different people (or by an LLM assisted curation step followed by
human review) land in the same representation. Fields marked required must be present for
a sample to enter a pooled model. Everything else is optional and recorded as missing.

Vocabularies follow, where one exists, an external standard:
  case definition      WHO 2021 (symptoms at 3 months lasting 2 months), NICE NG188 (12 weeks),
                       CDC (4 weeks), NASEM 2024, RECOVER PASC score, self report, study specific
  acute severity       WHO Clinical Progression Scale (0 to 10) collapsed to four bands
  time anchor          symptom onset, positive test, hospital admission, hospital discharge
  sample matrix        plasma, serum, whole blood, PBMC, saliva, urine, BAL, CSF, other
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import pandas as pd

LC_STATUS = {"LC", "COVID_recovered", "never_infected", "disease_comparator", "unknown"}
CASE_DEFINITIONS = {
    "WHO_2021",
    "NICE_NG188",
    "CDC_4wk",
    "NASEM_2024",
    "RECOVER_PASC_score",
    "self_report",
    "study_specific",
    "not_applicable",
}
TIME_ANCHORS = {"symptom_onset", "positive_test", "admission", "discharge", "unknown"}
SEVERITY_BANDS = {"ambulatory", "hospitalised_no_oxygen", "hospitalised_oxygen", "icu_or_ventilated", "unknown"}
SEX = {"female", "male", "other", "unknown"}
INFERRED_SEX = {"female", "male", "ambiguous"}  # from XIST and RPS4Y1 expression (chorus.sex)
MATRIX = {"plasma", "serum", "whole_blood", "PBMC", "saliva", "urine", "BAL", "CSF", "other"}
VARIANT_ERA = {"pre_alpha", "alpha", "delta", "omicron", "mixed", "unknown"}
LAYERS = {"bulk_rna", "sc_rna", "olink", "ms_proteomics", "metabolomics", "dnam", "genomics", "other"}


@dataclass
class ClinicalSchema:
    """Column names and vocabularies of the harmonised metadata table."""

    required: List[str] = field(
        default_factory=lambda: [
            "sample_id",
            "subject_id",
            "study_id",
            "lc_status",
            "case_definition",
            "sample_matrix",
            "omics_layer",
        ]
    )
    optional: List[str] = field(
        default_factory=lambda: [
            "days_since_infection",
            "time_anchor",
            "acute_severity",
            "hospitalised",
            "age",
            "sex",
            "inferred_sex",
            "inferred_sex_method",
            "bmi",
            "vaccinated_before_infection",
            "variant_era",
            "comorbidity_count",
            "symptom_domains",
            "lc_definition_strict_who",
            "platform",
            "collection_country",
            "overlaps_with_study",
            "notes",
        ]
    )
    vocab: Dict[str, set] = field(
        default_factory=lambda: {
            "lc_status": LC_STATUS,
            "case_definition": CASE_DEFINITIONS,
            "time_anchor": TIME_ANCHORS,
            "acute_severity": SEVERITY_BANDS,
            "sex": SEX,
            "inferred_sex": INFERRED_SEX,
            "sample_matrix": MATRIX,
            "variant_era": VARIANT_ERA,
            "omics_layer": LAYERS,
        }
    )

    @property
    def columns(self) -> List[str]:
        return self.required + self.optional

    def template(self, n: int = 0) -> pd.DataFrame:
        """Empty metadata table with all schema columns, for curators."""
        return pd.DataFrame({c: pd.Series([None] * n, dtype="object") for c in self.columns})


def validate_metadata(meta: pd.DataFrame, schema: Optional[ClinicalSchema] = None, strict: bool = True) -> List[str]:
    """Validate a metadata table against the schema.

    Returns a list of human readable problems. With ``strict=True`` an exception is raised
    when any required column is missing or a vocabulary is violated.
    """
    schema = schema or ClinicalSchema()
    problems: List[str] = []
    for col in schema.required:
        if col not in meta.columns:
            problems.append(f"missing required column: {col}")
    for col, vocab in schema.vocab.items():
        if col in meta.columns:
            bad = set(meta[col].dropna().astype(str)) - vocab
            if bad:
                problems.append(f"column {col} has values outside vocabulary: {sorted(bad)[:10]}")
    if "sample_id" in meta.columns and meta["sample_id"].duplicated().any():
        problems.append("duplicated sample_id values")
    if "days_since_infection" in meta.columns:
        d = pd.to_numeric(meta["days_since_infection"], errors="coerce")
        if (d < 0).any():
            problems.append("negative days_since_infection")
    if "age" in meta.columns:
        a = pd.to_numeric(meta["age"], errors="coerce")
        if ((a < 0) | (a > 120)).any():
            problems.append("age outside 0 to 120")
    if strict and problems:
        raise ValueError("metadata validation failed:\n  " + "\n  ".join(problems))
    return problems


def binary_label(meta: pd.DataFrame, positive: str = "LC", negative: str = "COVID_recovered") -> pd.Series:
    """Map lc_status to a binary outcome; samples outside the contrast are NaN.

    The primary contrast in CHORUS is LC versus COVID recovered, because never infected
    controls mix the Long COVID signal with the signal of having had COVID at all.
    """
    y = pd.Series(float("nan"), index=meta.index)
    y[meta["lc_status"] == positive] = 1.0
    y[meta["lc_status"] == negative] = 0.0
    return y


def strict_who_subset(meta: pd.DataFrame) -> pd.Series:
    """Boolean mask for a sensitivity analysis restricted to WHO compatible cases.

    A sample passes when it is a control, or when it is an LC case either labelled under
    WHO_2021 or sampled at least 90 days after infection under another definition.
    """
    is_case = meta["lc_status"] == "LC"
    who = meta.get("case_definition", pd.Series("unknown", index=meta.index)) == "WHO_2021"
    days = pd.to_numeric(meta.get("days_since_infection", pd.Series(float("nan"), index=meta.index)), errors="coerce")
    return (~is_case) | who | (days >= 90)
