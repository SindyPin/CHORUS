"""Tests of scripts/acute_covid_sweep.py on small synthetic inputs: label mapping, one sample per participant, overlap
detection and the verification checks."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _sweep():
    spec = importlib.util.spec_from_file_location("acute_covid_sweep", ROOT / "scripts" / "acute_covid_sweep.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_label_mapping_and_tables():
    s = _sweep()
    d = pd.DataFrame({"geo_accession": ["G1", "G2", "G3"], "ch_disease": ["COVID19", "Healthy control", "COVID19"], "ch_sex": ["Male", "F", "x"], "ch_age": ["60", "30", "na"]})
    L = s.LABELERS["GSE152641"](d)
    assert L["infection_class"].tolist() == ["acute_covid", "never_infected", "acute_covid"]
    assert L["sex"].tolist() == ["male", "female", "unknown"] and np.isnan(L["age"].iloc[2])
    t = s.mapping_table(L)
    assert set(t["label_original"]) == {"disease=COVID19", "disease=Healthy control"} and t["n_samples"].sum() == 3
    g = pd.DataFrame({"geo_accession": ["A", "B", "C", "D"], "ch_disease_state": ["COVID-19", "COVID-19", "Healthy", "COVID-19"], "ch_severity": ["Critical", "Non-critical", "Healthy", "Critical"],
                      "ch_individual": ["p1", "p2", "h1", "p1"], "ch_time": ["T1", "T1", "T0", "T2"]})
    L = s.LABELERS["GSE213313"](g)
    assert L["severity_class"].tolist()[:2] == ["severe_critical", "mild_moderate"] and pd.isna(L["severity_class"].iloc[2])
    # COMBAT: data dictionary 0 = male, 1 = female; COVID_SEV (oxygen) is not assigned a severity class
    c = pd.DataFrame({"RNASeq_sample_ID": ["r1", "r2", "r3", None], "Source": ["COVID_SEV", "HV", "COVID_CRIT", "HV"], "COMBAT_ID": ["a", "b", "c", "d"], "Death28": [0, 0, 0, 0],
                      "COMBAT_participant_timepoint_ID": ["a-Ja001", "b-Ha001", "c-Ja003", "d-Ha001"], "Sex": [0, 1, 0, 1], "Institute": ["Oxford"] * 4, "Age": [5, 5, 5, 5]})
    L = s.LABELERS["COMBAT2022_WB"](c)
    assert len(L) == 3 and L["sex"].tolist() == ["male", "female", "male"]
    assert pd.isna(L["severity_class"].iloc[0]) and L["severity_class"].iloc[2] == "severe_critical"


def test_one_sample_per_participant_and_counts():
    s = _sweep()
    L = s._frame(sample_id=["s1", "s2", "s3", "s4", "s5"], participant_id=["p1", "p1", "p2", "p2", "p3"], order=[0, 7, 0, 3, 0],
                 label_original="x", infection_class=["excluded", "acute_covid", "acute_covid", "acute_covid", "never_infected"],
                 severity_original="x", severity_class=[None, "severe_critical", "mild_moderate", "severe_critical", None])
    K = s.one_per_participant(L)
    assert sorted(K["sample_id"]) == ["s2", "s3", "s5"]  # p1: its acute sample, not the earlier excluded one; p2: the first acute
    c = s.class_counts(L)
    assert c["infection_acute_covid"] == 2 and c["infection_never_infected"] == 1 and c["severity_severe_critical"] == 1 and c["severity_mild_moderate"] == 1
    assert s.contrasts_available(c) == [] and s.contrasts_available({**c, "infection_acute_covid": 10, "infection_never_infected": 12}, 10) == ["infection"]


def test_overlap_detection_and_canonical(tmp_path, monkeypatch):
    s = _sweep()
    monkeypatch.setattr(s, "WORK", tmp_path)
    (tmp_path / "geo").mkdir()
    pd.DataFrame({"geo_accession": [f"GSM{i}" for i in range(10)], "biosample": [f"SAMN{i}" for i in range(10)]}).to_csv(tmp_path / "geo" / "GSE1_samples.csv", index=False)
    pd.DataFrame({"geo_accession": ["GSM3", "GSM4", "GSM99"], "biosample": ["SAMN3", "SAMN4", "SAMN99"]}).to_csv(tmp_path / "geo" / "GSE2_samples.csv", index=False)
    pd.DataFrame({"geo_accession": ["GSM50"], "biosample": ["SAMN50"]}).to_csv(tmp_path / "geo" / "GSE3_samples.csv", index=False)
    monkeypatch.setattr(s, "COHORTS", {"C1": {"accession": "GSE1", "layer": "RNA"}})
    ov = s.overlaps(["C1"])
    assert ov == {"GSE2": ("C1", 4)}  # two GSM and two BioSamples shared; GSE3 shares none
    assert s.canonical("10.5281/zenodo.123") == "zenodo.123" and s.canonical("10.17632/AbC9.2") == "10.17632/abc9" and s.canonical("GSE5") == "GSE5"


def _synthetic_cohort(tmp_path, s, monkeypatch, *, shuffle_sex=False, reverse_ifn=False, duplicate=False, confound=False, late=False, seed=0):
    rng = np.random.default_rng(seed)
    n = 40
    y = np.r_[np.ones(20), np.zeros(20)].astype(int)
    sex = np.array(["female", "male"] * 20)
    genes = [f"ENSG{i:011d}" for i in range(300)] + ["ENSG00000229807", "ENSG00000129824"] + list(s.IFN_GENES.values())
    X = rng.normal(5, 1, size=(len(genes), n))
    X[300] = np.where(sex == "female", 9, 1) + rng.normal(0, 0.2, n)  # XIST
    X[301] = np.where(sex == "male", 9, 1) + rng.normal(0, 0.2, n)  # RPS4Y1
    X[302:] += (-2.0 if reverse_ifn else 2.0) * y
    X[:50] += 1.5 * y  # a signal for the cross validated AUROC
    samples = [f"S{i}" for i in range(n)]
    M = pd.DataFrame(X, index=genes, columns=samples)
    if duplicate:
        M["S1"] = M["S0"] * 1.0001
    M.loc[genes[:150], "S5"] = np.nan  # an outlier on detected features
    rep_sex = np.roll(sex, 1) if shuffle_sex else sex  # shifted by one sample: the metadata no longer match the columns
    meta = pd.DataFrame({"sample_id": samples, "subject_id": [f"p{i}" for i in range(n)], "sex": rep_sex, "age": 50,
                         "infection_class": np.where(y == 1, "acute_covid", "never_infected"), "severity_class": None, "one_sample_per_participant": True,
                         "days_from_anchor_original": np.where(y == 1, 45 if late else 10, np.nan), "in_acute_window": np.where(y == 1, True, np.nan),
                         "time_anchor": np.where(y == 1, "symptom_onset", None), "acute_order": 0,
                         "plate": np.where(y == 1, "P1", "P2") if confound else rng.choice(["P1", "P2"], n)})
    d = tmp_path / "cohort"
    d.mkdir(exist_ok=True)
    (tmp_path / "data" / "meta").mkdir(parents=True, exist_ok=True)
    (tmp_path / "build").mkdir(exist_ok=True)
    M.to_csv(d / "SYN_matrix.csv.gz")
    meta.to_csv(tmp_path / "data" / "meta" / "SYN_meta.csv", index=False)
    info = {"matrix": "cohort/SYN_matrix.csv.gz", "value_kind": "normalised (log2)", "n_matrix_samples": n, "n_labelled_samples": n, "n_in_both": n, "n_duplicated_columns": 0,
            "labelled_not_in_matrix": [], "matrix_not_labelled": [], "mapping_rate": 1.0, "n_mapped": len(genes), "n_features_in": len(genes),
            "n_targets_with_duplicates": 0, "duplicate_rule": "feature with the highest mean kept", "notes": "synthetic"}
    (tmp_path / "build" / "SYN.json").write_text(json.dumps(info))
    monkeypatch.setattr(s, "ROOT", tmp_path)
    monkeypatch.setattr(s, "WORK", tmp_path)
    monkeypatch.setattr(s, "cohort_dir", lambda cid: d)
    monkeypatch.setattr(s, "COHORTS", {"SYN": {"layer": "RNA", "platform": "synthetic", "accession": "SYN"}})
    monkeypatch.setattr(s, "PAPER_COUNTS", {"SYN": {"source": "synthetic table", "numbers": {"infection_acute_covid": 20, "infection_never_infected": 20}}})
    return s.verify_cohort("SYN").set_index("check")


def test_verification_checks(tmp_path, monkeypatch):
    s = _sweep()
    V = _synthetic_cohort(tmp_path, s, monkeypatch)
    assert V.loc["a integrity", "status"] == "pass" and V.loc["b counts", "status"] == "pass" and V.loc["c value type", "status"] == "pass"
    assert V.loc["d sex", "status"] == "pass" and V.loc["g label direction", "status"] == "pass"
    assert V.loc["e duplicates (within)", "status"] == "pass" and V.loc["f outliers", "status"] == "flag" and "S5" in V.loc["f outliers", "details"]
    assert V.loc["i confounding", "status"] == "pass" and V.loc["j timing", "status"] == "pass"
    assert V.loc["h signal (infection)", "status"] == "report" and V.loc["k demographics", "status"] == "report"
    auc = float(V.loc["h signal (infection)", "value"].split("AUROC ")[1].split(";")[0])
    assert auc > 0.9
    st, flags = s.cohort_status(V.reset_index(), pd.DataFrame(columns=["cohort_a", "cohort_b", "n_pairs_above_0.99"]), "SYN")
    assert st == "verified with flags" and any(f.startswith("f outliers") for f in flags)


@pytest.mark.parametrize("kw,check,status", [
    (dict(shuffle_sex=True), "d sex", "fail"),
    (dict(reverse_ifn=True), "g label direction", "flag"),
    (dict(duplicate=True), "e duplicates (within)", "flag"),
    (dict(confound=True), "i confounding", "flag"),
    (dict(late=True), "j timing", "flag"),
])
def test_verification_checks_detect_problems(tmp_path, monkeypatch, kw, check, status):
    s = _sweep()
    V = _synthetic_cohort(tmp_path, s, monkeypatch, **kw)
    assert V.loc[check, "status"] == status
    if status == "fail":
        st, _ = s.cohort_status(V.reset_index(), pd.DataFrame(columns=["cohort_a", "cohort_b", "n_pairs_above_0.99"]), "SYN")
        assert st == "failed"


def test_auto_screen_rules():
    s = _sweep()
    C = pd.DataFrame([
        {"accession": "GSE1", "source": "GEO", "title": "Single-cell RNA-seq of COVID-19 PBMC", "description": "", "layer": "RNA", "tissue": "PBMC", "access_type": "open", "organism": "Homo sapiens", "n_samples": "40"},
        {"accession": "GSE2", "source": "GEO", "title": "Whole blood RNA-seq in COVID-19", "description": "hospitalized", "layer": "RNA", "tissue": "whole blood", "access_type": "open", "organism": "Homo sapiens", "n_samples": "12"},
        {"accession": "GSE3", "source": "GEO", "title": "Whole blood RNA-seq in COVID-19", "description": "hospitalized patients", "layer": "RNA", "tissue": "whole blood", "access_type": "open", "organism": "Homo sapiens", "n_samples": "80"},
        {"accession": "GSE4", "source": "GEO", "title": "Whole blood RNA-seq in tuberculosis", "description": "", "layer": "RNA", "tissue": "whole blood", "access_type": "open", "organism": "Homo sapiens", "n_samples": "80"},
        {"accession": "EGAS1", "source": "controlled access archive", "title": "COVID-19 blood", "description": "", "layer": "RNA", "tissue": "blood", "access_type": "controlled (application)", "organism": "", "n_samples": None},
    ])
    A = s.auto_screen(C).set_index("accession")
    assert A.loc["GSE1", "reason"] == "single cell data" and "12 samples" in A.loc["GSE2", "reason"]
    assert A.loc["GSE3", "decision"] == "needs decision" and "not a COVID-19 dataset" in A.loc["GSE4", "reason"]
    assert "not open access" in A.loc["EGAS1", "reason"]
