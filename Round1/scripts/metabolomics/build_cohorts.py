"""Build the public metabolomics cohorts with an LC versus COVID_recovered contrast
(docs/metabolomics_cohorts.md). Each builder reads data/metabolomics/<study>/raw and writes, through
scripts/proteomics/common.write_cohort,

  data/metabolomics/<study>/<study>_matrix.csv    metabolite names as deposited (rows) x samples
  data/metabolomics/<study>/<study>_features.csv  the metabolomics adapter's feature table
  data/meta/<study>_meta.csv                      ClinicalSchema metadata

Values are log2 of the deposited values (areas, concentrations or ratios); zero or negative values
are missing. Subject ids are prefixed by the study (``Martinez2024_C1``), because several deposits
number their subjects C1, L1 and so on. Names keep a ``<platform>::`` prefix where a study measured several platforms. The
shared identifiers are added by scripts/metabolomics/map_ids.py.

  python scripts/metabolomics/build_cohorts.py [Martinez2024 Lopez2023 Saito2024 Holmes2021]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proteomics"))
from common import ROOT, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

OUT = ROOT / "data" / "metabolomics"


def raw(study: str) -> Path:
    return OUT / study / "raw"


def base_meta(index, **cols) -> pd.DataFrame:
    m = pd.DataFrame(index=pd.Index(index, name="sample_id"))
    for k, v in cols.items():
        m[k] = v
    return m


# ------------------------------------------------------------------------------ Martinez 2024
MARTINEZ_ANALYSES = {  # analysis id -> platform prefix
    "AN005077": "LCMSPOS",
    "AN005078": "LCMSNEG",
    "AN005079": "GCMS",
    "AN005080": "CEMSPOS",
    "AN005081": "CEMSNEG",
}


def read_mwtab_data(path: Path) -> pd.DataFrame:
    """MS_METABOLITE_DATA block of an mwTab text file: samples x metabolite names."""
    lines = path.read_text().splitlines()
    a, b = lines.index("MS_METABOLITE_DATA_START"), lines.index("MS_METABOLITE_DATA_END")
    rows = [l.split("\t") for l in lines[a + 1 : b]]
    samples = rows[0][1:]
    data = [r for r in rows[1:] if r[0] != "Factors"]
    df = pd.DataFrame([r[1 : len(samples) + 1] for r in data], index=[r[0].strip() for r in data], columns=samples)
    return df.apply(pd.to_numeric, errors="coerce").T


def build_Martinez2024() -> dict:
    """Metabolomics Workbench ST003103 (Martinez S et al. 2024 J Proteome Res, PMID 38566450): plasma,
    40 Long COVID and 40 fully recovered from mild COVID-19, five MS platforms; QC samples dropped."""
    blocks = []
    for an, pref in MARTINEZ_ANALYSES.items():
        d = read_mwtab_data(raw("Martinez2024") / f"ST003103_{an}.mwtab.txt")
        d.columns = [f"{pref}::{c}" for c in d.columns]
        blocks.append(d)
    X = pd.concat(blocks, axis=1)
    X = X.loc[[s for s in X.index if s[0] in "CL" and not s.startswith("QC")]]
    X = X.loc[sorted(X.index, key=lambda s: (s[0], int(s[1:])))]
    adapted = to_feature_space("metabolomics", X, log2=True)
    lc = np.where(X.index.str.startswith("L"), "LC", "COVID_recovered")
    meta = base_meta(
        X.index,
        subject_id="Martinez2024_" + X.index,
        lc_status=lc,
        case_definition="study_specific",
        sample_matrix="plasma",
        sample_phase="post_acute",
        control_type=np.where(lc == "COVID_recovered", "COVID_recovered", "not_applicable"),
        platform="untargeted CE-ESI(+/-)-TOF-MS, GC-Q-MS and LC-ESI(+/-)-QTOF-MS lipidomics (Metabolomics Workbench ST003103)",
        notes="Workbench factor 'Long COVID' or 'Recovered from acute COVID-19'; the abstract compares Long COVID with individuals fully recovered from acute mild COVID-19",
    )
    n_missing_gc = int(X[[c for c in X.columns if c.startswith("GCMS::")]].isna().all(axis=1).sum())
    return write_cohort("Martinez2024", adapted, meta, {"n_analytes": X.shape[1], "samples_without_GCMS": n_missing_gc}, root=OUT)


# ------------------------------------------------------------------------------ Lopez 2023
LOPEZ_FIRST, LOPEZ_LAST = "Creatinine", "Octadecenoylcarnitine"  # the metabolite columns of the sheet


def build_Lopez2023() -> dict:
    """Mendeley Data 10.17632/8zfdjsypd8.1 (López-Hernández Y et al. 2023 Sci Rep, PMID 37528111):
    plasma, targeted LC-MS/MS and FIA-MS/MS (108 metabolites, µM). POST-COVID class A and class B
    are LC, POST-COVID RECOVERED is COVID_recovered (about two years after infection); COVID-19 are
    acute samples of the same patients (Folio <n>.1 pairs with post-COVID Folio <n>); CONTROL are
    negative (healthy) controls."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = pd.read_excel(raw("Lopez2023") / "DATA AVAILABILITY.xls", engine="xlrd", dtype={"Folio": str})
    cols = list(v.columns)
    mets = cols[cols.index(LOPEZ_FIRST) : cols.index(LOPEZ_LAST) + 1]
    folio = v["Folio"].astype(str).str.strip()
    group, sub = v["GROUP"].astype(str).str.strip(), v["SUBGROUP"].astype(str).str.strip()
    subject = ("Lopez2023_" + pd.Series(np.where(group == "CONTROL", folio, "P" + folio.str.split(".").str[0]))).values
    sid = np.where(group == "CONTROL", folio, np.where(group == "COVID-19", "P" + folio.str.replace(".", "_", regex=False) + "_acute", "P" + folio + "_post"))
    X = v[mets].set_index(pd.Index(sid, name="sample_id"))
    X.columns = [str(c).replace("\xa0", " ").strip() for c in X.columns]
    adapted = to_feature_space("metabolomics", X, log2=True)
    lc = np.select([group == "CONTROL", group == "COVID-19", sub.isin(["CLASS A", "CLASS B"]), sub == "RECOVERED"], ["never_infected", "LC", "LC", "COVID_recovered"], "unknown")
    post_status = dict(zip(subject[group == "POST-COVID"], lc[group == "POST-COVID"]))
    lc = np.where(group == "COVID-19", [post_status.get(s, "unknown") for s in subject], lc)  # acute samples carry the later outcome
    meta = base_meta(
        X.index,
        subject_id=subject,
        lc_status=lc,
        case_definition=np.where(group == "CONTROL", "not_applicable", "study_specific"),
        sample_matrix="plasma",
        sample_phase=np.select([group == "CONTROL", group == "COVID-19"], ["not_applicable", "acute"], "post_acute"),
        control_type=np.select([group == "CONTROL", sub == "RECOVERED"], ["never_infected_contemporary", "COVID_recovered"], "not_applicable"),
        timepoint=np.select([group == "COVID-19", group == "POST-COVID"], ["acute", "2y"], "single_visit"),
        age=pd.to_numeric(v["AGE"], errors="coerce").values,
        platform="targeted LC-MS/MS and FIA-MS/MS (TMIC), µM",
        collection_country="Mexico",
        lc_subgroup=np.where(group == "POST-COVID", sub, ""),
        group_deposited=group.values,
        sex_code_deposited=v["SEX"].values,
        notes="sex coded 1/2 in the deposit without a key; LC = post-COVID class A or class B (paper: arbitrarily defined classes by number and type of persisting symptoms), recovered = post-COVID recovered (non-long COVID)",
    )
    return write_cohort("Lopez2023", adapted, meta, {"n_analytes": len(mets)}, root=OUT)


# ------------------------------------------------------------------------------ Saito 2024
def build_Saito2024() -> dict:
    """Supplementary Table 2 (DataSheet_2.xlsx) of Saito S et al. 2024 Front Immunol (PMID 38304426,
    PMC10830702): plasma, dansyl chemical isotope labelling LC-MS, peak ratio of each sample to a
    13C-labelled pool. Peak pairs with a Tier 1 or Tier 2 identification are kept, named by the
    identified compound. Sample groups in the file: A (acute ICU, 15), C (15), HC (healthy, 15),
    L (Long COVID, 30). The clinical table (Supplementary Table 1) names the groups LC, R (recovered),
    HC and A with the same sizes; C is taken as R. Sample numbers in the two tables differ (L1 to L30
    versus LC1 to LC32 with gaps), so per sample clinical data are not linked."""
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        v = pd.read_excel(raw("Saito2024") / "DataSheet_2.xlsx", header=None)
    files = v.iloc[1]
    scols = [j for j in range(v.shape[1]) if isinstance(files.iloc[j], str) and files.iloc[j].startswith("IsoMS_results_")]
    names = [files.iloc[j].split("IsoMS_results_")[1].split("_")[0] for j in scols]
    body = v.iloc[4:]
    tier = body[9].astype(str)
    comp = np.where(tier == "Tier 1", body[12], np.where(tier == "Tier 2", body[17], None))
    ext = np.where(tier == "Tier 1", body[11], np.where(tier == "Tier 2", body[16], None))
    keep = pd.Series(comp, index=body.index).notna() & tier.isin(["Tier 1", "Tier 2"])
    vals = body.loc[keep, scols].apply(pd.to_numeric, errors="coerce")
    X = vals.T
    X.index = names
    X.columns = [str(c).strip() for c in pd.Series(comp, index=body.index)[keep]]
    X = X.loc[[n for n in X.index if not n.startswith("QC")]]
    X = X.T.groupby(level=0, sort=False).first().T if X.columns.duplicated().any() else X  # repeated compound: first peak pair
    adapted = to_feature_space("metabolomics", X, log2=True)
    grp = pd.Series(X.index).str.extract(r"^([A-Z]+)")[0].values
    lc = pd.Series(grp).map({"L": "LC", "C": "COVID_recovered", "HC": "never_infected", "A": "unknown"}).values
    meta = base_meta(
        X.index,
        subject_id="Saito2024_" + X.index,
        lc_status=lc,
        case_definition=np.where(grp == "L", "study_specific", "not_applicable"),
        sample_matrix="plasma",
        sample_phase=pd.Series(grp).map({"L": "post_acute", "C": "post_acute", "HC": "not_applicable", "A": "acute"}).values,
        control_type=pd.Series(grp).map({"C": "COVID_recovered", "HC": "never_infected_unspecified"}).fillna("not_applicable").values,
        timepoint=pd.Series(grp).map({"L": "12m", "C": "12m", "HC": "single_visit", "A": "acute"}).values,
        platform="dansyl chemical isotope labelling LC-MS (Tier 1 and Tier 2 identifications), ratio to 13C pool",
        collection_country="Canada",
        group_deposited=grp,
        notes="L = Long COVID meeting ME/CFS criteria about 12 months after infection; C taken as recovered (R) by group size and order; A = acute ICU patients, outcome unknown",
    )
    (OUT / "Saito2024").mkdir(parents=True, exist_ok=True)
    dep = pd.DataFrame({"name": [str(c).strip() for c in pd.Series(comp, index=body.index)[keep]], "external_identifier": pd.Series(ext, index=body.index)[keep].astype(str).values, "tier": tier[keep].values}).drop_duplicates("name")
    dep.to_csv(OUT / "Saito2024" / "Saito2024_deposited_ids.csv", index=False)
    return write_cohort("Saito2024", adapted, meta, {"n_analytes": int(keep.sum())}, root=OUT)


# ------------------------------------------------------------------------------ Holmes 2021
def build_Holmes2021() -> dict:
    """Zenodo 10.5281/zenodo.4713391 (Holmes E et al. 2021 J Proteome Res, PMID 34009992): plasma
    of 27 non-hospitalised patients about 3 months after COVID-19 (the follow-up cohort only; the
    controls and acute patients of the paper are not in the deposit), targeted MS amino acids
    (Quantity x 2 for the dilution, as the README says; the 13C/15N labelled internal standards dropped) and tryptophan pathway (conc). The deposit
    gives four symptom composites (composite, compositeBinary, compositeRec, compositeBinaryRec) but
    not their coding, so lc_status is left unknown and both binary composites are kept."""
    r = raw("Holmes2021")
    s = pd.read_csv(r / "sampleDescriptionList.txt", sep="\t", dtype=str)
    aa = pd.concat([pd.read_csv(r / f, sep="\t", dtype=str) for f in ("targetedMS_aminoAcids.txt", "targetedMS_aminoAcids_1172.txt")])
    aa = aa[aa["SampleType"] == "Sample"]
    aa = aa[~aa["Analyte Name"].str.contains(r"13C|15N", regex=True)]  # isotope labelled internal standards
    aa["q"] = pd.to_numeric(aa["Quantity [units]"].str.replace(",", ""), errors="coerce") * 2
    A = aa.pivot_table(index="sampleID", columns="Analyte Name", values="q", aggfunc="mean")
    tr = pd.read_csv(r / "targetedMS_tryptophan.txt", sep="\t")
    tr = tr[tr["sample_type"] == "Analyte"].dropna(subset=["sampleID"])
    conc = [c for c in tr.columns if c.endswith("_conc")]
    T = tr.groupby("sampleID")[conc].mean()
    T.columns = [c[: -len("_conc")].replace("_", " ") for c in T.columns]
    A.columns = [f"AA::{c}" for c in A.columns]
    T.columns = [f"TRP::{c}" for c in T.columns]
    X = A.join(T, how="outer").reindex(s["sampleID"])
    adapted = to_feature_space("metabolomics", X, log2=True)
    meta = base_meta(
        X.index,
        subject_id="Holmes2021_" + X.index,
        lc_status="unknown",
        case_definition="study_specific",
        sample_matrix="plasma",
        sample_phase="post_acute",
        timepoint="3m",
        age=pd.to_numeric(s["age"], errors="coerce").values,
        sex=s["gender"].map({"F": "female", "M": "male"}).fillna("unknown").values,
        platform="targeted LC-MS amino acids and tryptophan pathway (Zenodo 4713391)",
        composite_binary=s["compositeBinary"].values,
        composite_binary_rec=s["compositeBinaryRec"].values,
        notes="lc_status unknown: symptom composites are deposited without their coding (compositeBinary or compositeBinaryRec)",
    )
    meta["bmi"] = pd.to_numeric(s["bmi0"], errors="coerce").values
    trp_1172 = "not read (per sample file for COV01172 has a different layout)"
    return write_cohort("Holmes2021", adapted, meta, {"n_analytes": X.shape[1], "tryptophan_COV01172": trp_1172}, root=OUT)


BUILDERS = {"Martinez2024": build_Martinez2024, "Lopez2023": build_Lopez2023, "Saito2024": build_Saito2024, "Holmes2021": build_Holmes2021}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(BUILDERS):
        row = BUILDERS[name]()
        print({k: v for k, v in row.items()})
