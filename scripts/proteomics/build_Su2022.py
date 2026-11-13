"""Build the Su2022 Olink cohort (Su Y, ..., Heath JR 2022 Cell, "Multiple early factors anticipate
post-acute COVID-19 sequelae", doi:10.1016/j.cell.2022.01.014, PMID 35216672, PMC8786632;
Mendeley Data doi:10.17632/96v329bg7g.1, CC BY 4.0).

Raw files (data/proteomics/Su2022/raw):
  Table S1.xlsx   S1.1 INCOV (one row per blood draw: days since symptom onset, encounter location,
                  WHO ordinal scale, age, sex, BMI), S1.2 healthy controls, S1.3 PASC data (one row
                  per INCOV patient with a T3 symptom survey), S1.4 HAARVI
  Table S2.xlsx   S2.1 Proteomics: one row per sample, columns <gene symbol>_<panel> with
                  panels CVD2 (Cardiovascular II), INF (Inflammation), IRE (Immune Response),
                  MET (Metabolism), ODA (Organ Damage). Values are the authors' batch corrected,
                  age, sex and BMI adjusted Z scores (relative to matched healthy controls), not NPX.

Only S2.1 proteomics columns are used. INCOV age, sex and BMI come from Table S1.1 (clinical
baseline values); the BMI column of Table S2.1 differs for many draws and is filled where S1.1 has
none, so it is kept separately as ``bmi_table_s2``. The assays are named by gene symbol; they are mapped to
UniProt with the HGNC complete set (data/ref/hgnc_complete_set.txt): approved symbol first, then a
unique previous symbol, then a unique alias symbol. A symbol whose HGNC entry lists several UniProt
accessions is treated as ambiguous and dropped (not guessed). Assays with two symbols ("IL27,EBI3")
keep all accessions (multi-protein key). The mapping is passed as the ``annotation`` of the olink
adapter (wide form); the adapter averages assays of one protein measured on several panels.

No vendor QC flags are deposited, so no values are masked. No normalisation, imputation or
cross-plate correction is applied here (the deposited values are already batch corrected by the
authors).

Labels: LC = INCOV patient with >= 1 PASC symptom reported at T3 in Table S1.3 (the authors' "1-3
PASC" and ">=4 PASC" groups); COVID_recovered = patient with no PASC reported (the authors' "no
PASC" group: no individual symptom = 1 and all four symptom groups = 0). T1 and T2 draws are acute
and carry the patient's later outcome; T3 draws are post-acute. Samples PD1/PD2 are in S2.1 but not
in Table S1.1 and are not described in the paper: phase and label unknown. Healthy controls (Arivale
wellness program) are never_infected.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ROOT, raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Su2022"
HGNC = ROOT / "data" / "ref" / "hgnc_complete_set.txt"
PANELS = {
    "CVD2": "Olink Target 96 Cardiovascular II",
    "INF": "Olink Target 96 Inflammation",
    "IRE": "Olink Target 96 Immune Response",
    "MET": "Olink Target 96 Metabolism",
    "ODA": "Olink Target 96 Organ Damage",
}
PLATFORM = "Olink Target 96 (Cardiovascular II, Inflammation, Immune Response, Metabolism, Organ Damage); authors' Z scores"
PASC_GROUPS = ["Respiratory viral", "Neurological", "Anosmia/Dysgeusia", "Gastrointestinal"]

def read_xlsx_sheet(path: Path, name: str, header: int = 0) -> pd.DataFrame:
    """One worksheet as a DataFrame (row ``header`` as column names); "NA" and "None" are missing."""
    df = pd.read_excel(path, sheet_name=name, header=header, na_values=["NA", "None"], keep_default_na=True)
    df.columns = [str(c).strip() for c in df.columns]
    return df


# ---------------------------------------------------------------------------------- symbols
def hgnc_map(symbols) -> tuple:
    """symbol -> UniProt accession via HGNC; returns (dict, per-symbol route table)."""
    h = pd.read_csv(HGNC, sep="\t", dtype=str, usecols=["symbol", "alias_symbol", "prev_symbol", "uniprot_ids"])
    approved = h.set_index("symbol")["uniprot_ids"]

    def index(col):
        d = {}
        for sym, vals in zip(h["symbol"], h[col]):
            if isinstance(vals, str):
                for s in vals.split("|"):
                    d.setdefault(s, set()).add(sym)
        return d

    prev, alias = index("prev_symbol"), index("alias_symbol")
    out, rows = {}, []
    for s in sorted(set(symbols)):
        if s in approved.index:
            route, hits = "symbol", [s]
        elif s in prev:
            route, hits = "prev_symbol", sorted(prev[s])
        elif s in alias:
            route, hits = "alias_symbol", sorted(alias[s])
        else:
            route, hits = "fail", []
        accs = sorted({a for x in hits if isinstance(approved.get(x), str) for a in approved[x].split("|")})
        if route == "fail" or not accs:
            status = "fail"
        elif len(hits) > 1 or len(accs) > 1:
            status = "ambiguous"
        else:
            status = "mapped"
            out[s] = accs[0]
        rows.append({"symbol": s, "route": route, "hgnc_symbols": "|".join(hits), "uniprot_ids": "|".join(accs), "status": status})
    return out, pd.DataFrame(rows)


def assay_annotation(assays) -> tuple:
    syms = {a: a.rsplit("_", 1)[0].split(",") for a in assays}
    m, table = hgnc_map([s for v in syms.values() for s in v])
    ann = []
    for a, ss in syms.items():
        accs = [m.get(s) for s in ss]
        ann.append({"assay": a, "UniProt": ";".join(accs) if all(accs) else np.nan})
    return pd.DataFrame(ann), table


# ---------------------------------------------------------------------------------- labels
def severity_band(wos: float) -> str:
    if pd.isna(wos):
        return "unknown"
    if wos <= 2:
        return "ambulatory"
    if wos == 3:
        return "hospitalised_no_oxygen"
    if wos == 4:
        return "hospitalised_oxygen"
    return "icu_or_ventilated"


def wos_numeric(v) -> float:
    s = str(v).strip()
    if s in ("1 or 2", "<=2"):
        return 2.0
    return pd.to_numeric(s, errors="coerce")


def pasc_table(raw: Path) -> pd.DataFrame:
    s = read_xlsx_sheet(raw / "Table S1.xlsx", "S1.3 PASC data").set_index("Study subject ID")
    s = s.apply(pd.to_numeric, errors="coerce")
    sym = [c for c in s.columns if c not in PASC_GROUPS]
    out = pd.DataFrame(index=s.index)
    out["n_pasc_symptoms"] = (s[sym] == 1).sum(axis=1)
    out["pasc_symptoms"] = s[sym].apply(lambda r: ";".join(c for c in sym if r[c] == 1), axis=1)
    out["pasc_groups"] = s[PASC_GROUPS].apply(lambda r: ";".join(c for c in PASC_GROUPS if r[c] == 1), axis=1)
    no_pasc = (out["n_pasc_symptoms"] == 0) & (s[PASC_GROUPS] == 0).all(axis=1)
    out["lc_status"] = np.where(out["n_pasc_symptoms"] >= 1, "LC", np.where(no_pasc, "COVID_recovered", "unknown"))
    out["pasc_count_group"] = np.select(
        [out["n_pasc_symptoms"] >= 4, out["n_pasc_symptoms"] >= 1, no_pasc], [">=4 PASC", "1-3 PASC", "no PASC"], "unknown"
    )
    # binary PASC group of the authors' XGBoost model: 1 if respiratory viral, neurological or
    # anosmia/dysgeusia = 1; 0 if all three = 0
    g3 = s[PASC_GROUPS[:3]]
    out["pasc_group_binary_xgboost"] = np.where((g3 == 1).any(axis=1), "1", np.where((g3 == 0).all(axis=1), "0", ""))
    return out


# ---------------------------------------------------------------------------------- build
def build() -> list:
    raw = raw_dir(STUDY)
    P = read_xlsx_sheet(raw / "Table S2.xlsx", "S2.1 Proteomics")
    id_cols = ["Patient Subject ID", "Blood Draw", "Healthy or INCOV", "age", "sex", "BMI"]
    assays = [c for c in P.columns if c not in id_cols]
    is_healthy = P["Healthy or INCOV"] == "Healthy"
    hid = P["Patient Subject ID"].map(lambda x: f"{int(x)}" if isinstance(x, float) else str(x))
    P["sample_id"] = np.where(is_healthy, "HC-" + hid, P["Blood Draw"].astype(str))
    P["subject_id"] = np.where(is_healthy, "HC-" + hid, P["Patient Subject ID"].astype(str))
    if P["sample_id"].duplicated().any():
        raise ValueError("duplicated sample ids in S2.1")
    wide = P.set_index("sample_id")[assays].apply(pd.to_numeric, errors="coerce")

    ann, sym_table = assay_annotation(assays)
    adapted = to_feature_space("olink", wide, annotation=ann, assay_col="assay", uniprot_col="UniProt")

    s1 = read_xlsx_sheet(raw / "Table S1.xlsx", "S1.1 INCOV").set_index("Blood draw")
    s1["wos"] = s1["Who Ordinal Scale"].map(wos_numeric)
    s1["tp"] = s1.index.str.split("-").str[1]
    acute = s1[s1["tp"].isin(["T1", "T2"])]
    peak_wos = acute.groupby("Study subject ID")["wos"].max()
    hosp = acute.groupby("Study subject ID")["Encounter location"].apply(lambda s: s.isin(["Hospital", "ICU"]).any())
    pasc = pasc_table(raw)

    sex_map = {"female": "female", "male": "male"}
    rows = []
    for sid, r in P.set_index("sample_id").iterrows():
        subj = r["subject_id"]
        row = {"subject_id": subj, "sample_matrix": "plasma", "platform": PLATFORM, "collection_country": "United States"}
        notes = []
        if r["Healthy or INCOV"] == "Healthy":
            row.update(
                lc_status="never_infected",
                case_definition="not_applicable",
                control_type="never_infected_unspecified",
                sample_phase="not_applicable",
                timepoint="single_visit",
                age=r["age"],
                sex=sex_map.get(str(r["sex"]).lower(), "unknown"),
                bmi=r["BMI"],
                cohort_arm="healthy_Arivale",
                bmi_table_s2=r["BMI"],
            )
            notes.append("healthy control from the Arivale wellness program (paper Methods); collection dates not given")
        else:
            tp = sid.split("-")[1]
            row.update(timepoint=tp, cohort_arm="INCOV", bmi_table_s2=r["BMI"])
            outcome = pasc.loc[subj] if subj in pasc.index else None
            for c in ["n_pasc_symptoms", "pasc_symptoms", "pasc_groups", "pasc_count_group", "pasc_group_binary_xgboost"]:
                row[c] = outcome[c] if outcome is not None else np.nan
            row["peak_acute_wos"] = peak_wos.get(subj, np.nan)
            row["acute_severity"] = severity_band(peak_wos.get(subj, np.nan))
            row["hospitalised"] = bool(hosp[subj]) if subj in hosp.index else np.nan
            label = outcome["lc_status"] if outcome is not None else "unknown"
            if sid in s1.index:
                m = s1.loc[sid]
                row.update(age=m["Age at baseline"], sex=sex_map.get(str(m["Sex"]).lower(), "unknown"), bmi=m["BMI at baseline"])
                row["wos_at_draw"] = m["Who Ordinal Scale"]
                row["encounter_location"] = m["Encounter location"]
                d = pd.to_numeric(m["Observation days since onset of symptoms"], errors="coerce")
                if pd.notna(d) and d < 0:
                    notes.append(f"S1.1 days since symptom onset is negative ({d}); left empty")
                    d = np.nan
                row.update(days_since_infection=d, time_anchor="symptom_onset" if pd.notna(d) else None)
            else:
                m = s1[s1["Study subject ID"] == subj].iloc[0]
                row.update(age=m["Age at baseline"], sex=sex_map.get(str(m["Sex"]).lower(), "unknown"), bmi=m["BMI at baseline"])
                notes.append("blood draw not in Table S1.1; baseline age, sex, BMI of the patient from Table S1.1")
            if tp in ("T1", "T2"):
                row.update(sample_phase="acute", lc_status=label, control_type="not_applicable")
                row["case_definition"] = "study_specific" if label != "unknown" else "not_applicable"
                if label == "unknown":
                    notes.append("no T3 PASC survey for this patient in Table S1.3")
            elif tp == "T3":
                row.update(sample_phase="post_acute", lc_status=label)
                row["control_type"] = "COVID_recovered" if label == "COVID_recovered" else "not_applicable"
                row["case_definition"] = "study_specific" if label != "unknown" else "not_applicable"
                if label == "unknown":
                    notes.append("T3 draw without PASC status in Table S1.3")
            else:
                row.update(sample_phase="unknown", lc_status="unknown", control_type="not_applicable", case_definition="not_applicable")
                notes.append(f"draw {tp} not described in the paper or Table S1.1; phase and label unknown (patient outcome kept in pasc columns)")
        row["notes"] = "; ".join(notes)
        rows.append(row)

    meta = pd.DataFrame(rows, index=P["sample_id"].values)
    meta.index.name = "sample_id"
    meta = meta.loc[adapted.X.index]

    st = sym_table.copy()
    extra = {
        "n_assays": len(assays),
        "assays_without_uniprot": adapted.notes.get("assays_without_uniprot"),
        "assays_without_uniprot_names": ",".join(ann.loc[ann["UniProt"].isna(), "assay"]),
        "n_symbols": len(st),
        "symbols_by_route": ",".join(f"{k}={v}" for k, v in st["route"].value_counts().items()),
        "symbols_by_status": ",".join(f"{k}={v}" for k, v in st["status"].value_counts().items()),
        "symbols_not_approved": ",".join(f"{a}->{b}" for a, b in st.loc[st["route"] != "symbol", ["symbol", "hgnc_symbols"]].values),
        "symbols_ambiguous": ",".join(f"{a}({b})" for a, b in st.loc[st["status"] == "ambiguous", ["symbol", "uniprot_ids"]].values),
        "features_from_several_panels": int((adapted.features["n_sources"] > 1).sum()),
        "values_missing": int(adapted.X.isna().sum().sum()),
        "n_timepoint": ",".join(f"{k}={v}" for k, v in meta["timepoint"].value_counts().sort_index().items()),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
