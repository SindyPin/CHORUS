"""Build the Chowdhury2024 plasma proteomics cohort (Chowdhury MMH et al. 2024 Front Physiol 15:1416639,
"Impact of a tailored exercise regimen on physical capacity, systemic oxidative stress and plasma proteomic
profile in PCC patients", doi:10.3389/fphys.2024.1416639, PMID 39234310, PMC11371593; PRIDE PXD047623, CC0).

Raw files (data/proteomics/Chowdhury2024/raw):
  diann/report.pg_matrix.tsv   DIA-NN 1.8.1 protein group matrix (55 runs; extracted from diann.7z)
  diann/report.stats.tsv       per run identification counts (kept as cohort columns)
  metadata.sdrf.tsv            PRIDE SDRF (51 rows: 21 PCC pre, 21 PCC post, 9 no-PCC controls)

Values: DIA-NN 1.8.1 pg_matrix quantities (MaxLFQ over DIA-NN normalised precursor quantities, 1% FDR,
library-free FASTA search of UniProt UP000005640 incl. TrEMBL, --pg-level 1, MBR/--reanalyse). The id column is
``Protein.Group`` (not ``Protein.Ids``); it is passed unchanged to the ms_proteomics adapter, which keys each
group by its first accession (DIA-NN sorts accessions alphabetically, so a TrEMBL A0A... accession is often
first), sets zeros missing and log2 transforms. No normalisation, no imputation, no filtering beyond DIA-NN's.

Sample mapping: pg_matrix run stem (path and ".d" removed) to the SDRF ``comment[data file]`` stem (".7z"
removed). Two SDRF rows are wrong: pre-covex26 and pre-covex27 point to DIA_Post_26/27 (already used by
post-covex26/27); they are mapped to DIA_Pre_26/27 by source name (pre-covexNN -> DIA_Pre_NN). Runs absent
from the SDRF are ``unknown``: DIA_Pre_13 / DIA_Post_13 (named like the PCC runs of participant 13; arm not
given), DIA_C59_60a and DIA_Sheela (no label). SDRF sex "H" = homme (male): the per arm male counts match the
paper (Table 2 / section 3.6: exercise 6 of 12 in proteomics, no-exercise 2 of 10; controls 3 of 9).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, swissprot_accessions, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Chowdhury2024"
PLATFORM = "LC-MS/MS diaPASEF timsTOF Pro (DIA-NN 1.8.1 library-free, MaxLFQ pg_matrix)"
ID_COL = "Protein.Group"
ANNOT_COLS = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"]
RUN_RE = re.compile(r"^DIA_(?P<name>.+)_(?P<slot>Slot\d+-\d+)_1_(?P<idx>\d+)$")
SEX = {"F": "female", "H": "male"}
LC_NOTE = ("authors' PCC group (BQC19 biobank; >=1 WHO candidate symptom 3 months after a positive PCR test; "
           "8 week intervention, blood at recruitment and one week after the last session)")
CONTROL_NOTE = ("authors' no-PCC control (Table 1: 'convalescent individuals without PCC'; volunteers from hospital "
                "workers / healthy research co-workers; single sample; infection not documented per person)")


def stem(path: str) -> str:
    return re.sub(r"\.(d|7z)$", "", str(path).replace("\\", "/").split("/")[-1])


def sdrf_table(raw: Path) -> pd.DataFrame:
    s = pd.read_csv(raw / "metadata.sdrf.tsv", sep="\t", dtype=str)
    s["file_stem"] = s["comment[data file]"].map(stem)
    s["mismatch"] = False
    # pre-covexNN rows must point to DIA_Pre_NN; two rows point to the Post file
    for i, r in s.iterrows():
        m = re.match(r"^(pre|post)-covex(\d+)$", r["source name"])
        if m:
            want = f"DIA_{m.group(1).capitalize()}_{int(m.group(2)):02d}_"
            if not r["file_stem"].startswith(want):
                s.loc[i, "mismatch"] = True
    return s


def build() -> list:
    raw = raw_dir(STUDY)
    pg = pd.read_csv(raw / "diann" / "report.pg_matrix.tsv", sep="\t")
    stats = pd.read_csv(raw / "diann" / "report.stats.tsv", sep="\t")
    stats.index = stats["File.Name"].map(stem)
    runs = {c: stem(c) for c in pg.columns if c not in ANNOT_COLS}
    table = pg[[ID_COL] + list(runs)].rename(columns=runs)
    adapted = to_feature_space("ms_proteomics", table, id_col=ID_COL, sample_cols=list(runs.values()), preferred_accessions=swissprot_accessions())

    s = sdrf_table(raw)
    by_file = {r["file_stem"]: r for _, r in s[~s["mismatch"]].iterrows()}
    fixed = {}
    for _, r in s[s["mismatch"]].iterrows():
        m = re.match(r"^(pre|post)-covex(\d+)$", r["source name"])
        want = f"DIA_{m.group(1).capitalize()}_{int(m.group(2)):02d}_"
        hits = [x for x in runs.values() if x.startswith(want)]
        if len(hits) != 1:
            raise ValueError(f"cannot resolve SDRF row {r['source name']}")
        fixed[hits[0]] = r
    if set(by_file) & set(fixed):
        raise ValueError("SDRF correction collides with a direct match")
    by_file.update(fixed)
    # arm of each PCC participant from the post sample label
    arm_of = {}
    for _, r in s.iterrows():
        if r["characteristics[disease]"].startswith("PCC post"):
            arm_of[r["characteristics[individual]"]] = "no_exercise" if "no-exercice" in r["characteristics[disease]"] else "exercise"

    rows = []
    for run in adapted.X.index:
        m = RUN_RE.match(run)
        if not m:
            raise ValueError(f"unexpected run name {run}")
        name = m.group("name")
        row = {"sample_id": run, "raw_run": run, "ms_slot": m.group("slot"), "ms_run_index": int(m.group("idx")),
               "precursors_identified": int(stats.loc[run, "Precursors.Identified"]),
               "proteins_identified": int(stats.loc[run, "Proteins.Identified"])}
        pp = re.match(r"^(Pre|Post)_(\d+)$", name)
        r = by_file.get(run)
        if r is not None:
            src = r["source name"]
            row.update({"sdrf_source_name": src, "sdrf_disease": r["characteristics[disease]"],
                        "sdrf_individual": r["characteristics[individual]"], "sdrf_assay_name": r["assay name"],
                        "sdrf_data_file": r["comment[data file]"], "sdrf_data_file_corrected": run in fixed,
                        "age": float(r["characteristics[age]"]), "sex": SEX[r["characteristics[sex]"]]})
            if src.startswith(("pre-covex", "post-covex")):
                num = int(re.search(r"(\d+)$", src).group(1))
                row.update({"subject_id": f"{STUDY}_PCC{num:02d}", "lc_status": "LC", "control_type": "not_applicable",
                            "sample_phase": "post_acute", "case_definition": "study_specific",
                            "timepoint": "pre_intervention" if src.startswith("pre") else "post_intervention",
                            "intervention_arm": arm_of[r["characteristics[individual]"]],
                            "time_anchor": "unknown", "acute_severity": "unknown", "notes": LC_NOTE})
                if run in fixed:
                    row["notes"] += f"; SDRF data file '{r['comment[data file]']}' corrected to this run by source name"
            elif src.startswith("control"):
                row.update({"subject_id": f"{STUDY}_{src}", "lc_status": "COVID_recovered", "control_type": "COVID_recovered",
                            "sample_phase": "post_acute", "case_definition": "study_specific", "timepoint": "single_visit",
                            "intervention_arm": "not_applicable", "time_anchor": "unknown", "acute_severity": "unknown",
                            "notes": CONTROL_NOTE})
            else:
                raise ValueError(f"unexpected SDRF source name {src}")
        else:
            subj = f"{STUDY}_PCC{int(pp.group(2)):02d}" if pp else f"{STUDY}_run_{name}"
            note = "run not in the SDRF; no group label deposited"
            if pp:
                note += (f"; named like the PCC {pp.group(1).lower()}-intervention runs (participant {int(pp.group(2))});"
                         " arm, age and sex not deposited")
            row.update({"subject_id": subj, "lc_status": "unknown", "case_definition": "not_applicable",
                        "sample_phase": "unknown",
                        "timepoint": f"{pp.group(1).lower()}_intervention" if pp else "unknown",
                        "intervention_arm": "unknown", "notes": note})
        rows.append(row)
    meta = pd.DataFrame(rows).set_index("sample_id")
    meta["sample_matrix"] = "plasma"
    meta["platform"] = PLATFORM
    meta["collection_country"] = "Canada"
    meta["value_type"] = "log2 DIA-NN 1.8.1 MaxLFQ protein group quantity (pg_matrix)"
    meta["diann_version"] = "1.8.1"

    first = ["subject_id", "lc_status", "case_definition", "sample_matrix", "platform", "control_type", "sample_phase",
             "timepoint", "time_anchor", "acute_severity", "age", "sex", "collection_country", "notes"]
    meta = meta[first + [c for c in meta.columns if c not in first]]

    missing_sdrf = sorted(set(s["file_stem"]) - set(by_file) - set(runs.values()))
    extra = {
        "id_col": ID_COL,
        "n_protein_groups_diann": int(len(pg)),
        "adapter_notes": adapted.notes,
        "n_missing_values": int(adapted.X.isna().sum().sum()),
        "n_runs_not_in_sdrf": int((meta["sdrf_source_name"].isna()).sum()) if "sdrf_source_name" in meta else len(meta),
        "sdrf_rows_corrected": sorted(fixed),
        "sdrf_files_without_run": missing_sdrf,
        "n_by_status_timepoint": meta.groupby(["lc_status", "timepoint"]).size().to_dict(),
        "n_by_arm": meta.loc[meta["lc_status"] == "LC"].groupby(["intervention_arm", "timepoint"]).size().to_dict(),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        print(s)
