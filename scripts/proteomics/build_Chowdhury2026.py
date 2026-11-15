"""Build the Chowdhury2026 plasma proteomics cohort (Chowdhury MMH, Quenum AJI, Rioux-Perreault C, Lucier JF,
Ilangumaran S, Piche A, Allard-Chamard H, Ramanathan S 2026 Sci Rep, "Distinct plasma proteome signature at 3 months
post-COVID-19 infection irrespective of post-COVID condition", doi:10.1038/s41598-026-46180-y, PMID 41998026,
PMC13260905; PRIDE PXD066724, released 2026-03-24).

Raw files (data/proteomics/Chowdhury2026/raw):
  diann/report.pg_matrix.tsv   DIA-NN 1.8.1 protein group matrix (173 runs; extracted from diann.7z)
  diann/report.stats.tsv       per run identification counts (kept as cohort columns)
  metadata.sdrf.tsv            PRIDE SDRF (150 rows, one per person: disease LC 53, C 62, N 35)
  PXD066724_community_annotated.sdrf.tsv   PRIDE community SDRF (same disease values for every row)

Values: DIA-NN 1.8.1 pg_matrix quantities (library-free search of UniProt UP000005640, MBR / --reanalyse,
--pg-level 1, report.log.txt). Same processing as Chowdhury2024: ``Protein.Group`` passed to the ms_proteomics
adapter with the HGNC Swiss-Prot accessions preferred, zeros missing, log2. No normalisation, imputation or filtering.

Samples: the 150 runs named in the SDRF ``comment[data file]`` (".d.7z" removed). Runs of the search that the SDRF
does not name are not used: 20 reruns of first batch samples ("sample_NNNb"), one pool ("mix_of_8") and the second
runs of samples 213 and 215 ("_NE_"; the SDRF names their "_Nc_" runs).

Labels: SDRF ``characteristics[disease]`` (paper Table 1: no infection 35, convalescent 62, PCC 53, plasma 3 months
after PCR confirmed infection, BQC19 biobank, Sherbrooke, infections before October 2021). LC = PCC (symptoms at 3
months consistent with the WHO criteria); COVID_recovered = C (PCR negative 1 month after infection, free of
infection related symptoms at 3 months); never_infected = N (PCR negative at sampling; control_type
never_infected_contemporary). Second batch file names carry a group code (Lc, c, cR, N, NE, NT, Nc); for 5 runs
(samples 226, 227, 231, 232, 233) the code is "Lc" while both SDRFs say C. The SDRF is used (its totals equal the
paper's); the code is kept in ``file_group_code`` and the 5 runs are named in the notes.

Batches: runs "Sheela_DIA_sample_NNN" (first batch) and "Sheela_sample_NNN_<code>_DIA_TRAP_400nl_flow" (second batch,
trap column at 400 nl/min) are kept in ``ms_batch``. Chowdhury2024 (PXD047623) is from the same laboratory and the same
BQC19 Sherbrooke biobank; whether people overlap cannot be resolved from the deposits.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, swissprot_accessions, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Chowdhury2026"
PLATFORM = "LC-MS/MS diaPASEF timsTOF Pro (DIA-NN 1.8.1 library-free, MaxLFQ pg_matrix)"
ID_COL = "Protein.Group"
ANNOT_COLS = ["Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"]
SEX = {"F": "female", "M": "male"}
STATUS = {"LC": "LC", "C": "COVID_recovered", "N": "never_infected"}
CONTROL = {"LC": "not_applicable", "C": "COVID_recovered", "N": "never_infected_contemporary"}
NOTE = {
    "LC": "PCC: symptoms 3 months after PCR confirmed infection consistent with the WHO criteria (paper Results)",
    "C": "convalescent: PCR negative 1 month after infection and free of infection related symptoms at 3 months",
    "N": "uninfected: SARS-CoV-2 RT-PCR negative at sample collection",
}


def stem(path: str) -> str:
    return re.sub(r"\.(d|7z)$", "", re.sub(r"\.d\.7z$", "", str(path).replace("\\", "/").split("/")[-1]))


def build() -> list:
    raw = raw_dir(STUDY)
    pg = pd.read_csv(raw / "diann" / "report.pg_matrix.tsv", sep="\t")
    stats = pd.read_csv(raw / "diann" / "report.stats.tsv", sep="\t")
    stats.index = stats["File.Name"].map(stem)
    s = pd.read_csv(raw / "metadata.sdrf.tsv", sep="\t", dtype=str)
    s["run"] = s["comment[data file]"].map(stem)
    runs = {c: stem(c) for c in pg.columns if c not in ANNOT_COLS}
    keep = {c: r for c, r in runs.items() if r in set(s["run"])}
    missing = set(s["run"]) - set(keep.values())
    if missing:
        raise ValueError(f"SDRF runs absent from the pg_matrix: {sorted(missing)}")
    table = pg[[ID_COL] + list(keep)].rename(columns=keep)
    adapted = to_feature_space("ms_proteomics", table, id_col=ID_COL, sample_cols=list(keep.values()), preferred_accessions=swissprot_accessions())

    s = s.set_index("run")
    rows = []
    for run in adapted.X.index:
        r = s.loc[run]
        d = r["characteristics[disease]"]
        code = re.search(r"sample_\d+_([A-Za-z]+)_DIA", run)
        code = code.group(1) if code else None
        note = NOTE[d] + "; BQC19 Sherbrooke, infection before October 2021"
        if code == "Lc" and d != "LC":
            note += "; file name group code Lc, SDRF disease C (SDRF used)"
        rows.append({
            "subject_id": f"{STUDY}_{r['source name']}",
            "lc_status": STATUS[d],
            "case_definition": "WHO_2021" if d == "LC" else "not_applicable",
            "control_type": CONTROL[d],
            "sample_phase": "post_acute" if d != "N" else "not_applicable",
            "sample_matrix": "plasma",
            "platform": PLATFORM,
            "timepoint": "3_months" if d != "N" else "single_visit",
            "days_since_infection": 90 if d != "N" else None,
            "time_anchor": "positive_test" if d != "N" else None,
            "acute_severity": "unknown",
            "sex": SEX.get(r["characteristics[sex]"], "unknown"),
            "collection_country": "Canada",
            "overlaps_with_study": "Chowdhury2024 (same laboratory and biobank; not resolvable)",
            "notes": note,
            "sdrf_source_name": r["source name"],
            "sdrf_individual": r["characteristics[individual]"],
            "sdrf_disease": d,
            "file_group_code": code,
            "ms_batch": "second_trap_400nl" if "TRAP" in run else "first",
            "precursors_identified": int(stats.loc[run, "Precursors.Identified"]),
            "proteins_identified": int(stats.loc[run, "Proteins.Identified"]),
            "value_type": "log2 DIA-NN 1.8.1 MaxLFQ protein group quantity (pg_matrix)",
        })
    meta = pd.DataFrame(rows, index=adapted.X.index)
    meta.index.name = "sample_id"
    extra = {
        "id_col": ID_COL,
        "n_protein_groups_diann": int(len(pg)),
        "adapter_notes": adapted.notes,
        "n_missing_values": int(adapted.X.isna().sum().sum()),
        "runs_not_in_sdrf": sorted(set(runs.values()) - set(keep.values())),
        "runs_code_Lc_sdrf_C": sorted(meta.index[(meta["file_group_code"] == "Lc") & (meta["sdrf_disease"] != "LC")]),
        "n_by_status_batch": meta.groupby(["lc_status", "ms_batch"]).size().to_dict(),
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
