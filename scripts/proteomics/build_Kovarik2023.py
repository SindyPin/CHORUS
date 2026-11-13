"""Build the Kovarik2023 untargeted plasma proteomics cohort (Kovarik JJ, Bileck A, Hagn G, ...,
Schmetterer KG, Gerner C 2023 iScience 26:105717, "A multi-omics based anti-inflammatory immune
signature characterizes long COVID-19 syndrome", doi:10.1016/j.isci.2022.105717, PMID 36507225,
PMC9719844; PRIDE PXD036969 (plasma), CC0).

Raw file (data/proteomics/Kovarik2023_plasma/raw):
  MQ_Search.zip   MaxQuant 1.6.17.0 output (timsTOF Pro, PASEF, DDA; LFQ, match between runs).
      Only MQ_Search/proteinGroups.txt, summary.txt and parameters.txt are read, directly from the
      zip (nothing is extracted).

Values: MaxQuant "LFQ intensity <experiment>" columns (39 runs, one per participant). Rows flagged
"+" in "Only identified by site", "Reverse" or "Potential contaminant" are dropped before the
adapter; the adapter also drops groups whose id contains its contaminant / decoy patterns. The id
column is "Majority protein IDs"; the adapter keys each group by its first accession (isoform
suffix removed), sets zero LFQ to missing and log2 transforms. No normalisation beyond MaxQuant
LFQ, no imputation.

Samples: MaxQuant experiment names "Covid_<code>_Protifi_Mini" with code LC<n>, i<n> or pc<n>
(13 each). The group of each prefix was established by matching the per-sample log2 LFQ columns
of the paper's Table S1 (mmc2.xlsx, columns "long-covid_k", "healthy_k", "recovered_k") to
log2 of these LFQ columns over the 278 proteins of Table S1: every Table S1 column equals exactly
one MaxQuant column (mean absolute difference 0.0; next best > 0.55). Result:
  LC<n>  -> long COVID (LC)              long-covid_1..13 = LC1..LC13
  i<n>   -> healthy vaccinated (H)       healthy_1..13   = i3,i6,i7,i8,i9,i10,i11,i12,i13,i15,i17,i19,i20
  pc<n>  -> recovered (R)                recovered_1..13 = pc1,pc2,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc14,pc15,pc16
i.e. in each group Table S1 orders the codes by their number. The supplementary file is not read
here (it is not in raw/); the mapping is hard coded in PREFIX_GROUP and the Table S1 label is kept
in ``table_s1_label``.

Clinical data: only group level (Table 1: sex f:m, median age and range, months after exposure,
symptoms); no per-sample sex, age, time since infection, severity or vaccination is deposited.
"""

from __future__ import annotations

import io
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import raw_dir, swissprot_accessions, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Kovarik2023_plasma"
PLATFORM = "LC-MS/MS DDA-PASEF timsTOF Pro (MaxQuant 1.6.17.0 LFQ, match between runs)"
ID_COL = "Majority protein IDs"
FLAG_COLS = ("Only identified by site", "Reverse", "Potential contaminant")
EXP_RE = re.compile(r"^Covid_((LC|pc|i)(\d+))_Protifi_Mini$")
# experiment code prefix -> (lc_status, control_type, sample_phase, Table S1 group label)
PREFIX_GROUP = {
    "LC": ("LC", "not_applicable", "post_acute", "long-covid"),
    "pc": ("COVID_recovered", "COVID_recovered", "post_acute", "recovered"),
    "i": ("never_infected", "never_infected_contemporary", "not_applicable", "healthy"),
}
GROUP_NOTE = {
    "LC": "authors' long COVID syndrome group (PCR-confirmed infection >=3 months before; chronic fatigue and/or lack of concentration plus >=1 further chronic symptom; anti-N and anti-S positive)",
    "pc": "authors' recovered group (infection >=3 months before, symptom free at blood draw; anti-N and anti-S positive)",
    "i": "authors' healthy group (no SARS-CoV-2 infection history, anti-N negative / anti-S positive, sampled ~3 months after two vaccine doses)",
}


def read_member(z: zipfile.ZipFile, name: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(z.read(f"MQ_Search/{name}")), sep="\t", low_memory=False)


def build() -> list:
    raw = raw_dir(STUDY)
    with zipfile.ZipFile(raw / "MQ_Search.zip") as z:
        pg = read_member(z, "proteinGroups.txt")
        summ = read_member(z, "summary.txt")
        params = read_member(z, "parameters.txt").set_index("Parameter")["Value"]

    summ = summ[summ["Raw file"] != "Total"].dropna(subset=["Experiment"])
    raw_of = summ.set_index("Experiment")["Raw file"].to_dict()

    lfq = {c: c[len("LFQ intensity "):] for c in pg.columns if c.startswith("LFQ intensity ")}
    if len(lfq) != len(raw_of) or set(lfq.values()) != set(raw_of):
        raise ValueError("LFQ columns do not match the experiments of summary.txt")

    n_groups = len(pg)
    flagged = pd.Series(False, index=pg.index)
    flag_counts = {}
    for c in FLAG_COLS:
        f = pg[c].astype(str).str.strip() == "+"
        flag_counts[c] = int(f.sum())
        flagged |= f
    table = pg.loc[~flagged, [ID_COL] + list(lfq)].rename(columns=lfq)
    adapted = to_feature_space("ms_proteomics", table, id_col=ID_COL, sample_cols=list(lfq.values()), preferred_accessions=swissprot_accessions())

    rows = []
    for exp in adapted.X.index:
        m = EXP_RE.match(exp)
        if not m:
            raise ValueError(f"unexpected experiment name {exp}")
        code, prefix, num = m.group(1), m.group(2), int(m.group(3))
        status, control, phase, s1 = PREFIX_GROUP[prefix]
        rows.append({
            "sample_id": exp, "subject_id": f"Kovarik2023_{code}", "lc_status": status, "control_type": control,
            "sample_phase": phase, "group_code": code, "group_prefix": prefix, "group_number": num, "s1_group": s1,
            "raw_file": raw_of[exp], "notes": GROUP_NOTE[prefix],
        })
    meta = pd.DataFrame(rows).set_index("sample_id")
    # Table S1 label: within a group, codes ordered by number
    meta["table_s1_label"] = meta["s1_group"] + "_" + meta.groupby("s1_group")["group_number"].rank(method="first").astype(int).astype(str)
    meta = meta.drop(columns=["s1_group", "group_number"])

    meta["case_definition"] = "study_specific"
    meta["sample_matrix"] = "plasma"
    meta["platform"] = PLATFORM
    meta["collection_country"] = "Austria"
    meta["timepoint"] = "single visit, May-July 2021"
    infected = meta["group_prefix"] != "i"
    meta["time_anchor"] = None
    meta.loc[infected, "time_anchor"] = "unknown"
    meta["acute_severity"] = None
    meta.loc[infected, "acute_severity"] = "unknown"
    meta["variant_era"] = None
    meta.loc[infected, "variant_era"] = "unknown"
    meta["mq_version"] = params.get("Version")
    meta["fasta"] = Path(str(params.get("Fasta file", "")).replace("\\", "/")).name
    meta["value_type"] = "log2 MaxQuant LFQ intensity"

    first = ["subject_id", "lc_status", "case_definition", "sample_matrix", "platform", "control_type", "sample_phase", "timepoint",
             "time_anchor", "acute_severity", "variant_era", "collection_country", "notes"]
    meta = meta[first + [c for c in meta.columns if c not in first]]

    extra = {
        "id_col": ID_COL,
        "n_protein_groups_maxquant": int(n_groups),
        "n_groups_flagged": int(flagged.sum()),
        "flag_counts": flag_counts,
        "adapter_notes": adapted.notes,
        "n_missing_values": int(adapted.X.isna().sum().sum()),
        "macrophage_cohort": "not written: PXD036970 samples are U937 cell line M1/M2 supernatants, not participants",
    }
    return [write_cohort(STUDY, adapted, meta, extra)]


if __name__ == "__main__":
    for s in build():
        print(s)
