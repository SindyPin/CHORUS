"""Infer sex per sample inside each cohort from XIST and RPS4Y1 and store it in the metadata.

For every enabled cohort of the config, the raw matrix is loaded (all samples, before the
one sample per subject and expression filters), converted to log2 CPM with the library size
over all genes, and passed to chorus.sex.infer_sex (two component Gaussian mixture; ambiguous
samples flagged). Outputs:

  data/meta/inferred_sex.csv         one row per sample: study_id, sample_id, inferred_sex,
                                     posterior, mahalanobis2, reason, xist_logcpm,
                                     rps4y1_logcpm, recorded sex, provenance
  data/meta/<study_id>_meta.csv      columns inferred_sex and inferred_sex_method (the
                                     provenance string) added or replaced
  <copy>                             optional copy of the per sample table (--copy)

scripts/build_cohort_meta.py merges data/meta/inferred_sex.csv back when it rebuilds the tables.

  python scripts/infer_sex.py configs/lc_transcriptomics_blood.yaml --copy results/control_sex_label/inferred_sex.csv
"""

from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import pandas as pd
import yaml

from chorus.harmonise.transcriptomics import log_cpm
from chorus.io import Cohort, _read_table
from chorus.sex import infer_sex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--meta-dir", default="data/meta")
    ap.add_argument("--copy", help="also write the per sample table here")
    ap.add_argument("--min-posterior", type=float, default=0.95)
    ap.add_argument("--outlier-p", type=float, default=0.001)
    a = ap.parse_args()
    cfg = yaml.safe_load(Path(a.config).read_text())
    today = dt.date.today().isoformat()
    tables = []
    for c in cfg["cohorts"]:
        if not c.get("enabled", True):
            continue
        co = Cohort.from_files(c["study_id"], cfg.get("layer", "bulk_rna"), c["matrix"], c["meta"], is_counts=c.get("is_counts", False))
        t = infer_sex(log_cpm(co.X), min_posterior=a.min_posterior, outlier_p=a.outlier_p)
        t["provenance"] = t["inferred_sex_method"] + f"; posterior >= {a.min_posterior}, outlier p {a.outlier_p}; scripts/infer_sex.py {today}"
        t.insert(0, "sample_id", t.index)
        t.insert(0, "study_id", c["study_id"])
        t["recorded_sex"] = co.meta["sex"].reindex(t.index)
        tables.append(t.drop(columns="inferred_sex_method"))

        meta = _read_table(c["meta"]).reset_index()
        meta["inferred_sex"] = meta["sample_id"].map(t["inferred_sex"])
        meta["inferred_sex_method"] = meta["sample_id"].map(t["provenance"])
        meta.to_csv(c["meta"], index=False)

        rec = t.dropna(subset=["recorded_sex"])
        rec = rec[rec["inferred_sex"] != "ambiguous"]
        agree = f"; agrees with recorded sex in {int((rec['inferred_sex'] == rec['recorded_sex']).sum())} of {len(rec)}" if len(rec) else ""
        print(f"{c['study_id']}: {t['inferred_sex'].value_counts().to_dict()} ({t['provenance'].iloc[0].split(';')[0]}){agree}")
    out = pd.concat(tables, ignore_index=True)
    out.to_csv(Path(a.meta_dir) / "inferred_sex.csv", index=False)
    if a.copy:
        Path(a.copy).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(a.copy, index=False)
    print(out[out["inferred_sex"] == "ambiguous"].to_string(index=False))


if __name__ == "__main__":
    main()
