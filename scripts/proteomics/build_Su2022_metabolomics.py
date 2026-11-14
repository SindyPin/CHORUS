"""Build the Su2022 plasma metabolomics layer (Su Y, ..., Heath JR 2022 Cell, INCOV; Mendeley Data
10.17632/96v329bg7g.1, Table S2.xlsx sheet "S2.2 Metabolomics", CC BY 4.0), the second layer measured on
the same blood draws as the Olink cohort Su2022.

Values: the authors' "batch-corrected, and age-, sex-, BMI-adjusted Z-scores plasma metabolomics data"
(Table S2.0 catalogue), one row per draw, one column per metabolite name; passed through the metabolomics
adapter unchanged (names as deposited, no RefMet mapping, no transform).

Sample ids follow build_Su2022.py: "HC-<id>" for healthy donors, the Blood Draw (INCOV<nnn>-T<k>) for patients,
so a draw has the same sample_id in both layers. The metadata are built with build_Su2022.sample_metadata from
the same clinical tables. Writes data/metabolomics/Su2022_metabolomics/ and data/meta/Su2022_metabolomics_meta.csv.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_Su2022 import read_xlsx_sheet, sample_metadata  # noqa: E402
from common import ROOT, raw_dir, write_cohort  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402

STUDY = "Su2022_metabolomics"
PLATFORM = "untargeted plasma metabolomics (authors' batch corrected, age, sex and BMI adjusted Z scores)"
OUT = ROOT / "data" / "metabolomics"


def build() -> list:
    raw = raw_dir("Su2022")
    M = read_xlsx_sheet(raw / "Table S2.xlsx", "S2.2 Metabolomics")
    id_cols = ["Patient Subject ID", "Blood Draw", "Healthy or INCOV", "age", "sex", "BMI"]
    analytes = [c for c in M.columns if c not in id_cols]
    healthy = M["Healthy or INCOV"] == "Healthy"
    hid = M["Patient Subject ID"].map(lambda x: f"{int(x)}" if isinstance(x, float) else str(x))
    M = M.assign(sample_id=np.where(healthy, "HC-" + hid, M["Blood Draw"].astype(str)), subject_id=np.where(healthy, "HC-" + hid, M["Patient Subject ID"].astype(str)))
    if M["sample_id"].duplicated().any():
        raise ValueError("duplicated sample ids in S2.2")
    adapted = to_feature_space("metabolomics", M.set_index("sample_id")[analytes])
    meta = sample_metadata(M, raw, PLATFORM).loc[adapted.X.index]
    extra = {"n_analytes": len(analytes), "values_missing": int(adapted.X.isna().sum().sum()), "n_timepoint": ",".join(f"{k}={v}" for k, v in meta["timepoint"].value_counts().sort_index().items())}
    return [write_cohort(STUDY, adapted, meta, extra, root=OUT)]


if __name__ == "__main__":
    for s in build():
        for k, v in s.items():
            print(f"{k}: {v}")
