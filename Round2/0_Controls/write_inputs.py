"""inputs.csv of each step folder: every data file and config the step reads, with its SHA-256.

  python Round2/0_Controls/write_inputs.py
"""

import glob
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as K  # noqa: E402

STEP01 = K.CTRL / "01_one_cohort" / "output" / "control_anchors.csv"
VERTICAL = ["Su2022", "Su2022_metabolomics", "Cervia2024_SomaScan", "Cervia2024_MS"]


def cohort_files(studies: list) -> list:
    E = K.config_entries()
    rows = []
    for s in sorted(set(studies)):
        parent = K.PARENT.get(s, s)
        for role in ("matrix", "meta"):
            rows.append((E[parent][role], f"{role} of {parent}"))
    return rows


def write(step: str, rows: list):
    rows = rows + [(str(p.relative_to(K.CODE)), "Round1 config (cohort entries)") for p in map(Path, sorted(glob.glob(str(K.CODE / "Round1/configs/*.yaml"))))]
    rows = list(dict.fromkeys(rows))
    out = pd.DataFrame([dict(path=p, role=r, sha256=K.file_sha(K.CODE / p)) for p, r in rows])
    out.to_csv(K.CTRL / step / "inputs.csv", index=False)
    print(step, len(out))


def main():
    A = pd.read_csv(STEP01).fillna("")
    cands = []
    for control in ("sex", "age"):
        cands += list(K.candidates(control)["study_id"])
    # step 01 reads the metadata of every candidate (to choose) and the anchor and partner matrices (A_all and shared features)
    meta_only = [(K.config_entries()[K.PARENT.get(s, s)]["meta"], f"meta of {K.PARENT.get(s, s)} (candidate)") for s in sorted(set(cands))]
    used = [s for s in list(A["anchor"]) + list(A["partner"]) if s]
    write("01_one_cohort", cohort_files(used) + meta_only)
    write("02_horizontal_one_plus_one", cohort_files(used) + [("Round2/0_Controls/01_one_cohort/output/control_anchors.csv", "anchors and partners (step 01)")])
    write("03_vertical_one_plus_one", cohort_files(VERTICAL))
    S = [p for s in ("01_one_cohort", "02_horizontal_one_plus_one", "03_vertical_one_plus_one")
         for p in sorted(glob.glob(str(K.CTRL / s / "output" / "raw" / "*.csv"))) + sorted(glob.glob(str(K.CTRL / s / "output" / "*.csv")))]
    rows = [(str(Path(p).relative_to(K.CODE)), "step output") for p in S]
    out = pd.DataFrame([dict(path=p, role=r, sha256=K.file_sha(K.CODE / p)) for p, r in rows])
    out.to_csv(K.CTRL / "04_summary" / "inputs.csv", index=False)
    print("04_summary", len(out))


if __name__ == "__main__":
    main()
