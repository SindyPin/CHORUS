"""Run every proteomics cohort builder and write the generated tables of docs/proteomics_cohorts.md.

  python scripts/proteomics/build_all.py              # run the builders, then the tables
  python scripts/proteomics/build_all.py --tables-only # tables from the files already written

The tables go between the markers ``<!-- generated:start -->`` and ``<!-- generated:end -->`` of
docs/proteomics_cohorts.md; the rest of that file is written by hand. The per cohort summaries are
also written to data/proteomics/cohort_summary.csv (local).

Proteins shared: number of UniProt keys two cohorts have in common. A multi-protein key (Olink and
SomaScan assays of complexes, "P29459_P29460") counts as shared only when the same key is present.
"""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, META, ROOT  # noqa: E402

# builder module -> study ids it writes
BUILDERS = {
    "build_Buggert2025": ["Buggert2025"],
    "build_Zoodsma2022": ["Zoodsma2022_Hannover", "Zoodsma2022_Breda", "Zoodsma2022_Nijmegen"],
    "build_Cervia2024": ["Cervia2024_SomaScan", "Cervia2024_MS"],
    "build_Su2022": ["Su2022"],
    "build_Vijayakumar2022": ["Vijayakumar2022"],
    "build_Kovarik2023": ["Kovarik2023_plasma"],
    "build_Captur2022": ["Captur2022"],
    "build_Chowdhury2024": ["Chowdhury2024"],
}
DOC = ROOT / "docs" / "proteomics_cohorts.md"
START, END = "<!-- generated:start -->", "<!-- generated:end -->"


def run_builders() -> list:
    rows = []
    for mod in BUILDERS:
        path = Path(__file__).resolve().parent / f"{mod}.py"
        if not path.exists():
            print(f"[skip] {mod}: no builder")
            continue
        try:
            rows += importlib.import_module(mod).build()
        except FileNotFoundError as e:
            print(f"[skip] {mod}: {e}")
    return rows


def written_studies() -> list:
    return sorted(p.parent.name for p in DATA.glob("*/*_matrix.csv") if p.name == f"{p.parent.name}_matrix.csv")


def summarise(study: str) -> dict:
    """Summary row from the written files (the same fields as common.summary_row)."""
    X = pd.read_csv(DATA / study / f"{study}_matrix.csv", index_col=0)
    meta = pd.read_csv(META / f"{study}_meta.csv", low_memory=False)
    phase = meta.get("sample_phase", pd.Series("unknown", index=meta.index)).astype(str)
    st = meta["lc_status"].astype(str)
    post = phase != "acute"
    return {
        "study_id": study,
        "layer": meta["omics_layer"].iloc[0],
        "n_samples": len(meta),
        "n_subjects": meta["subject_id"].nunique(),
        "n_LC": int(((st == "LC") & post).sum()),
        "n_COVID_recovered": int(((st == "COVID_recovered") & post).sum()),
        "n_never_infected": int((st == "never_infected").sum()),
        "n_acute": int((phase == "acute").sum()),
        "n_other": int((~st.isin(["LC", "COVID_recovered", "never_infected"]) & post).sum()),
        "n_LC_subjects": meta.loc[(st == "LC") & post, "subject_id"].nunique(),
        "n_COVID_recovered_subjects": meta.loc[(st == "COVID_recovered") & post, "subject_id"].nunique(),
        "n_proteins": X.shape[0],
        "sample_matrix": ",".join(sorted(meta["sample_matrix"].dropna().astype(str).unique())),
        "platform": "; ".join(sorted(meta.get("platform", pd.Series(dtype=str)).dropna().astype(str).unique())),
        "proteins": set(X.index.astype(str)),
    }


def tables(other_layers: dict) -> str:
    studies = written_studies()
    rows = [summarise(s) for s in studies]
    head = "| cohort | layer | samples (subjects) | n LC | n COVID recovered | n never infected | n acute | n other or unknown | sample matrix | platform | proteins | other layer on the same people |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n"
    body = ""
    for r in rows:
        body += (
            f"| {r['study_id']} | {r['layer']} | {r['n_samples']} ({r['n_subjects']}) | {r['n_LC']} ({r['n_LC_subjects']}) | "
            f"{r['n_COVID_recovered']} ({r['n_COVID_recovered_subjects']}) | {r['n_never_infected']} | {r['n_acute']} | {r['n_other']} | "
            f"{r['sample_matrix']} | {r['platform']} | {r['n_proteins']} | {other_layers.get(r['study_id'], 'not checked')} |\n"
        )
    ov = pd.DataFrame(index=studies, columns=studies, dtype=object)
    for a in rows:
        for b in rows:
            ov.loc[a["study_id"], b["study_id"]] = len(a["proteins"] & b["proteins"])
    ovt = "| | " + " | ".join(studies) + " |\n|---|" + "---|" * len(studies) + "\n"
    for s in studies:
        ovt += f"| {s} | " + " | ".join(str(ov.loc[s, t]) for t in studies) + " |\n"
    out = pd.DataFrame([{k: v for k, v in r.items() if k != "proteins"} for r in rows])
    out.to_csv(DATA / "cohort_summary.csv", index=False)
    return (
        "### Cohorts\n\nSamples are matrix columns (runs); subjects in parentheses. n LC and n COVID recovered count post-acute and "
        "unknown-phase samples only (acute samples are in n acute whatever their later outcome).\n\n"
        + head
        + body
        + "\n### Proteins shared between cohorts\n\nNumber of UniProt keys in common (diagonal: proteins in the cohort).\n\n"
        + ovt
    )


def write_doc(text: str) -> None:
    doc = DOC.read_text()
    if START not in doc or END not in doc:
        raise ValueError(f"{DOC} lacks the generated markers")
    pre, rest = doc.split(START, 1)
    _, post = rest.split(END, 1)
    DOC.write_text(pre + START + "\n" + text + END + post)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-only", action="store_true")
    a = ap.parse_args()
    if not a.tables_only:
        for r in run_builders():
            print({k: v for k, v in r.items() if k in ("study_id", "n_samples", "n_proteins")})
    from other_layers import OTHER_LAYERS

    write_doc(tables(OTHER_LAYERS))


if __name__ == "__main__":
    main()
