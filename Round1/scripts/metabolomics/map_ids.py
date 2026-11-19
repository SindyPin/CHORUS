"""Map the metabolite names of every metabolomics cohort to the shared identifier
(chorus.metabolites, adapter metabolomics_id) and write the mapping report.

Per cohort, in data/metabolomics/<study>/ (local):
  <study>_id_map.csv         every deposited name: query, RefMet name and ids, shared_id, id_source
  <study>_ids_matrix.csv     shared identifiers (rows) x samples, the matrix of configs/lc_metabolomics*.yaml
  <study>_ids_features.csv   per identifier, the names that mapped to it and the one kept

Report (committed), in results/lc_metabolomics/mapping/:
  mapping_summary.csv        per cohort: names, unannotated, matched to RefMet, with HMDB, mapped
                             (by id_source), identifiers after duplicates are dropped
  n_shared.csv               cohort x cohort number of shared identifiers
  shared_by_k.csv            number of identifiers present in at least k cohorts, and the
                             triples of cohorts (among the labelled cohorts) with their shared counts

RefMet match answers (name conversion service, batch POST) are cached in data/ref/refmet/refmet_match_cache.csv; the RefMet table is
data/ref/refmet/refmet.csv (docs/provenance.md).

  python scripts/metabolomics/map_ids.py [--offline]
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "proteomics"))
from common import ROOT  # noqa: E402

from chorus.adapters import to_feature_space  # noqa: E402
from chorus.metabolites import UNANNOTATED, fetch_match, load_refmet, match_names, metabolite_ids, pairwise_shared, strip_prefix  # noqa: E402

DATA = ROOT / "data" / "metabolomics"
REFMET = ROOT / "data" / "ref" / "refmet" / "refmet.csv"
CACHE = ROOT / "data" / "ref" / "refmet" / "refmet_match_cache.csv"
STUDIES = ["Martinez2024", "Lopez2023", "Saito2024", "Holmes2021", "Su2022_metabolomics"]
LABELLED = ["Martinez2024", "Lopez2023", "Saito2024", "Su2022_metabolomics"]  # with LC and COVID_recovered samples


def deposited_ids(study: str):
    p = DATA / study / f"{study}_deposited_ids.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p, dtype=str).set_index("name")
    d["hmdb_id"] = d["external_identifier"].where(d["external_identifier"].str.startswith("HMDB", na=False))
    return d


def map_study(study: str, refmet: pd.DataFrame, offline: bool = False) -> dict:
    X = pd.read_csv(DATA / study / f"{study}_matrix.csv", index_col=0).T
    names = list(X.columns)
    queries = [strip_prefix(n) for n in names if not UNANNOTATED.match(strip_prefix(n))]
    matches = match_names(queries, CACHE, fetch=None if offline else fetch_match)
    T = metabolite_ids(names, refmet, matches, deposited_ids(study))
    T.to_csv(DATA / study / f"{study}_id_map.csv", index=False)
    ad = to_feature_space("metabolomics_id", X, ids=T.set_index("name")["shared_id"], log2=False)
    M = ad.X.T
    M.index.name = "metabolite_id"
    M.to_csv(DATA / study / f"{study}_ids_matrix.csv")
    ad.features.to_csv(DATA / study / f"{study}_ids_features.csv", index=False)
    src = T["id_source"].value_counts()
    row = {
        "study_id": study,
        "n_names": len(names),
        "n_unannotated": int(src.get("unannotated", 0)),
        "n_matched_refmet": int(T["refmet_name"].notna().sum()),
        "n_with_hmdb": int(T["shared_id"].fillna("").str.startswith("HMDB:").sum()),
        **{f"n_{k}": int(src.get(k, 0)) for k in ("hmdb_refmet", "hmdb_deposited", "refmet", "inchikey_deposited", "unmatched")},
        "n_names_mapped": int(T["shared_id"].notna().sum()),
        "n_identifiers": M.shape[0],
        "n_names_dropped_as_duplicate": int(ad.notes["names_dropped_as_duplicate"]),
    }
    return {"row": row, "ids": list(M.index)}


DOC = ROOT / "docs" / "metabolomics_cohorts.md"
START, END = "<!-- generated:start -->", "<!-- generated:end -->"


def cohort_rows() -> pd.DataFrame:
    rows = []
    for s in STUDIES:
        p = ROOT / "data" / "meta" / f"{s}_meta.csv"
        if not p.exists():
            continue
        m = pd.read_csv(p, low_memory=False)
        st, ph = m["lc_status"].astype(str), m.get("sample_phase", pd.Series("unknown", index=m.index)).astype(str)
        post = ph != "acute"
        rows.append({"cohort": s, "samples (subjects)": f"{len(m)} ({m['subject_id'].nunique()})", "LC": int(((st == "LC") & post).sum()), "COVID_recovered": int(((st == "COVID_recovered") & post).sum()),
                     "never_infected": int((st == "never_infected").sum()), "acute": int((ph == "acute").sum()), "label unknown": int(((st == "unknown") & post).sum()),
                     "names": len(pd.read_csv(DATA / s / f"{s}_matrix.csv", usecols=[0]))})
    return pd.DataFrame(rows)


def md(df: pd.DataFrame, index: bool = False) -> str:
    d = df.reset_index() if index else df
    return "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * d.shape[1] + "\n" + "".join("| " + " | ".join(map(str, r)) + " |\n" for r in d.itertuples(index=False))


def write_doc(summ: pd.DataFrame, shared: pd.DataFrame, by_k: pd.DataFrame) -> None:
    """The generated tables of docs/metabolomics_cohorts.md (between the generated markers)."""
    if not DOC.exists() or START not in DOC.read_text():
        return
    text = (
        "### Cohorts\n\nSamples per class after the builders (post-acute samples for LC and COVID_recovered; acute samples are counted under acute whatever their later outcome).\n\n"
        + md(cohort_rows())
        + "\n### Mapping per cohort\n\nNames: metabolite columns as deposited. unannotated: m/z and RT features, not sent to RefMet. matched_refmet: names standardised by RefMet. "
        "hmdb_refmet, hmdb_deposited, refmet, inchikey_deposited: names mapped by each rule. identifiers: after duplicate names of one identifier are dropped.\n\n"
        + md(summ)
        + "\n### Identifiers shared between cohorts\n\nDiagonal: identifiers in the cohort.\n\n"
        + md(shared, index=True).replace("| index |", "| |", 1)
        + "\n### Identifiers shared by several cohorts\n\n"
        + md(by_k)
    )
    doc = DOC.read_text()
    pre, rest = doc.split(START, 1)
    _, post = rest.split(END, 1)
    DOC.write_text(pre + START + "\n" + text + END + post)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="use the match cache only")
    ap.add_argument("--out", default=str(ROOT / "results" / "lc_metabolomics" / "mapping"))
    a = ap.parse_args()
    refmet = load_refmet(REFMET)
    res = {s: map_study(s, refmet, a.offline) for s in STUDIES if (DATA / s / f"{s}_matrix.csv").exists()}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    summ = pd.DataFrame([r["row"] for r in res.values()])
    summ.to_csv(out / "mapping_summary.csv", index=False)
    ids = {s: r["ids"] for s, r in res.items()}
    pairwise_shared(ids).to_csv(out / "n_shared.csv")
    rows = []
    allids = pd.Series([i for v in ids.values() for i in set(v)]).value_counts()
    for k in range(1, len(ids) + 1):
        rows.append({"what": f"identifiers in at least {k} cohorts", "cohorts": "all", "n": int((allids >= k).sum())})
    lab = [s for s in LABELLED if s in ids]
    for r in range(3, len(lab) + 1):
        for combo in itertools.combinations(lab, r):
            n = len(set.intersection(*(set(ids[c]) for c in combo)))
            rows.append({"what": f"identifiers shared by all {r} cohorts", "cohorts": ";".join(combo), "n": n})
    pd.DataFrame(rows).to_csv(out / "shared_by_k.csv", index=False)
    write_doc(summ, pairwise_shared(ids), pd.DataFrame(rows))
    print(summ.to_string(index=False))
    print(pairwise_shared(ids).to_string())
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
