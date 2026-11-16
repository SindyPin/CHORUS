"""Concordance of the LC effect across cohorts.

Inside each cohort the standardised LC effect of every gene is Hedges g (LC minus
COVID_recovered, chorus.features.per_study_effects) on log CPM, after the preprocessing of the
config (one sample per subject, deduplication) and on the common expressed genes (the genes
that pass the expression filter in every cohort). Outputs, in --out:

  effects.csv                   gene x cohort Hedges g, with the gene symbol
  spearman_all_genes.csv        cohort x cohort Spearman correlation of g over all common genes
  spearman_all_genes_p.csv      the matching two sided p values
  spearman_top500.csv           the same over the 500 strongest genes (largest mean |g| across
                                cohorts), and spearman_top500_p.csv
  top500_genes.csv              those genes, with mean |g|
  cohort_annotation.csv         per cohort: case_definition, negative_class_lc_status (lc_status of
                                the negative class), control_type (the schema control_type of the
                                negative class, for example convalescent_unverified), time_anchor,
                                sample matrix, class counts and days_since_infection per class, all
                                from the metadata
  n_shared.csv                  cohort x cohort number of features with an effect in both
  n_shared_top500.csv           the same among the 500 strongest genes
  pairs.csv                     one row per pair of cohorts and feature set (all, top500): n_shared,
                                spearman, p and fewer_than_10_shared (n_shared below MIN_SHARED)

In spearman_top500.csv (and its p values) a pair that shares fewer than MIN_SHARED (10) of the
strongest features is NaN. spearman_all_genes.csv keeps every pair; pairs with fewer than 10 shared
features are flagged in pairs.csv.

Metabolomics configs (feature_space metabolite_id) are keyed by the shared metabolite identifier and
labelled with the RefMet name. Proteomics configs (feature_space uniprot) are keyed by UniProt accession; the symbol comes from
the uniprot_ids column of the HGNC complete set (a multi-protein key "P1_P2" gets the symbols
joined by "_"). With ``preprocessing: features: union`` in the config, every cohort keeps its own
features and each pair of cohorts is correlated over the features both have (n_shared.csv). The
values are those after the config's preprocessing, including its imputation.

  python scripts/concordance.py configs/lc_transcriptomics_blood.yaml --out results/concordance
  python scripts/concordance.py configs/lc_proteomics_affinity.yaml --out results/lc_proteomics_affinity/concordance --top 50
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from chorus.cli import load_config
from chorus.features import per_study_effects
from chorus.schema import binary_label

HGNC = "data/ref/hgnc_complete_set.txt"
MIN_SHARED = 10  # fewest shared features for a correlation over the strongest features; fewer are flagged everywhere
REFMET = "data/ref/refmet/refmet.csv"


def cohort_effects(X: pd.DataFrame, M: pd.DataFrame, positive: str = "LC", negative: str = "COVID_recovered") -> pd.DataFrame:
    """Hedges g per gene (rows) and cohort (columns), cohorts in the order of ``M``."""
    y = binary_label(M, positive, negative)
    m = y.notna()
    s = M["study_id"].astype(str)
    eff = per_study_effects(X[m], y[m], s[m])
    order = [c for c in pd.unique(s) if c in eff]
    return pd.DataFrame({c: eff[c]["g"] for c in order})


def spearman_matrix(effects: pd.DataFrame, min_shared: int = 0) -> tuple:
    """Pairwise Spearman correlation and p value between the columns (genes with both values).
    Pairs with fewer than ``min_shared`` genes with both values are NaN."""
    cols = list(effects.columns)
    r = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols)
    p = pd.DataFrame(0.0, index=cols, columns=cols)
    for i, a in enumerate(cols):
        for b in cols[i + 1 :]:
            ok = effects[a].notna() & effects[b].notna()
            if ok.sum() < max(min_shared, 1):
                r.loc[a, b] = r.loc[b, a] = p.loc[a, b] = p.loc[b, a] = np.nan
                continue
            rho, pv = stats.spearmanr(effects.loc[ok, a], effects.loc[ok, b])
            r.loc[a, b] = r.loc[b, a] = rho
            p.loc[a, b] = p.loc[b, a] = pv
    return r, p


def shared_counts(effects: pd.DataFrame) -> pd.DataFrame:
    """Number of features with a value in both columns, for every pair of columns."""
    ok = effects.notna().astype(int)
    return ok.T @ ok


def pair_table(r: pd.DataFrame, p: pd.DataFrame, n: pd.DataFrame, feature_set: str, min_shared: int = MIN_SHARED) -> pd.DataFrame:
    """One row per pair of columns: shared features, correlation, p value and the flag for fewer
    than ``min_shared`` shared features."""
    cols = list(r.columns)
    rows = [
        {"feature_set": feature_set, "cohort_a": a, "cohort_b": b, "n_shared": int(n.loc[a, b]), "spearman": r.loc[a, b], "p": p.loc[a, b],
         f"fewer_than_{min_shared}_shared": bool(n.loc[a, b] < min_shared)}
        for i, a in enumerate(cols) for b in cols[i + 1 :]
    ]
    return pd.DataFrame(rows)


def feature_symbols(features: pd.Index, hgnc: str, feature_space: str = "ensembl_gene") -> pd.Series:
    """HGNC symbol per feature: by Ensembl gene id, or by UniProt accession for ``uniprot``. For
    ``metabolite_id`` (``HMDB:<id>`` or ``RefMet:<name>``) the RefMet name, from the RefMet table."""
    if feature_space == "metabolite_id":
        t = pd.read_csv(REFMET, dtype=str, usecols=["refmet_name", "hmdb_id"]).dropna()
        by_hmdb = t.drop_duplicates("hmdb_id").set_index("hmdb_id")["refmet_name"]
        out = [str(f).split(":", 1)[1] if str(f).startswith("RefMet:") else by_hmdb.get(str(f).split(":", 1)[1]) for f in features]
        return pd.Series(out, index=features, dtype=object)
    if feature_space != "uniprot":
        t = pd.read_csv(hgnc, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
        return t.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"].reindex(features)
    t = pd.read_csv(hgnc, sep="\t", usecols=["symbol", "uniprot_ids"], dtype=str).dropna()
    acc = {}
    for sym, ids in zip(t["symbol"], t["uniprot_ids"]):
        for a in ids.split("|"):
            acc.setdefault(a.strip(), sym)
    out = [("_".join(acc.get(a, "") for a in str(f).split("_")) if all(a in acc for a in str(f).split("_")) else None) for f in features]
    return pd.Series(out, index=features, dtype=object)


def strongest_genes(effects: pd.DataFrame, n: int = 500) -> pd.Series:
    """The ``n`` genes with the largest mean absolute effect across cohorts."""
    return effects.abs().mean(axis=1).sort_values(ascending=False).head(n)


def cohort_annotation(M: pd.DataFrame, positive: str = "LC", negative: str = "COVID_recovered") -> pd.DataFrame:
    rows = []

    def values(x):
        return ";".join(sorted(x.dropna().astype(str).unique())) or "missing"

    for sid, m in M.groupby("study_id", sort=False):
        pos, neg = m[m["lc_status"] == positive], m[m["lc_status"] == negative]
        d_pos = pd.to_numeric(pos.get("days_since_infection"), errors="coerce")
        d_neg = pd.to_numeric(neg.get("days_since_infection"), errors="coerce")
        rows.append(
            {
                "study_id": sid,
                "case_definition": values(pos["case_definition"]),
                "negative_class_lc_status": values(neg["lc_status"]),
                "control_type": values(neg["control_type"]) if "control_type" in neg else "missing",
                "control_case_definition": values(neg["case_definition"]),
                "time_anchor": values(m.loc[m["lc_status"].isin([positive, negative]), "time_anchor"]) if "time_anchor" in m else "missing",
                "sample_matrix": values(m["sample_matrix"]),
                "n_LC": len(pos),
                "n_control": len(neg),
                "other_lc_status": values(m.loc[~m["lc_status"].isin([positive, negative]), "lc_status"]),
                "days_LC_median": d_pos.median(),
                "days_LC_min": d_pos.min(),
                "days_LC_max": d_pos.max(),
                "days_control_median": d_neg.median(),
                "days_control_min": d_neg.min(),
                "days_control_max": d_neg.max(),
            }
        )
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=500)
    ap.add_argument("--hgnc", default=HGNC)
    a = ap.parse_args()
    cs, lc, cfg = load_config(a.config)
    feats = cfg["filtered_features"]
    X, M = cs.pooled(feats)
    eff = cohort_effects(X, M, lc.positive, lc.negative)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    space = cs.cohorts[0].feature_space
    sym = feature_symbols(eff.index, a.hgnc, space)
    eff.insert(0, "symbol", sym.values)
    eff.rename_axis("gene_id").to_csv(out / "effects.csv")
    g = eff.drop(columns="symbol")
    r, p = spearman_matrix(g)
    r.to_csv(out / "spearman_all_genes.csv")
    p.to_csv(out / "spearman_all_genes_p.csv")
    n = shared_counts(g)
    n.to_csv(out / "n_shared.csv")
    top = strongest_genes(g, a.top)
    pd.DataFrame({"symbol": sym.reindex(top.index).values, "mean_abs_g": top.values}, index=top.index.rename("gene_id")).to_csv(out / f"top{a.top}_genes.csv")
    rt, pt = spearman_matrix(g.loc[top.index], MIN_SHARED)
    nt = shared_counts(g.loc[top.index])
    rt.to_csv(out / f"spearman_top{a.top}.csv")
    pt.to_csv(out / f"spearman_top{a.top}_p.csv")
    nt.to_csv(out / f"n_shared_top{a.top}.csv")
    pairs = pd.concat([pair_table(r, p, n, "all"), pair_table(rt, pt, nt, f"top{a.top}")], ignore_index=True)
    pairs.to_csv(out / "pairs.csv", index=False)
    ann = cohort_annotation(M, lc.positive, lc.negative)
    ann.to_csv(out / "cohort_annotation.csv", index=False)
    print(f"{g.shape[0]} common genes, {g.shape[1]} cohorts")
    print(ann.to_string(index=False))
    print("\nSpearman, all common genes\n" + r.round(3).to_string())
    print(f"\nSpearman, top {a.top} genes\n" + rt.round(3).to_string())
    low = pairs[pairs[f"fewer_than_{MIN_SHARED}_shared"]]
    print(f"\nPairs with fewer than {MIN_SHARED} shared features: {len(low)}\n" + (low.to_string(index=False) if len(low) else ""))


if __name__ == "__main__":
    main()
