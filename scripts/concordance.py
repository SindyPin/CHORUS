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
  cohort_annotation.csv         per cohort: case_definition, control type (lc_status of the
                                negative class), time_anchor, sample matrix, class counts and
                                days_since_infection per class, all from the metadata

  python scripts/concordance.py configs/lc_transcriptomics_blood.yaml --out results/concordance
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


def cohort_effects(X: pd.DataFrame, M: pd.DataFrame, positive: str = "LC", negative: str = "COVID_recovered") -> pd.DataFrame:
    """Hedges g per gene (rows) and cohort (columns), cohorts in the order of ``M``."""
    y = binary_label(M, positive, negative)
    m = y.notna()
    s = M["study_id"].astype(str)
    eff = per_study_effects(X[m], y[m], s[m])
    order = [c for c in pd.unique(s) if c in eff]
    return pd.DataFrame({c: eff[c]["g"] for c in order})


def spearman_matrix(effects: pd.DataFrame) -> tuple:
    """Pairwise Spearman correlation and p value between the columns (genes with both values)."""
    cols = list(effects.columns)
    r = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols)
    p = pd.DataFrame(0.0, index=cols, columns=cols)
    for i, a in enumerate(cols):
        for b in cols[i + 1 :]:
            ok = effects[a].notna() & effects[b].notna()
            rho, pv = stats.spearmanr(effects.loc[ok, a], effects.loc[ok, b])
            r.loc[a, b] = r.loc[b, a] = rho
            p.loc[a, b] = p.loc[b, a] = pv
    return r, p


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
                "control_type": values(neg["lc_status"]),
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
    sym = pd.read_csv(a.hgnc, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna().drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
    eff.insert(0, "symbol", sym.reindex(eff.index).values)
    eff.rename_axis("gene_id").to_csv(out / "effects.csv")
    g = eff.drop(columns="symbol")
    r, p = spearman_matrix(g)
    r.to_csv(out / "spearman_all_genes.csv")
    p.to_csv(out / "spearman_all_genes_p.csv")
    top = strongest_genes(g, a.top)
    pd.DataFrame({"symbol": sym.reindex(top.index).values, "mean_abs_g": top.values}, index=top.index.rename("gene_id")).to_csv(out / f"top{a.top}_genes.csv")
    rt, pt = spearman_matrix(g.loc[top.index])
    rt.to_csv(out / f"spearman_top{a.top}.csv")
    pt.to_csv(out / f"spearman_top{a.top}_p.csv")
    ann = cohort_annotation(M, lc.positive, lc.negative)
    ann.to_csv(out / "cohort_annotation.csv", index=False)
    print(f"{g.shape[0]} common genes, {g.shape[1]} cohorts")
    print(ann.to_string(index=False))
    print("\nSpearman, all common genes\n" + r.round(3).to_string())
    print(f"\nSpearman, top {a.top} genes\n" + rt.round(3).to_string())


if __name__ == "__main__":
    main()
