"""Map a genes x samples matrix keyed by HGNC symbols to Ensembl gene identifiers.

Used for author matrices deposited with gene symbols (GSE270045), so that they share the
Ensembl feature space of the other cohorts. Symbols are matched to the approved HGNC symbol
first, then to previous symbols that point to a single approved gene. Rows that map to the
same Ensembl gene are summed (the input is on a linear scale); unmapped rows are dropped and
counted in the report.

  python scripts/map_symbols_to_ensembl.py data/counts/GSE270045/GSE270045_LC_counts.tsv.gz \
      data/ref/hgnc_complete_set.txt data/counts/GSE270045/GSE270045_counts_ensembl.tsv.gz
"""

from __future__ import annotations

import argparse

import pandas as pd


def symbol_map(hgnc: pd.DataFrame) -> pd.Series:
    hgnc = hgnc[hgnc["ensembl_gene_id"].notna()]
    approved = hgnc.set_index("symbol")["ensembl_gene_id"]
    prev = hgnc[["prev_symbol", "ensembl_gene_id"]].dropna()
    prev = prev.assign(prev_symbol=prev["prev_symbol"].str.split("|")).explode("prev_symbol")
    prev = prev[~prev["prev_symbol"].isin(approved.index)]
    unique_prev = prev.groupby("prev_symbol")["ensembl_gene_id"].nunique() == 1
    prev = prev[prev["prev_symbol"].isin(unique_prev[unique_prev].index)].drop_duplicates("prev_symbol")
    return pd.concat([approved, prev.set_index("prev_symbol")["ensembl_gene_id"]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("hgnc")
    ap.add_argument("out")
    a = ap.parse_args()
    X = pd.read_csv(a.matrix, sep="\t", index_col=0)
    m = symbol_map(pd.read_csv(a.hgnc, sep="\t", dtype=str))
    ens = X.index.to_series().map(m)
    print(f"{ens.notna().sum()} of {len(X)} symbols mapped to Ensembl; {ens.isna().sum()} dropped")
    Y = X[ens.notna()].groupby(ens[ens.notna()]).sum()
    Y.index.name = "gene_id"
    print(f"{len(Y)} Ensembl genes x {Y.shape[1]} samples written to {a.out}")
    Y.to_csv(a.out, sep="\t")


if __name__ == "__main__":
    main()
