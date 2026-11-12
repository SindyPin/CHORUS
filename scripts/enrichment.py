"""Enrichment of each LOSO fold's invariant gene set in the prior gene lists.

For every key ``<held_out>|<harmoniser>|invariant`` of selected_features.json, the selected
genes are tested against each list with a one sided hypergeometric test (overlap at least as
large as observed) over the expressed universe: the genes that pass the expression filter of
the config in every cohort, which is the set the selector chose from. List genes outside the
universe are dropped before the test. Lists:

  TACO_411        priors/taco_411_dce_genes.txt (Ensembl)
  MRCONTROL_32    priors/mrcontrol_32_genes.txt (Ensembl)
  review_complement  priors/review_complement_9.txt (UniProt accessions, mapped to Ensembl
                  through the uniprot_ids column of the HGNC complete set)

q_bh is the Benjamini-Hochberg adjustment over all rows of the output.

  python scripts/enrichment.py configs/lc_transcriptomics_blood.yaml \
      results/lc_transcriptomics_blood/selected_features.json --out results/enrichment.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
from scipy import stats

from chorus.cli import load_config
from chorus.features.priors import load_gene_list

HGNC = "data/ref/hgnc_complete_set.txt"
LISTS = {
    "TACO_411": "priors/taco_411_dce_genes.txt",
    "MRCONTROL_32": "priors/mrcontrol_32_genes.txt",
    "review_complement": "priors/review_complement_9.txt",
}


def uniprot_to_ensembl(accessions: Iterable[str], hgnc: pd.DataFrame) -> List[str]:
    """Ensembl gene ids of the HGNC genes whose uniprot_ids (pipe separated) contain an accession."""
    acc = set(accessions)
    out = []
    for _, r in hgnc.dropna(subset=["uniprot_ids", "ensembl_gene_id"]).iterrows():
        if acc & set(str(r["uniprot_ids"]).split("|")):
            out.append(r["ensembl_gene_id"])
    return sorted(set(out))


def hypergeom_enrichment(selected: Iterable[str], gene_list: Iterable[str], universe: Iterable[str]) -> Dict[str, float]:
    """One sided hypergeometric test of the overlap between ``selected`` and ``gene_list`` within ``universe``."""
    U = set(universe)
    S = set(selected) & U
    L = set(gene_list) & U
    k = len(S & L)
    expected = len(S) * len(L) / len(U) if U else np.nan
    p = float(stats.hypergeom.sf(k - 1, len(U), len(L), len(S))) if U else np.nan
    return {
        "universe": len(U),
        "list_in_universe": len(L),
        "n_selected": len(S),
        "overlap": k,
        "expected": expected,
        "fold_enrichment": k / expected if expected else np.nan,
        "p_hypergeom": p,
    }


def bh(p: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjusted p values (NaN kept)."""
    ok = p.dropna()
    n = len(ok)
    order = ok.sort_values()
    q = (order * n / np.arange(1, n + 1))[::-1].cummin()[::-1].clip(upper=1.0)
    return q.reindex(p.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("selected")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hgnc", default=HGNC)
    ap.add_argument("--selector", default="invariant")
    a = ap.parse_args()
    _, _, cfg = load_config(a.config)
    universe = cfg["filtered_features"]
    hgnc = pd.read_csv(a.hgnc, sep="\t", usecols=["symbol", "ensembl_gene_id", "uniprot_ids"], dtype=str)
    sym = hgnc.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]
    lists = {}
    for name, path in LISTS.items():
        genes = load_gene_list(path)
        if not all(g.startswith("ENSG") for g in genes):
            genes = uniprot_to_ensembl(genes, hgnc)
        lists[name] = genes
    sel = json.loads(Path(a.selected).read_text())
    rows = []
    for key, genes in sel.items():
        held_out, harmoniser, selector = key.split("|")
        if selector != a.selector:
            continue
        for name, lg in lists.items():
            r = hypergeom_enrichment(genes, lg, universe)
            ov = sorted(set(genes) & set(lg) & set(universe))
            rows.append({"held_out": held_out, "harmoniser": harmoniser, "selector": selector, "gene_list": name, "list_size": len(lg), **r, "overlap_genes": ";".join(str(sym.get(g, g)) for g in ov)})
    res = pd.DataFrame(rows)
    res.insert(res.columns.get_loc("p_hypergeom") + 1, "q_bh", bh(res["p_hypergeom"]))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out, index=False)
    print(res.drop(columns="overlap_genes").to_string(index=False))


if __name__ == "__main__":
    main()
