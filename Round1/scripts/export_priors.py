"""Export the TACO 411 DCE genes and the MRCONTROL 32 causal genes as Ensembl prior lists.

Sources (clone both repositories outside this one):
  TACO                     Data/411_DCEG_Paths_average.csv: the 411 gene columns (Ensembl, versioned)
  Causal-Multiomics-Method SM/SM1/Alpha1_Summary_Statistics.tsv (16 genes at alpha = 1, Mendelian
                           randomisation) plus the 16 genes of the alpha = 0 table in README.md
                           (control theory network drivers); together the 32 putative causal genes

Identifiers are mapped through the HGNC complete set: Ensembl to symbol for the 411, symbol to
Ensembl for the 32 (approved symbol, then a previous or alias symbol that points to one gene).
Each line is an Ensembl id without version followed by "# SYMBOL".

  python scripts/export_priors.py ../TACO ../Causal-Multiomics-Method data/ref/hgnc_complete_set.txt
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd


def commit(repo: Path) -> str:
    return subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()


def symbol_to_ensembl(hgnc: pd.DataFrame, symbols: list) -> dict:
    h = hgnc.dropna(subset=["ensembl_gene_id"])
    approved = h.set_index("symbol")["ensembl_gene_id"].to_dict()
    other: dict = {}
    for col in ["prev_symbol", "alias_symbol"]:
        x = h[[col, "ensembl_gene_id"]].dropna()
        x = x.assign(**{col: x[col].str.split("|")}).explode(col)
        n = x.groupby(col)["ensembl_gene_id"].nunique()
        for s, e in x[x[col].isin(n[n == 1].index)].drop_duplicates(col).itertuples(index=False):
            other.setdefault(s, e)
    return {s: approved.get(s) or other.get(s) for s in symbols}


def write(path: Path, header: list, rows: list) -> None:
    lines = [f"# {h}" for h in header] + [f"{e}  # {s}" for e, s in rows]
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("taco")
    ap.add_argument("mrcontrol")
    ap.add_argument("hgnc")
    ap.add_argument("--out", default="priors")
    a = ap.parse_args()
    taco, mrc, out = Path(a.taco), Path(a.mrcontrol), Path(a.out)
    hgnc = pd.read_csv(a.hgnc, sep="\t", dtype=str)

    cols = pd.read_csv(taco / "Data" / "411_DCEG_Paths_average.csv", nrows=0).columns
    ens = [c.split(".")[0] for c in cols if c.startswith("ENSG")]
    e2s = hgnc.dropna(subset=["ensembl_gene_id"]).set_index("ensembl_gene_id")["symbol"].to_dict()
    missing = [e for e in ens if e not in e2s]
    write(
        out / "taco_411_dce_genes.txt",
        [
            f"TACO 411 DCE genes: gene columns of Data/411_DCEG_Paths_average.csv, SindyPin/TACO commit {commit(taco)}.",
            "Ensembl gene ids without version; symbol from the HGNC complete set (scripts/export_priors.py).",
        ],
        [(e, e2s.get(e, "no HGNC symbol")) for e in ens],
    )
    print(f"taco_411_dce_genes.txt: {len(ens)} ids ({len(set(ens))} unique); without HGNC symbol: {missing or 'none'}")

    a1 = pd.read_csv(mrc / "SM" / "SM1" / "Alpha1_Summary_Statistics.tsv", sep="\t")
    readme = (mrc / "README.md").read_text()
    a0 = re.findall(r"^\|\s+\d+\|(\S+)\s*\|", readme[readme.index("alpha=0.00"):], re.M)[:16]
    symbols = list(a1["Gene_Symbol"]) + a0
    m = symbol_to_ensembl(hgnc, symbols)
    source = dict(zip(a1["Gene_Symbol"], a1["Ensembl_Gene_ID"]))
    disagree = {s: (source[s], m[s]) for s in source if m.get(s) and m[s] != source[s]}
    unmapped = [s for s in symbols if not m.get(s)]
    write(
        out / "mrcontrol_32_genes.txt",
        [
            f"MRCONTROL 32 putative causal genes, SindyPin/Causal-Multiomics-Method commit {commit(mrc)}: 16 at alpha = 1",
            "(SM/SM1/Alpha1_Summary_Statistics.tsv, Mendelian randomisation) and 16 at alpha = 0 (README.md, control theory).",
            "Ensembl gene ids without version mapped from the symbols through the HGNC complete set (scripts/export_priors.py).",
        ],
        [(m[s], s) for s in symbols if m.get(s)],
    )
    print(f"mrcontrol_32_genes.txt: {len(symbols)} symbols ({len(set(symbols))} unique), {len(symbols) - len(unmapped)} mapped; unmapped: {unmapped or 'none'}")
    print(f"HGNC mapping disagrees with the Ensembl ids in Alpha1_Summary_Statistics.tsv: {disagree or 'none'}")


if __name__ == "__main__":
    main()
