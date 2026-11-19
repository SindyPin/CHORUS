"""External causal evidence lists used as feature priors.

CHORUS ships no gene lists of its own. Point the loaders at the published supplementary
tables of the group's prior work, or any other list:

  MRCONTROL  32 putative causal genes (TWMR + control theory), Pinero et al. 2025 PLoS Comput Biol
  TACO       411 DCE genes and the 6 gene nested CV signature (TYMP, MAPK8, PLA2G4B, PIM2, TK1, STAT3),
             Pinero et al. 2026 Mach Learn Sci Technol
  SPLIT      consensus driver genes supported by two or more of TWMR, CT and DCE

Identifiers are harmonised to the matrix feature space (Ensembl gene without version, or
HGNC symbol) through an optional mapping table.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


def load_gene_list(path: str | Path, column: Optional[str] = None, strip_version: bool = True) -> List[str]:
    """Load a one column text file or a csv/tsv/xlsx column of gene identifiers."""
    path = Path(path)
    if path.suffix in {".txt", ".lst"}:
        genes = [ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()]  # "#" starts a comment, also after an id
        genes = [g for g in genes if g]
    else:
        if path.suffix in {".xlsx", ".xls"}:
            df = pd.read_excel(path)
        else:
            df = pd.read_csv(path, sep="\t" if path.suffix in {".tsv", ".txt"} else ",")
        col = column or df.columns[0]
        genes = df[col].dropna().astype(str).tolist()
    if strip_version:
        genes = [g.split(".")[0] if g.startswith("ENSG") else g for g in genes]
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


@dataclass
class CausalPriors:
    lists: Dict[str, List[str]] = field(default_factory=dict)

    def add(self, name: str, genes: List[str]) -> "CausalPriors":
        self.lists[name] = list(genes)
        return self

    def union(self, *names: str) -> List[str]:
        names = names or tuple(self.lists)
        out, seen = [], set()
        for n in names:
            for g in self.lists[n]:
                if g not in seen:
                    seen.add(g)
                    out.append(g)
        return out

    def intersection(self, *names: str) -> List[str]:
        names = names or tuple(self.lists)
        sets = [set(self.lists[n]) for n in names]
        common = set.intersection(*sets) if sets else set()
        return [g for g in self.lists[names[0]] if g in common]

    def evidence_count(self) -> pd.Series:
        """Number of lists supporting each gene (a crude evidence score)."""
        counts: Dict[str, int] = {}
        for genes in self.lists.values():
            for g in set(genes):
                counts[g] = counts.get(g, 0) + 1
        return pd.Series(counts).sort_values(ascending=False)

    def to_feature_space(self, mapping: pd.DataFrame, src_col: str, dst_col: str) -> "CausalPriors":
        m = dict(zip(mapping[src_col].astype(str), mapping[dst_col].astype(str)))
        return CausalPriors({k: [m[g] for g in v if g in m] for k, v in self.lists.items()})
