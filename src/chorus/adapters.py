"""Layer adapters: one declared function per omics layer that turns a raw matrix into the
shared feature space of that layer.

Every adapter returns an ``Adapted`` record: a samples by features matrix keyed by the shared
identifier, a feature table that says which source identifiers went into each feature, and
the name of the feature space. The shared spaces are

  bulk_rna        Ensembl gene id without version (``ENSG00000123456``); duplicates after the
                  version is removed are summed (counts are additive)
  olink           UniProt accession; an assay measured on several panels is averaged (NPX is
                  log2); a multi-protein assay is keyed by its accessions sorted and joined by "_"
  somascan        UniProt accession; human protein aptamers only, log2 of RFU, several aptamers of
                  one protein averaged on the log2 scale
  ms_proteomics   UniProt accession of the leading protein of each protein group, isoform suffix
                  removed; contaminant and decoy groups dropped; zero intensities missing; log2;
                  groups with the same leading accession averaged on the log2 scale

Adapters only rename, filter and aggregate. Normalisation and cross study harmonisation stay in
``chorus.harmonise``. Register a new layer with ``@layer_adapter("name", "feature_space")`` and
call it through ``to_feature_space(layer, raw, ...)``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, Optional

import numpy as np
import pandas as pd

UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
# above this maximum a matrix is on the linear scale and is log2 transformed
LINEAR_SCALE_MAX = 50.0


@dataclass
class Adapted:
    X: pd.DataFrame  # samples x features, shared identifiers as columns
    features: pd.DataFrame  # one row per output feature: feature, source_ids (joined by " | "), n_sources
    feature_space: str
    layer: str
    notes: Dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class LayerAdapter:
    layer: str
    feature_space: str
    fn: Callable[..., Adapted]


ADAPTERS: Dict[str, LayerAdapter] = {}


def layer_adapter(layer: str, feature_space: str):
    """Register ``fn`` as the adapter of ``layer``. One adapter per layer."""

    def deco(fn: Callable[..., Adapted]) -> Callable[..., Adapted]:
        if layer in ADAPTERS:
            raise ValueError(f"an adapter for layer {layer} is already registered")
        ADAPTERS[layer] = LayerAdapter(layer, feature_space, fn)
        return fn

    return deco


def to_feature_space(layer: str, raw, **kwargs) -> Adapted:
    """Turn a raw matrix of ``layer`` into its shared feature space."""
    if layer not in ADAPTERS:
        raise KeyError(f"no adapter for layer {layer}; registered: {sorted(ADAPTERS)}")
    out = ADAPTERS[layer].fn(raw, **kwargs)
    if out.X.columns.duplicated().any():
        raise AssertionError(f"{layer} adapter returned duplicated features")
    return out


# ---------------------------------------------------------------------------------- helpers
def canonical_uniprot(acc: str) -> Optional[str]:
    """UniProt accession without isoform suffix, or None when the string is not one."""
    a = str(acc).strip().split("-")[0]
    return a if UNIPROT_RE.match(a) else None


def uniprot_key(accessions: Iterable[str]) -> Optional[str]:
    """Key of an assay or aptamer with one or several targets: canonical accessions, sorted,
    de-duplicated and joined by "_". None when no accession is valid."""
    accs = sorted({a for a in (canonical_uniprot(x) for x in accessions) if a})
    return "_".join(accs) if accs else None


def split_ids(s, seps: str = r"[;,|_\s]+") -> list:
    return [x for x in re.split(seps, str(s)) if x and x.lower() != "nan"]


def _aggregate(wide: pd.DataFrame, key: pd.Series, how: str = "mean") -> tuple:
    """Columns of ``wide`` renamed by ``key`` (source id -> shared id, NaN dropped) and
    aggregated per shared id. Returns (matrix, feature table)."""
    key = key.reindex(wide.columns).dropna()
    W = wide[key.index]
    grouped = W.T.groupby(key.values, sort=True)
    X = (grouped.sum(min_count=1) if how == "sum" else grouped.mean()).T
    feats = key.groupby(key.values).apply(lambda s: " | ".join(map(str, s.index))).rename("source_ids").to_frame()
    feats["n_sources"] = key.value_counts().reindex(feats.index).values
    feats.index.name = "feature"
    return X, feats.reset_index()


def _log2_if_linear(X: pd.DataFrame) -> tuple:
    if np.nanmax(X.values) > LINEAR_SCALE_MAX:
        return np.log2(X.where(X > 0)), "log2"
    return X, "as_provided"


# ---------------------------------------------------------------------------------- adapters
@layer_adapter("bulk_rna", "ensembl_gene")
def bulk_rna(raw: pd.DataFrame) -> Adapted:
    """Samples x genes (counts or continuous). Ensembl version suffixes are removed and the
    columns that then coincide are summed; other identifiers are kept as they are."""
    cols = pd.Series([c.split(".")[0] if str(c).startswith("ENSG") else c for c in raw.columns], index=raw.columns)
    if cols.duplicated().any():
        X, feats = _aggregate(raw.apply(pd.to_numeric, errors="coerce"), cols, how="sum")
    else:
        X = raw.copy()
        X.columns = cols.values
        feats = pd.DataFrame({"feature": cols.values, "source_ids": [str(c) for c in raw.columns], "n_sources": 1})
    return Adapted(X, feats, "ensembl_gene", "bulk_rna")


@layer_adapter("olink", "uniprot")
def olink(
    raw: pd.DataFrame,
    sample_col: str = "SampleID",
    assay_col: str = "OlinkID",
    uniprot_col: str = "UniProt",
    value_col: str = "NPX",
    qc_col: Optional[str] = None,
    qc_fail: Iterable[str] = ("FAIL", "WARN", "Warning"),
    keep_samples: Optional[Iterable[str]] = None,
    annotation: Optional[pd.DataFrame] = None,
) -> Adapted:
    """Olink NPX into UniProt space.

    ``raw`` is either the long NPX export (one row per sample and assay, with ``sample_col``,
    ``assay_col``, ``uniprot_col`` and ``value_col``), or a wide samples x assays matrix with an
    ``annotation`` table mapping ``assay_col`` to ``uniprot_col``. Values whose ``qc_col`` is in
    ``qc_fail`` become missing. Assays of the same protein on several panels are averaged."""
    notes = {}
    if annotation is None:
        L = raw.copy()
        if keep_samples is not None:
            L = L[L[sample_col].astype(str).isin(set(map(str, keep_samples)))]
        if qc_col is not None and qc_col in L.columns:
            bad = L[qc_col].astype(str).isin(set(qc_fail))
            notes["qc_masked_values"] = str(int(bad.sum()))
            L.loc[bad, value_col] = np.nan
        wide = L.pivot_table(index=sample_col, columns=assay_col, values=value_col, aggfunc="mean")
        ann = L[[assay_col, uniprot_col]].drop_duplicates(assay_col).set_index(assay_col)[uniprot_col]
    else:
        wide = raw.copy()
        ann = annotation.drop_duplicates(assay_col).set_index(assay_col)[uniprot_col]
    key = ann.map(lambda s: uniprot_key(split_ids(s)))
    notes["assays_without_uniprot"] = str(int(key.reindex(wide.columns).isna().sum()))
    X, feats = _aggregate(wide.apply(pd.to_numeric, errors="coerce"), key)
    return Adapted(X, feats, "uniprot", "olink", notes)


@layer_adapter("somascan", "uniprot")
def somascan(
    raw: pd.DataFrame,
    annotation: pd.DataFrame,
    aptamer_col: str = "SeqId",
    uniprot_col: str = "UniProt",
    keep: Optional[pd.Series] = None,
) -> Adapted:
    """SomaScan RFU (samples x aptamers) into UniProt space. ``annotation`` maps ``aptamer_col``
    to ``uniprot_col``; ``keep`` (boolean, indexed by aptamer) selects the human protein
    aptamers when the annotation has a type or organism column. RFU are log2 transformed when on
    the linear scale; several aptamers of one protein are averaged on the log2 scale."""
    ann = annotation.drop_duplicates(aptamer_col).set_index(aptamer_col)
    if keep is not None:
        ann = ann[keep.reindex(ann.index).fillna(False).astype(bool)]
    key = ann[uniprot_col].map(lambda s: uniprot_key(split_ids(s)))
    wide = raw.apply(pd.to_numeric, errors="coerce")
    wide, transform = _log2_if_linear(wide)
    X, feats = _aggregate(wide, key)
    return Adapted(X, feats, "uniprot", "somascan", {"transform": transform, "aptamers_without_uniprot": str(int(key.reindex(raw.columns).isna().sum()))})


@layer_adapter("ms_proteomics", "uniprot")
def ms_proteomics(
    raw: pd.DataFrame,
    id_col: str,
    sample_cols: Iterable[str],
    contaminant_patterns: Iterable[str] = ("CON__", "REV__", "Cont_", "contam", "DECOY", "rev_"),
    log2: Optional[bool] = None,
    preferred_accessions: Optional[Iterable[str]] = None,
) -> Adapted:
    """Protein group table (one row per group) into UniProt space: samples x leading accession.

    The leading protein is the first accession of ``id_col`` (MaxQuant ``Protein IDs`` or
    ``Majority protein IDs``, DIA-NN ``Protein.Group``, FASTA headers ``sp|P12345|NAME_HUMAN``
    are parsed). With ``preferred_accessions`` (for example the reviewed Swiss-Prot accessions of
    the HGNC table) the first accession of the group in that set is used instead, when there is one;
    DIA-NN sorts the accessions of a group alphabetically, so its first accession is often an
    unreviewed TrEMBL entry. Groups containing a contaminant or decoy pattern are dropped. Zero
    intensities are missing. With ``log2`` None the values are log2 transformed when on the linear
    scale."""
    t = raw.copy()
    ids = t[id_col].astype(str)
    bad = ids.apply(lambda s: any(p.lower() in s.lower() for p in contaminant_patterns))
    t = t[~bad]

    preferred = set(preferred_accessions) if preferred_accessions is not None else None

    def accession(tok: str) -> Optional[str]:
        if "|" in tok:  # sp|P12345|NAME
            parts = tok.split("|")
            tok = parts[1] if len(parts) > 1 else parts[0]
        return canonical_uniprot(tok)

    def leading(s: str) -> Optional[str]:
        accs = [accession(t) for t in re.split(r"[;\s]+", s.strip()) if t]
        if preferred is not None:
            for a in accs:
                if a in preferred:
                    return a
        return accs[0] if accs else None

    gid = pd.Index(t[id_col].astype(str))
    if gid.duplicated().any():
        gid = pd.Index([f"{g}#{i}" for i, g in enumerate(gid)])
    key = pd.Series(t[id_col].astype(str).map(leading).values, index=gid)
    wide = t[list(sample_cols)].apply(pd.to_numeric, errors="coerce").replace(0, np.nan).T
    wide.columns = gid
    transform = "as_provided"
    if log2 or (log2 is None and np.nanmax(wide.values) > LINEAR_SCALE_MAX):
        wide, transform = np.log2(wide.where(wide > 0)), "log2"
    X, feats = _aggregate(wide, key)
    notes = {"transform": transform, "groups_contaminant_or_decoy": str(int(bad.sum())), "groups_without_uniprot": str(int(key.isna().sum()))}
    if preferred is not None:
        notes["groups_keyed_outside_preferred"] = str(int((~key.dropna().isin(preferred)).sum()))
    return Adapted(X, feats, "uniprot", "ms_proteomics", notes)
