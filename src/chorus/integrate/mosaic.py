"""Mosaic multi-omics integration across cohorts with partially overlapping layers.

Three routes, in increasing ambition:

1. ``LateFusion``: one harmonised, LOSO validated model per layer; a subject's prediction is
   the missingness aware weighted mean of the layer models that exist for that subject
   (stacking with missingness indicators, as in MSB style late fusion). Works with any
   pattern of present and absent layers and is the primary CHORUS route for small cohorts.

2. ``knowledge_link``: diagonal integration of a transcriptomic and a proteomic cohort at the
   gene product level through an Ensembl to UniProt mapping, so that a model trained on one
   layer can be scored on the other at rank level. Expect modest transfer (median RNA to
   protein correlation about 0.6 for shared products).

3. ``mofa_factors``: MOFA+ with cohorts as groups and layers as views (optional dependency
   mofapy2). Factors are learned on training cohorts; held-out samples are scored by least
   squares projection onto the loadings of the views they have. Reported as a comparator.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


class LateFusion:
    """Combine per layer out-of-fold probabilities for the same subjects.

    ``fit`` takes a table with columns subject_id, layer, p (probability), y (label) built
    from LOSO predictions of the layer models; the meta learner is logistic regression on the
    layer probabilities plus indicators of which layers are present. Subjects with one layer
    are scored with that layer alone.
    """

    def __init__(self, layers: List[str]):
        self.layers = layers

    def _design(self, wide: pd.DataFrame) -> pd.DataFrame:
        cols = {}
        for l in self.layers:
            p = wide.get(l, pd.Series(np.nan, index=wide.index))
            cols[f"p_{l}"] = p.fillna(0.5)
            cols[f"has_{l}"] = p.notna().astype(float)
        return pd.DataFrame(cols, index=wide.index)

    def fit(self, table: pd.DataFrame) -> "LateFusion":
        wide = table.pivot_table(index="subject_id", columns="layer", values="p", aggfunc="mean")
        y = table.groupby("subject_id")["y"].first().loc[wide.index]
        self.meta_ = LogisticRegression(max_iter=2000).fit(self._design(wide), y)
        return self

    def predict_proba(self, table: pd.DataFrame) -> pd.Series:
        wide = table.pivot_table(index="subject_id", columns="layer", values="p", aggfunc="mean")
        return pd.Series(self.meta_.predict_proba(self._design(wide))[:, 1], index=wide.index)


def knowledge_link(X_rna: pd.DataFrame, X_prot: pd.DataFrame, mapping: pd.DataFrame, gene_col: str = "ensembl_gene", protein_col: str = "uniprot") -> Dict[str, pd.DataFrame]:
    """Reindex both matrices to shared gene products and rank transform within sample.

    Returns dict with keys ``rna`` and ``protein`` on a common ``product`` column space.
    """
    m = mapping[[gene_col, protein_col]].dropna().drop_duplicates()
    m = m[m[gene_col].isin(X_rna.columns) & m[protein_col].isin(X_prot.columns)]
    if m.empty:
        raise ValueError("no shared gene products between the two matrices")
    rna = X_rna[m[gene_col].values]
    rna.columns = m[protein_col].values
    rna = rna.T.groupby(level=0).mean().T
    prot = X_prot[sorted(set(m[protein_col]))]
    rna = rna[prot.columns]
    rank = lambda D: (D.rank(axis=1) - 1) / (D.shape[1] - 1)  # noqa: E731
    return {"rna": rank(rna), "protein": rank(prot)}


def mofa_factors(views: Dict[str, pd.DataFrame], groups: pd.Series, n_factors: int = 10, seed: int = 0, train_mask: Optional[pd.Series] = None) -> pd.DataFrame:
    """MOFA+ factors with cohorts as groups; held-out samples projected onto training loadings.

    ``views`` maps layer name to a samples by features matrix (samples missing from a view
    are allowed). ``groups`` maps sample to study. ``train_mask`` marks training samples; when
    given, the model is trained on them and all samples are scored by least squares
    projection using the views each sample has.
    """
    try:
        from mofapy2.run.entry_point import entry_point
    except ImportError as e:  # pragma: no cover
        raise ImportError("mofa_factors needs the optional dependency mofapy2") from e
    all_samples = pd.Index(sorted(set().union(*[set(v.index) for v in views.values()])))
    train = all_samples if train_mask is None else all_samples[train_mask.reindex(all_samples).fillna(False).astype(bool).values]
    ent = entry_point()
    data = []
    for name, V in views.items():
        Vt = V.reindex(train)
        data.append([Vt.loc[groups.loc[train] == g].values for g in pd.unique(groups.loc[train])])
    ent.set_data_options(scale_views=True)
    ent.set_data_matrix(data, views_names=list(views), groups_names=list(pd.unique(groups.loc[train])), samples_names=[list(train[groups.loc[train] == g]) for g in pd.unique(groups.loc[train])], features_names=[list(V.columns) for V in views.values()])
    ent.set_model_options(factors=n_factors)
    ent.set_train_options(seed=seed, verbose=False, convergence_mode="fast")
    ent.build()
    ent.run()
    W = {name: pd.DataFrame(ent.model.nodes["W"].getExpectation()[i], index=V.columns) for i, (name, V) in enumerate(views.items())}
    # project every sample with least squares on the concatenated loadings of its present views
    Z = pd.DataFrame(np.nan, index=all_samples, columns=[f"F{k+1}" for k in range(n_factors)])
    for s in all_samples:
        rows, ws = [], []
        for name, V in views.items():
            if s in V.index and V.loc[s].notna().any():
                x = V.loc[s]
                mu, sd = V.reindex(train).mean(), V.reindex(train).std().replace(0, 1)
                rows.append(((x - mu) / sd).fillna(0).values)
                ws.append(W[name].values)
        if rows:
            z, *_ = np.linalg.lstsq(np.vstack(ws), np.concatenate(rows), rcond=None)
            Z.loc[s] = z
    return Z
