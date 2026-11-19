"""Harmonisation operators with leave-one-study-out safe fit/transform semantics.

Design rule. ``fit`` receives only training studies. ``transform`` may receive a study that
was never seen during ``fit``. In that case the operator must adjust the new study with
information from that study alone and without its labels (the addon setting of Hornung et al.
2017, bapred). Operators therefore never store anything about a held-out study and never
look at ``y`` of the study being transformed.

Operators (string names used in configs):
  identity          no adjustment (pool as is)
  center_scale      per study, per feature z-score (location and scale)
  center            per study, per feature mean centring only
  quantile_ref      sample-wise quantile normalisation to a reference distribution learned on
                    the training studies (addon quantile normalisation)
  fsqn              feature-specific quantile normalisation of each study to the training
                    distribution of that feature (Franks et al. 2018)
  rank              within-sample rank transform (single sample, no fit needed)
  gene_pair         top scoring gene pair indicators selected on training studies
  combat_addon      location/scale ComBat style adjustment with addon estimation for a new
                    study (no covariates used for the new study)
  combat_pooled     ComBat (inmoose) fitted on all studies passed to fit; transductive, labels
                    never used; provided as a comparator and flagged as such in reports

All operators act on a continuous matrix (log CPM for RNA, NPX or log intensity for
proteins). Count matrices must be converted first (see ``transcriptomics.log_cpm``).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Optional

import numpy as np
import pandas as pd


class Harmoniser(ABC):
    name: str = "base"
    transductive: bool = False  # True when the operator must see the test study at fit time

    @abstractmethod
    def fit(self, X: pd.DataFrame, studies: pd.Series, y: Optional[pd.Series] = None) -> "Harmoniser":
        ...

    @abstractmethod
    def transform(self, X: pd.DataFrame, studies: pd.Series) -> pd.DataFrame:
        ...

    def fit_transform(self, X: pd.DataFrame, studies: pd.Series, y: Optional[pd.Series] = None) -> pd.DataFrame:
        return self.fit(X, studies, y).transform(X, studies)


class Identity(Harmoniser):
    name = "identity"

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        return self

    def transform(self, X, studies):
        return X.reindex(columns=self.features_)


class StudyCenterScale(Harmoniser):
    """Per study z-score of every feature. Addon: a new study is standardised with its own
    label free mean and standard deviation. ``scale=False`` gives mean centring only."""

    name = "center_scale"

    def __init__(self, scale: bool = True, eps: float = 1e-8):
        self.scale = scale
        self.eps = eps

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        return self

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        out = X.copy()
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            block = X.loc[idx]
            mu = block.mean(axis=0)
            out.loc[idx] = block - mu
            if self.scale:
                sd = block.std(axis=0, ddof=1).replace(0, np.nan).fillna(1.0)
                out.loc[idx] = out.loc[idx] / (sd + self.eps)
        return out


class QuantileToReference(Harmoniser):
    """Sample-wise quantile normalisation to a reference quantile curve learned on training
    samples. Each new sample is mapped independently, so it is addon safe."""

    name = "quantile_ref"

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        sorted_vals = np.sort(X.values, axis=1)
        self.reference_ = np.nanmean(sorted_vals, axis=0)
        return self

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        vals = X.values.astype(float)
        out = np.empty_like(vals)
        ref = self.reference_
        for i in range(vals.shape[0]):
            row = vals[i]
            order = np.argsort(row, kind="mergesort")
            ranks = np.empty_like(order)
            ranks[order] = np.arange(len(row))
            out[i] = ref[ranks]
        return pd.DataFrame(out, index=X.index, columns=X.columns)


class FeatureQuantileMap(Harmoniser):
    """Feature-specific quantile normalisation (FSQN). For every feature, the empirical
    distribution within a study is mapped onto the distribution of that feature in the
    training studies. Uses only the study's own values, so it is addon safe."""

    name = "fsqn"

    def __init__(self, n_quantiles: int = 101):
        self.n_quantiles = n_quantiles

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        q = np.linspace(0, 1, self.n_quantiles)
        self.q_ = q
        self.ref_quantiles_ = np.nanquantile(X.values.astype(float), q, axis=0)  # (nq, p)
        return self

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        out = X.copy().astype(float)
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            block = X.loc[idx].values.astype(float)
            n = block.shape[0]
            if n < 3:
                continue
            # empirical cdf position of each value within its study, then map to reference quantile
            ranks = block.argsort(axis=0).argsort(axis=0) / max(n - 1, 1)
            mapped = np.empty_like(block)
            for j in range(block.shape[1]):
                mapped[:, j] = np.interp(ranks[:, j], self.q_, self.ref_quantiles_[:, j])
            out.loc[idx] = mapped
        return out


class WithinSampleRank(Harmoniser):
    """Within-sample rank transform scaled to [0, 1]. Removes per study location and scale
    effects by construction. No fit needed."""

    name = "rank"

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        return self

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        r = X.rank(axis=1, method="average")
        return (r - 1) / (X.shape[1] - 1)


class GenePairRank(Harmoniser):
    """Top scoring pair style indicators: for selected pairs (a, b), feature = 1[x_a > x_b].

    Pairs are chosen on training data by the absolute difference in P(x_a > x_b) between
    classes, pooled across training studies with each study weighted equally. Requires ``y``
    at fit time (training labels only). Output is a new feature matrix.
    """

    name = "gene_pair"
    requires_y = True

    def __init__(self, n_pairs: int = 200, n_candidate_features: int = 300, random_state: int = 0):
        self.n_pairs = n_pairs
        self.n_candidate_features = n_candidate_features
        self.random_state = random_state

    def fit(self, X, studies, y=None):
        if y is None:
            raise ValueError("gene_pair harmoniser needs training labels")
        yv = y.loc[X.index].values
        ok = ~np.isnan(yv)
        Xv = X.values[ok].astype(float)
        yv = yv[ok]
        sv = studies.loc[X.index].values[ok]
        # candidate features: highest variance of within-sample rank
        ranks = X.loc[X.index[ok]].rank(axis=1).values
        var = ranks.var(axis=0)
        cand = np.argsort(-var)[: self.n_candidate_features]
        scores: Dict[tuple, float] = {}
        for s in np.unique(sv):
            m = sv == s
            if (yv[m] == 1).sum() < 3 or (yv[m] == 0).sum() < 3:
                continue
            Xs = Xv[m][:, cand]
            pos, neg = Xs[yv[m] == 1], Xs[yv[m] == 0]
            # P(x_a > x_b) per class, vectorised over pairs
            p_pos = (pos[:, :, None] > pos[:, None, :]).mean(axis=0)
            p_neg = (neg[:, :, None] > neg[:, None, :]).mean(axis=0)
            d = p_pos - p_neg
            iu = np.triu_indices(len(cand), k=1)
            for a, b, val in zip(iu[0], iu[1], d[iu]):
                scores[(a, b)] = scores.get((a, b), 0.0) + val
        n_env = max(1, len(np.unique(sv)))
        ranked = sorted(scores.items(), key=lambda kv: -abs(kv[1]) / n_env)[: self.n_pairs]
        self.pairs_ = [(X.columns[cand[a]], X.columns[cand[b]]) for (a, b), _ in ranked]
        self.features_ = X.columns
        return self

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        cols = {f"{a}>{b}": (X[a] > X[b]).astype(float) for a, b in self.pairs_}
        return pd.DataFrame(cols, index=X.index)


class AddonComBat(Harmoniser):
    """Location/scale ComBat style adjustment with addon estimation.

    Training studies: for each feature, the study location and variance are estimated and each
    study is shifted and scaled to the grand mean and pooled variance of the training pool. The
    outcome and any protected covariates enter the per study design centred within the study,
    so the location (the intercept) is the study mean and the covariates only remove their
    variation from the variance estimate; a covariate constant within a study is dropped from
    that study's design. Optional shrinkage of the batch parameters toward their across-feature
    average with weight ``shrink`` (0 = none, the default, approximating limma
    removeBatchEffect; 1 = full pooling); with ``shrink`` > 0 part of each feature's study
    offset is kept, so the study means are no longer equal after transform.

    New study at transform time: location and scale are estimated from that study alone,
    without covariates, exactly as the bapred "addon" variant. This keeps the held-out study
    label free. When the held-out study has a very different case to control ratio, this
    operator can remove part of the disease signal; the evaluation layer reports this
    through the comparison with ``center_scale`` and ``identity``.
    """

    name = "combat_addon"

    def __init__(self, protect: Optional[pd.DataFrame] = None, shrink: float = 0.0, eps: float = 1e-8):
        self.protect = protect
        self.shrink = shrink
        self.eps = eps

    def _batch_params(self, block: np.ndarray, Z: Optional[np.ndarray]) -> tuple:
        if Z is not None:
            # a covariate constant within the study (for example sample_matrix in a single matrix
            # cohort) is collinear with the intercept, so it is dropped
            Z = Z[:, np.ptp(Z, axis=0) > 0] if Z.shape[0] else Z
            if Z.shape[1] == 0:
                Z = None
        if Z is not None and Z.shape[0] == block.shape[0] and Z.shape[0] > Z.shape[1] + 2:
            # covariates centred within the study: the intercept is then the study mean
            Zc = np.column_stack([np.ones(len(block)), Z - Z.mean(axis=0)])
            beta, *_ = np.linalg.lstsq(Zc, block, rcond=None)
            resid = block - Zc @ beta
            mu = beta[0]
            var = resid.var(axis=0, ddof=Zc.shape[1])
        else:
            mu = block.mean(axis=0)
            var = block.var(axis=0, ddof=1)
        return mu, var

    def fit(self, X, studies, y=None):
        self.features_ = X.columns
        Xv = X.values.astype(float)
        self.grand_mean_ = Xv.mean(axis=0)
        self.pooled_var_ = Xv.var(axis=0, ddof=1) + self.eps
        self.batch_mu_: Dict[str, np.ndarray] = {}
        self.batch_var_: Dict[str, np.ndarray] = {}
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            Z = None
            if self.protect is not None:
                Zdf = self.protect.reindex(idx)
                if y is not None:
                    Zdf = pd.concat([Zdf, y.reindex(idx).rename("_y")], axis=1)
                Z = Zdf.apply(pd.to_numeric, errors="coerce").fillna(Zdf.apply(pd.to_numeric, errors="coerce").mean()).values
            elif y is not None:
                Z = y.reindex(idx).fillna(y.mean()).values.reshape(-1, 1)
            mu, var = self._batch_params(X.loc[idx].values.astype(float), Z)
            self.batch_mu_[s] = self._shrink(mu - self.grand_mean_) + self.grand_mean_
            self.batch_var_[s] = np.exp(self._shrink(np.log(var + self.eps) - np.log(self.pooled_var_))) * self.pooled_var_
        return self

    def _shrink(self, dev: np.ndarray) -> np.ndarray:
        return (1 - self.shrink) * dev + self.shrink * dev.mean()

    def transform(self, X, studies):
        X = X.reindex(columns=self.features_)
        out = X.copy().astype(float)
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            block = X.loc[idx].values.astype(float)
            if s in self.batch_mu_:
                mu, var = self.batch_mu_[s], self.batch_var_[s]
            else:  # addon: label free estimate from the new study alone
                mu, var = self._batch_params(block, None)
                mu = self._shrink(mu - self.grand_mean_) + self.grand_mean_
                var = np.exp(self._shrink(np.log(var + self.eps) - np.log(self.pooled_var_))) * self.pooled_var_
            out.loc[idx] = (block - mu) / np.sqrt(var + self.eps) * np.sqrt(self.pooled_var_) + self.grand_mean_
        return out


class PooledComBat(Harmoniser):
    """Transductive ComBat via inmoose, fitted on every study passed to ``fit`` (which, in the
    evaluation layer, includes the unlabelled held-out study). Labels are never used. Kept as
    a comparator because it is what most published pooled analyses do."""

    name = "combat_pooled"
    transductive = True

    def fit(self, X, studies, y=None):
        try:
            from inmoose.pycombat import pycombat_norm
        except ImportError as e:  # pragma: no cover
            raise ImportError("combat_pooled needs the optional dependency inmoose (pip install inmoose)") from e
        self.features_ = X.columns
        adj = pycombat_norm(X.values.T.astype(float), studies.loc[X.index].astype(str).values)
        self.adjusted_ = pd.DataFrame(np.asarray(adj).T, index=X.index, columns=X.columns)
        return self

    def transform(self, X, studies):
        missing = X.index.difference(self.adjusted_.index)
        if len(missing):
            raise RuntimeError("combat_pooled is transductive: all samples must be passed to fit")
        return self.adjusted_.loc[X.index]


HARMONISERS = {
    "identity": Identity,
    "center_scale": StudyCenterScale,
    "center": lambda: StudyCenterScale(scale=False),
    "quantile_ref": QuantileToReference,
    "fsqn": FeatureQuantileMap,
    "rank": WithinSampleRank,
    "gene_pair": GenePairRank,
    "combat_addon": AddonComBat,
    "combat_pooled": PooledComBat,
}


def make_harmoniser(name: str, **kwargs) -> Harmoniser:
    if name not in HARMONISERS:
        raise KeyError(f"unknown harmoniser {name}; choose from {sorted(HARMONISERS)}")
    factory = HARMONISERS[name]
    return factory(**kwargs) if kwargs else factory()
