"""Cohort-as-environment feature selection.

The study indicator is treated as an environment variable (Peters et al. 2016; Pfister et
al. 2021). For each feature we estimate the standardised association with the outcome within
every training study, then combine the estimates with a random effects meta-analysis
(DerSimonian and Laird). A feature is *invariant* when its pooled effect is non-zero and the
between study heterogeneity is small (low I squared, Cochran Q not significant) and its sign
is consistent across studies. This is the single feature screen of stabilised regression
rather than full invariant causal prediction, which is intractable for thousands of genes.

Selectors (string names used in configs):
  variance      top k most variable features (the TACO "most variable genes" baseline)
  differential  top k by pooled effect size ignoring heterogeneity (meta-analytic DE)
  invariant     top k by pooled effect among features passing the heterogeneity screen
  prior         features from external causal evidence (MRCONTROL, TACO DCE, SPLIT lists)
  prior_invariant  intersection or union of prior genes with the invariant screen
  all           no selection

Every selector is fitted on training studies only.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
from scipy import stats


def _hedges_g(a: np.ndarray, b: np.ndarray) -> tuple:
    """Hedges g (bias corrected standardised mean difference) and its variance, columnwise."""
    n1, n2 = a.shape[0], b.shape[0]
    m1, m2 = np.nanmean(a, axis=0), np.nanmean(b, axis=0)
    v1, v2 = np.nanvar(a, axis=0, ddof=1), np.nanvar(b, axis=0, ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / max(n1 + n2 - 2, 1))
    sp = np.where(sp <= 0, np.nan, sp)
    d = (m1 - m2) / sp
    J = 1 - 3 / (4 * (n1 + n2) - 9)
    g = J * d
    var = (n1 + n2) / (n1 * n2) + g**2 / (2 * (n1 + n2))
    return g, var


def per_study_effects(X: pd.DataFrame, y: pd.Series, studies: pd.Series, min_per_class: int = 3) -> Dict[str, pd.DataFrame]:
    """Hedges g (cases minus controls) and variance per feature for each study."""
    out: Dict[str, pd.DataFrame] = {}
    for s in pd.unique(studies):
        idx = studies.index[studies == s]
        ys = y.loc[idx]
        pos = idx[(ys == 1).values]
        neg = idx[(ys == 0).values]
        if len(pos) < min_per_class or len(neg) < min_per_class:
            continue
        g, v = _hedges_g(X.loc[pos].values.astype(float), X.loc[neg].values.astype(float))
        out[s] = pd.DataFrame({"g": g, "var": v}, index=X.columns)
    return out


def meta_analyse(effects: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Random effects (DerSimonian and Laird) meta-analysis per feature.

    Returns pooled effect, standard error, z, p, Cochran Q, Q p-value, I squared, tau squared,
    number of studies, and the fraction of studies whose effect has the pooled sign.
    """
    if not effects:
        raise ValueError("no study had enough cases and controls for effect estimation")
    G = pd.concat({s: e["g"] for s, e in effects.items()}, axis=1)
    V = pd.concat({s: e["var"] for s, e in effects.items()}, axis=1)
    W = 1.0 / V
    k = G.notna().sum(axis=1)
    fixed = (W * G).sum(axis=1) / W.sum(axis=1)
    Q = (W * (G.sub(fixed, axis=0)) ** 2).sum(axis=1)
    df = (k - 1).clip(lower=1)
    C = W.sum(axis=1) - (W**2).sum(axis=1) / W.sum(axis=1)
    tau2 = ((Q - df) / C).clip(lower=0)
    Wr = 1.0 / (V.add(tau2, axis=0))
    pooled = (Wr * G).sum(axis=1) / Wr.sum(axis=1)
    se = np.sqrt(1.0 / Wr.sum(axis=1))
    z = pooled / se
    p = 2 * stats.norm.sf(np.abs(z))
    qp = stats.chi2.sf(Q, df)
    i2 = ((Q - df) / Q).clip(lower=0).fillna(0) * 100
    sign_consistency = (np.sign(G).eq(np.sign(pooled), axis=0) & G.notna()).sum(axis=1) / k
    return pd.DataFrame(
        {"effect": pooled, "se": se, "z": z, "p": p, "Q": Q, "Q_p": qp, "I2": i2, "tau2": tau2, "k": k, "sign_consistency": sign_consistency}
    )


class BaseSelector:
    name = "base"

    def fit(self, X: pd.DataFrame, y: pd.Series, studies: pd.Series) -> "BaseSelector":
        raise NotImplementedError

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        return X.reindex(columns=self.selected_)

    def fit_transform(self, X, y, studies):
        return self.fit(X, y, studies).transform(X)


class AllFeatures(BaseSelector):
    name = "all"

    def fit(self, X, y, studies):
        self.selected_ = list(X.columns)
        return self


class VarianceSelector(BaseSelector):
    name = "variance"

    def __init__(self, k: int = 500):
        self.k = k

    def fit(self, X, y, studies):
        # variance computed after per study centring so that study offsets do not drive it
        Xc = X.copy()
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            Xc.loc[idx] = X.loc[idx] - X.loc[idx].mean(axis=0)
        v = Xc.var(axis=0)
        self.selected_ = list(v.sort_values(ascending=False).index[: self.k])
        return self


class DifferentialSelector(BaseSelector):
    """Top k by absolute pooled effect (meta-analytic differential expression)."""

    name = "differential"

    def __init__(self, k: int = 500):
        self.k = k

    def fit(self, X, y, studies):
        m = y.notna()
        self.table_ = meta_analyse(per_study_effects(X.loc[m], y.loc[m], studies.loc[m]))
        self.selected_ = list(self.table_["effect"].abs().sort_values(ascending=False).index[: self.k])
        return self


class InvariantSelector(BaseSelector):
    """Top k by pooled effect among features whose association is stable across studies.

    Screen: I squared <= ``max_i2`` AND Q_p >= ``min_q_p`` AND sign consistency >= ``min_sign``
    AND at least ``min_k`` studies. When fewer than ``k`` features pass, the screen is relaxed
    by ranking on a stability score (|z| penalised by I squared) so that a model can still be
    fitted; the number passing the strict screen is stored in ``n_strict_``.
    """

    name = "invariant"

    def __init__(self, k: int = 200, max_i2: float = 50.0, min_q_p: float = 0.05, min_sign: float = 1.0, min_k: int = 2, max_p: float = 0.05):
        self.k, self.max_i2, self.min_q_p, self.min_sign, self.min_k, self.max_p = k, max_i2, min_q_p, min_sign, min_k, max_p

    def fit(self, X, y, studies):
        m = y.notna()
        t = meta_analyse(per_study_effects(X.loc[m], y.loc[m], studies.loc[m]))
        strict = t[(t["p"] <= self.max_p) & (t["I2"] <= self.max_i2) & (t["Q_p"] >= self.min_q_p) & (t["sign_consistency"] >= self.min_sign) & (t["k"] >= self.min_k)]
        self.n_strict_ = int(strict.shape[0])
        ranked_strict = strict["z"].abs().sort_values(ascending=False)
        if len(ranked_strict) >= self.k:
            sel = list(ranked_strict.index[: self.k])
        else:
            score = (t["z"].abs() * (1 - t["I2"] / 100.0) * t["sign_consistency"]).sort_values(ascending=False)
            sel = list(ranked_strict.index) + [f for f in score.index if f not in set(ranked_strict.index)][: self.k - len(ranked_strict)]
        self.table_ = t
        self.selected_ = sel
        return self


class PriorSelector(BaseSelector):
    """Features from an external causal evidence list, optionally intersected or unioned with
    the invariant screen. Missing prior genes (not in the matrix) are dropped silently."""

    name = "prior"

    def __init__(self, genes: Iterable[str], mode: str = "prior_only", k: int = 200, **invariant_kwargs):
        self.genes = list(genes)
        self.mode = mode  # prior_only | intersect | union
        self.k = k
        self.invariant_kwargs = invariant_kwargs

    def fit(self, X, y, studies):
        prior = [g for g in self.genes if g in X.columns]
        if self.mode == "prior_only":
            self.selected_ = prior
        else:
            inv = InvariantSelector(k=self.k, **self.invariant_kwargs).fit(X, y, studies)
            inv_set = set(inv.selected_)
            if self.mode == "intersect":
                self.selected_ = [g for g in prior if g in inv_set]
            elif self.mode == "union":
                self.selected_ = prior + [g for g in inv.selected_ if g not in set(prior)]
            else:
                raise ValueError(self.mode)
            self.invariant_table_ = inv.table_
        if len(self.selected_) == 0:
            raise ValueError("prior selector produced no features; check identifier space")
        return self


SELECTORS = {
    "all": AllFeatures,
    "variance": VarianceSelector,
    "differential": DifferentialSelector,
    "invariant": InvariantSelector,
    "prior": PriorSelector,
}


def make_selector(name: str, label: str | None = None, **kwargs) -> BaseSelector:
    """``label`` only names the configuration in results (see ``chorus.eval.loso.selector_label``)."""
    if name not in SELECTORS:
        raise KeyError(f"unknown selector {name}; choose from {sorted(SELECTORS)}")
    return SELECTORS[name](**kwargs)
