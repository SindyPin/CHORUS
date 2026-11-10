"""Invariance anchored harmonisation (the CHORUS operator).

Motivation. Every location/scale harmoniser estimates a study's technical shift from all
genes of that study. When the case to control ratio differs between studies, that estimate
absorbs part of the disease signal: removing it removes biology, protecting it with the
outcome covariate fabricates biology under imbalance (Nygaard et al. 2016), and in addon
mode on an unlabelled cohort there is no covariate to protect with at all.

Idea. With cohorts as environments, a per gene random effects meta-analysis of the disease
association partitions the genes into three sets:

  invariant   consistent association in every training cohort (candidate causal predictors)
  null        no association in any training cohort and low heterogeneity (negative controls
              that are null in every environment, not only in one dataset)
  unstable    association present in some cohorts and absent or reversed in others (shortcuts)

The technical factors of each sample are estimated from the null genes only (a RUV style
factor model, Gagnon-Bartsch and Speed 2012) and removed from all genes. For a held out
cohort the factors are estimated from that cohort's own null gene values, which is label
free and by construction cannot remove disease signal however unbalanced the cohort is.
The classifier is then trained on the invariant set, and the unstable set supports a per
cohort concordance diagnostic (``sign_agreement``) that flags cohorts whose signal runs
against the training environments.

Model. Let Y be log expression centred by the training grand mean, C the null gene set and
k the number of technical factors (default: number of training studies). On the training
samples, W (n x k) is the leading left singular subspace of Y_C scaled by the singular
values, and A (k x p) = argmin ||Y - W A||. For any sample block Y_s (training or new),
W_s = Y_{s,C} A_C^T (A_C A_C^T)^{-1} and the harmonised values are (Y_s - W_s A) / sd,
where sd is the per gene residual standard deviation on the training samples. A is estimated
jointly with the training outcome (RUV-4 style) so that disease effects are not absorbed into
the loadings of disease genes under case imbalance. By default k is the number of training
studies plus ``extra_factors``, and a residual label free per study centring and scaling is
applied after factor removal (``post_center_scale``).
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from ..features.invariant import meta_analyse, per_study_effects
from .operators import Harmoniser


class InvarianceAnchored(Harmoniser):
    name = "chorus"
    requires_y = True

    def __init__(
        self,
        n_controls: int = 800,
        max_abs_g: float = 0.35,
        null_p_min: float = 0.2,
        max_i2: float = 50.0,
        k_factors: Optional[int] = None,
        min_variance_quantile: float = 0.2,
        invariant_max_p: float = 0.05,
        unstable_min_g: float = 0.5,
        extra_factors: int = 4,
        post_center_scale: bool = True,
        eps: float = 1e-8,
    ):
        self.n_controls = n_controls
        self.max_abs_g = max_abs_g
        self.null_p_min = null_p_min
        self.max_i2 = max_i2
        self.k_factors = k_factors
        self.min_variance_quantile = min_variance_quantile
        self.invariant_max_p = invariant_max_p
        self.unstable_min_g = unstable_min_g
        self.extra_factors = extra_factors
        self.post_center_scale = post_center_scale
        self.eps = eps

    # ------------------------------------------------------------------ partition
    def partition(self, X: pd.DataFrame, y: pd.Series, studies: pd.Series) -> Dict[str, List[str]]:
        m = y.notna()
        effects = per_study_effects(X.loc[m], y.loc[m], studies.loc[m])
        if len(effects) < 2:
            raise ValueError("invariance anchored harmonisation needs at least two labelled training studies")
        t = meta_analyse(effects)
        G = pd.concat({s: e["g"] for s, e in effects.items()}, axis=1)
        max_abs = G.abs().max(axis=1)
        var = X.loc[m].var(axis=0)
        var_ok = var >= var.quantile(self.min_variance_quantile)
        null_mask = (t["p"] > self.null_p_min) & (max_abs < self.max_abs_g) & (t["I2"] <= self.max_i2) & var_ok & max_abs.notna()
        controls = list(max_abs[null_mask].sort_values().index[: self.n_controls])
        inv_mask = (t["p"] <= self.invariant_max_p) & (t["I2"] <= self.max_i2) & (t["Q_p"] >= 0.05) & (t["sign_consistency"] >= 1.0)
        invariant = list(t.loc[inv_mask, "z"].abs().sort_values(ascending=False).index)
        unstable_mask = (max_abs >= self.unstable_min_g) & ~inv_mask
        unstable = list(max_abs[unstable_mask].sort_values(ascending=False).index)
        self.table_ = t
        self.per_study_g_ = G
        return {"invariant": invariant, "null": controls, "unstable": unstable}

    # ------------------------------------------------------------------ fit / transform
    def fit(self, X: pd.DataFrame, studies: pd.Series, y: Optional[pd.Series] = None) -> "InvarianceAnchored":
        if y is None:
            raise ValueError("chorus harmoniser needs training labels to choose the null gene set")
        self.features_ = X.columns
        parts = self.partition(X, y, studies)
        self.controls_ = parts["null"]
        self.invariant_ = parts["invariant"]
        self.unstable_ = parts["unstable"]
        if len(self.controls_) < 20:
            raise ValueError(f"only {len(self.controls_)} null genes passed the screen; relax max_abs_g or null_p_min")
        Xv = X.values.astype(float)
        self.mu_ = np.nanmean(Xv, axis=0)
        Y = Xv - self.mu_
        cidx = X.columns.get_indexer(self.controls_)
        k = self.k_factors or int(studies.loc[X.index].nunique()) + self.extra_factors
        k = max(1, min(k, len(self.controls_) - 1, Y.shape[0] - 1))
        U, S, _ = np.linalg.svd(Y[:, cidx], full_matrices=False)
        W = U[:, :k] * S[:k]
        # loadings are estimated jointly with the training outcome (RUV-4 style), so that the
        # disease effect is not absorbed into the technical loadings of disease genes when case
        # fractions differ between studies; the outcome column is dropped at transform time
        yv = pd.to_numeric(y.reindex(X.index), errors="coerce").values.astype(float)
        yv = np.where(np.isnan(yv), np.nanmean(yv), yv) - np.nanmean(yv)
        D = np.column_stack([W, yv])
        coef, *_ = np.linalg.lstsq(D, Y, rcond=None)  # (k + 1) x p
        A = coef[:k]
        self.beta_ = coef[k]
        self.A_ = A
        self.A_c_ = A[:, cidx]
        self.k_ = k
        resid = Y - D @ coef
        self.sd_ = resid.std(axis=0, ddof=1) + self.eps
        self.sd_[~np.isfinite(self.sd_)] = 1.0
        return self

    def _factors(self, Yc: np.ndarray) -> np.ndarray:
        # least squares W = Y_C A_C^T (A_C A_C^T)^-1, computed per sample block
        G = self.A_c_ @ self.A_c_.T
        return Yc @ self.A_c_.T @ np.linalg.pinv(G)

    def transform(self, X: pd.DataFrame, studies: pd.Series) -> pd.DataFrame:
        X = X.reindex(columns=self.features_)
        Y = X.values.astype(float) - self.mu_
        cidx = X.columns.get_indexer(self.controls_)
        out = np.empty_like(Y)
        for s in pd.unique(studies):
            rows = np.where(studies.loc[X.index].values == s)[0]
            W = self._factors(Y[rows][:, cidx])
            Z = (Y[rows] - W @ self.A_) / self.sd_
            if self.post_center_scale and len(rows) > 2:
                # residual per study location and scale, label free; the factor step has already
                # removed the bulk of the study shift from anchored (null gene) information
                Z = Z - Z.mean(axis=0)
                sd = Z.std(axis=0, ddof=1)
                Z = Z / np.where(sd > 0, sd, 1.0)
            out[rows] = Z
        return pd.DataFrame(out, index=X.index, columns=X.columns)


def sign_agreement(pooled_effect: pd.Series, X_te: pd.DataFrame, y_te: pd.Series, genes: List[str], min_per_class: int = 3) -> Dict[str, float]:
    """Concordance of a held out cohort with the training environments.

    For the given genes (normally the training invariant set), the fraction whose Hedges g in
    the held out cohort has the same sign as the pooled training effect, and the Spearman
    correlation of the two effect vectors. Values near 0.5 and 0 mean the cohort's signal is
    unrelated to the training cohorts; values well below 0.5 mean it runs against them, as
    happens when group is confounded with a technical batch in that cohort.
    """
    genes = [g for g in genes if g in X_te.columns and g in pooled_effect.index]
    yv = y_te.loc[X_te.index]
    pos, neg = X_te.index[(yv == 1).values], X_te.index[(yv == 0).values]
    if len(genes) == 0 or len(pos) < min_per_class or len(neg) < min_per_class:
        return {"n_genes": len(genes), "sign_agreement": np.nan, "effect_spearman": np.nan}
    a, b = X_te.loc[pos, genes].values.astype(float), X_te.loc[neg, genes].values.astype(float)
    sp = np.sqrt(((len(pos) - 1) * a.var(axis=0, ddof=1) + (len(neg) - 1) * b.var(axis=0, ddof=1)) / max(len(pos) + len(neg) - 2, 1))
    g_te = (a.mean(axis=0) - b.mean(axis=0)) / np.where(sp <= 0, np.nan, sp)
    g_tr = pooled_effect.loc[genes].values
    ok = np.isfinite(g_te) & np.isfinite(g_tr)
    if ok.sum() == 0:
        return {"n_genes": len(genes), "sign_agreement": np.nan, "effect_spearman": np.nan}
    agree = float((np.sign(g_te[ok]) == np.sign(g_tr[ok])).mean())
    rho = float(pd.Series(g_te[ok]).corr(pd.Series(g_tr[ok]), method="spearman")) if ok.sum() > 2 else np.nan
    return {"n_genes": int(ok.sum()), "sign_agreement": agree, "effect_spearman": rho}
