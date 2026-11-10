"""Leave-one-study-out (LOSO) benchmark with harmonisation fitted inside the training fold.

For every held-out study S:
  1. the harmoniser is fitted on the training studies only (or, for a transductive
     comparator, on all studies without labels and flagged ``leak_risk``), then applied to S
     in addon mode;
  2. the feature selector is fitted on the training studies, with those studies as the
     environments of the invariance screen;
  3. optional clinical covariates (age, sex, ...) are appended;
  4. the classifier is fitted on the training studies and evaluated on S.

A nested option chooses the harmoniser by an inner LOSO over the training studies, so that
the harmonisation strategy itself is selected without touching S (Van et al. 2024 found no
universal winner, which is why the choice is validated rather than fixed).

``run_within_study_cv`` gives the internal cross validation estimate per study with the same
selector and classifier, so the report can show the optimism gap (Bernau et al. 2014).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from ..features import make_selector
from ..harmonise import make_harmoniser
from ..io import CohortSet
from ..models import make_classifier
from ..schema import binary_label
from .metrics import binary_metrics, bootstrap_ci


@dataclass
class LosoConfig:
    harmonisers: List[str] = field(default_factory=lambda: ["identity", "center_scale", "fsqn", "rank", "combat_addon"])
    selectors: List[Dict] = field(default_factory=lambda: [{"name": "variance", "k": 500}, {"name": "differential", "k": 200}, {"name": "invariant", "k": 200}])
    classifiers: List[str] = field(default_factory=lambda: ["logreg_en", "rf"])
    positive: str = "LC"
    negative: str = "COVID_recovered"
    covariates: List[str] = field(default_factory=list)  # metadata columns appended as features
    protect: List[str] = field(default_factory=list)  # metadata columns protected in combat_addon
    nested_select_harmoniser: bool = False
    inner_classifier: str = "logreg"
    min_test_per_class: int = 5
    min_train_studies: int = 2
    n_boot: int = 500
    random_state: int = 0
    sample_mask: Optional[pd.Series] = None  # optional boolean mask over pooled index (sensitivity analyses)


def _covariate_matrix(meta: pd.DataFrame, covariates: List[str]) -> pd.DataFrame:
    cols = {}
    for c in covariates:
        if c not in meta.columns:
            continue
        v = meta[c]
        if c == "sex":
            cols["cov_sex_female"] = (v.astype(str) == "female").astype(float).where(v.notna())
        elif v.dtype == object:
            for lvl in sorted(v.dropna().astype(str).unique())[:10]:
                cols[f"cov_{c}_{lvl}"] = (v.astype(str) == lvl).astype(float)
        else:
            cols[f"cov_{c}"] = pd.to_numeric(v, errors="coerce")
    return pd.DataFrame(cols, index=meta.index)


def _fit_harmoniser(name: str, X_tr, s_tr, y_tr, X_te, s_te, protect_df: Optional[pd.DataFrame]):
    kwargs = {}
    if name == "combat_addon" and protect_df is not None and protect_df.shape[1] > 0:
        kwargs["protect"] = protect_df
    H = make_harmoniser(name, **kwargs)
    if getattr(H, "transductive", False):
        X_all = pd.concat([X_tr, X_te])
        s_all = pd.concat([s_tr, s_te])
        H.fit(X_all, s_all)
        return H.transform(X_tr, s_tr), H.transform(X_te, s_te), True
    y_fit = y_tr if getattr(H, "requires_y", False) or name == "combat_addon" else None
    H.fit(X_tr, s_tr, y_fit)
    return H.transform(X_tr, s_tr), H.transform(X_te, s_te), False


def _inner_select_harmoniser(cands: List[str], X_tr, s_tr, y_tr, selector_cfg: Dict, cfg: LosoConfig, protect_df) -> Tuple[str, Dict[str, float]]:
    """Inner LOSO over the training studies to pick a harmoniser."""
    scores: Dict[str, List[float]] = {h: [] for h in cands}
    for inner_out in pd.unique(s_tr):
        m_te = s_tr == inner_out
        if (y_tr[m_te] == 1).sum() < cfg.min_test_per_class or (y_tr[m_te] == 0).sum() < cfg.min_test_per_class:
            continue
        if s_tr[~m_te].nunique() < 1:
            continue
        for h in cands:
            try:
                Xa, Xb, _ = _fit_harmoniser(h, X_tr[~m_te], s_tr[~m_te], y_tr[~m_te], X_tr[m_te], s_tr[m_te], None if protect_df is None else protect_df.loc[s_tr.index[~m_te]])
                sel = make_selector(**selector_cfg).fit(Xa, y_tr[~m_te], s_tr[~m_te])
                clf = make_classifier(cfg.inner_classifier, cfg.random_state)
                clf.fit(sel.transform(Xa), y_tr[~m_te])
                scores[h].append(roc_auc_score(y_tr[m_te], clf.predict_proba(sel.transform(Xb))[:, 1]))
            except Exception:
                scores[h].append(np.nan)
    means = {h: float(np.nanmean(v)) if len(v) else np.nan for h, v in scores.items()}
    best = max(means, key=lambda h: (np.nan_to_num(means[h], nan=-1)))
    return best, means


def run_loso(cs: CohortSet, cfg: LosoConfig, features: Optional[pd.Index] = None, verbose: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """Run the benchmark. Returns (results, predictions, selected_features)."""
    X_all, M_all = cs.pooled(features)
    y_all = binary_label(M_all, cfg.positive, cfg.negative)
    s_all = M_all["study_id"].astype(str)
    if cfg.sample_mask is not None:
        keep = cfg.sample_mask.reindex(X_all.index).fillna(False).astype(bool)
        X_all, M_all, y_all, s_all = X_all[keep], M_all[keep], y_all[keep], s_all[keep]
    labelled = y_all.notna()
    X_all, M_all, y_all, s_all = X_all[labelled], M_all[labelled], y_all[labelled], s_all[labelled]
    cov_all = _covariate_matrix(M_all, cfg.covariates)
    prot_all = _covariate_matrix(M_all, cfg.protect) if cfg.protect else None

    rows, preds, selected = [], [], {}
    harmonisers = list(cfg.harmonisers) + (["nested"] if cfg.nested_select_harmoniser else [])
    for held_out in cs.study_ids:
        if cs[held_out].source.get("train_only") == "true":  # contributes to training, never held out
            continue
        te = s_all == held_out
        tr = ~te
        if (y_all[te] == 1).sum() < cfg.min_test_per_class or (y_all[te] == 0).sum() < cfg.min_test_per_class:
            if verbose:
                print(f"[skip] {held_out}: too few labelled samples per class")
            continue
        if s_all[tr].nunique() < cfg.min_train_studies:
            continue
        X_tr, X_te, y_tr, s_tr, s_te = X_all[tr], X_all[te], y_all[tr], s_all[tr], s_all[te]
        protect_tr = None if prot_all is None else prot_all.loc[X_tr.index]
        for sel_cfg in cfg.selectors:
            for h in harmonisers:
                t0 = time.time()
                chosen, inner = h, {}
                if h == "nested":
                    chosen, inner = _inner_select_harmoniser(cfg.harmonisers, X_tr, s_tr, y_tr, sel_cfg, cfg, protect_tr)
                try:
                    Xh_tr, Xh_te, leak = _fit_harmoniser(chosen, X_tr, s_tr, y_tr, X_te, s_te, protect_tr)
                    sel = make_selector(**sel_cfg).fit(Xh_tr, y_tr, s_tr)
                except Exception as e:  # keep the benchmark running, record the failure
                    rows.append({"held_out": held_out, "harmoniser": h, "harmoniser_chosen": chosen, "selector": sel_cfg["name"], "classifier": None, "error": str(e)[:200]})
                    continue
                Xs_tr = pd.concat([sel.transform(Xh_tr), cov_all.loc[X_tr.index]], axis=1)
                Xs_te = pd.concat([sel.transform(Xh_te), cov_all.loc[X_te.index]], axis=1)
                selected[(held_out, h, sel_cfg["name"])] = list(sel.selected_)
                for clf_name in cfg.classifiers:
                    try:
                        clf = make_classifier(clf_name, cfg.random_state)
                        if clf_name == "cross_study":
                            clf.fit(Xs_tr, y_tr, studies=s_tr)
                        else:
                            clf.fit(Xs_tr, y_tr)
                        p = clf.predict_proba(Xs_te)[:, 1]
                    except Exception as e:
                        rows.append({"held_out": held_out, "harmoniser": h, "harmoniser_chosen": chosen, "selector": sel_cfg["name"], "classifier": clf_name, "error": str(e)[:200]})
                        continue
                    m = binary_metrics(y_te := y_all[te].values, p)
                    est, lo, hi = bootstrap_ci(y_te, p, n_boot=cfg.n_boot, seed=cfg.random_state)
                    rows.append(
                        {
                            "held_out": held_out,
                            "harmoniser": h,
                            "harmoniser_chosen": chosen,
                            "leak_risk": leak,
                            "selector": sel_cfg["name"],
                            "n_features": Xs_tr.shape[1],
                            "n_strict_invariant": getattr(sel, "n_strict_", np.nan),
                            "classifier": clf_name,
                            "n_train": int(tr.sum()),
                            "n_train_studies": int(s_tr.nunique()),
                            **m,
                            "auroc_lo": lo,
                            "auroc_hi": hi,
                            "inner_scores": inner or None,
                            "seconds": round(time.time() - t0, 1),
                        }
                    )
                    preds.append(pd.DataFrame({"sample": X_te.index, "held_out": held_out, "harmoniser": h, "selector": sel_cfg["name"], "classifier": clf_name, "y": y_te, "p": p}))
                    if verbose:
                        print(f"[{held_out}] {h:>13} | {sel_cfg['name']:>12} | {clf_name:>11} | AUROC {m['auroc']:.3f} ({lo:.2f} to {hi:.2f})")
    results = pd.DataFrame(rows)
    predictions = pd.concat(preds, ignore_index=True) if preds else pd.DataFrame()
    return results, predictions, selected


def run_within_study_cv(cs: CohortSet, cfg: LosoConfig, features: Optional[pd.Index] = None, n_splits: int = 5) -> pd.DataFrame:
    """Internal stratified CV per study with the same selectors and classifiers (optimism reference)."""
    X_all, M_all = cs.pooled(features)
    y_all = binary_label(M_all, cfg.positive, cfg.negative)
    s_all = M_all["study_id"].astype(str)
    cov_all = _covariate_matrix(M_all, cfg.covariates)
    rows = []
    for sid in cs.study_ids:
        m = (s_all == sid) & y_all.notna()
        X, y, s = X_all[m], y_all[m], s_all[m]
        if (y == 1).sum() < n_splits or (y == 0).sum() < n_splits:
            continue
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=cfg.random_state)
        for sel_cfg in cfg.selectors:
            for clf_name in cfg.classifiers:
                if clf_name == "cross_study":
                    continue
                p = np.zeros(len(y))
                ok = True
                for tr_i, te_i in skf.split(X, y):
                    try:
                        sel_name = sel_cfg["name"]
                        # within one study there is a single environment; invariant falls back to differential
                        sc = dict(sel_cfg)
                        if sel_name == "invariant":
                            sc = {"name": "differential", "k": sel_cfg.get("k", 200)}
                        sel = make_selector(**sc).fit(X.iloc[tr_i], y.iloc[tr_i], s.iloc[tr_i])
                        Xtr = pd.concat([sel.transform(X.iloc[tr_i]), cov_all.loc[X.index[tr_i]]], axis=1)
                        Xte = pd.concat([sel.transform(X.iloc[te_i]), cov_all.loc[X.index[te_i]]], axis=1)
                        clf = make_classifier(clf_name, cfg.random_state).fit(Xtr, y.iloc[tr_i])
                        p[te_i] = clf.predict_proba(Xte)[:, 1]
                    except Exception:
                        ok = False
                        break
                if not ok:
                    continue
                mt = binary_metrics(y.values, p)
                rows.append({"study": sid, "selector": sel_cfg["name"], "classifier": clf_name, **mt})
    return pd.DataFrame(rows)


def summarise(results: pd.DataFrame, by: Tuple[str, ...] = ("harmoniser", "selector", "classifier")) -> pd.DataFrame:
    """Mean, median and spread of AUROC across held-out studies for each configuration."""
    ok = results[results.get("error").isna()] if "error" in results.columns else results
    g = ok.groupby(list(by))
    out = g["auroc"].agg(["mean", "median", "std", "min", "max", "count"]).rename(columns={"count": "n_studies"})
    out["auprc_mean"] = g["auprc"].mean()
    out["brier_mean"] = g["brier"].mean()
    return out.sort_values("mean", ascending=False)
