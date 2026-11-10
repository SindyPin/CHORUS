"""Classifiers used in the CHORUS benchmark.

Names (configs):
  logreg_en       elastic net logistic regression (strong linear baseline)
  rf              random forest
  extratrees      extremely randomised trees
  svc             RBF support vector classifier with probability outputs
  mlp             multilayer perceptron
  xgb             gradient boosting (optional dependency xgboost)
  tabpfn          TabPFN v2 foundation model (optional dependency tabpfn)
  taco            soft voting ensemble of tabpfn + rf + extratrees + svc + mlp, as in TACO
                  (falls back to the four model ensemble when tabpfn is not installed or
                  the environment sets CHORUS_DISABLE_TABPFN=1)
  cross_study     ensemble of per study learners weighted by cross study AUROC
                  (Patil and Parmigiani 2018); the merging versus ensembling comparator

All estimators expose ``fit(X, y)`` and ``predict_proba(X)``. Missing values are imputed with
the training median inside a pipeline so that proteomics matrices with NaN are accepted.
"""

from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier, VotingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


def _pipe(est, scale: bool = True) -> Pipeline:
    steps = [("impute", SimpleImputer(strategy="median"))]
    if scale:
        steps.append(("scale", StandardScaler()))
    steps.append(("clf", est))
    return Pipeline(steps)


def _tabpfn(random_state: int):
    if os.environ.get("CHORUS_DISABLE_TABPFN") == "1":  # for example when the TabPFN licence token is unavailable
        return None
    try:
        from tabpfn import TabPFNClassifier
    except ImportError:
        return None
    try:
        return TabPFNClassifier(random_state=random_state, ignore_pretraining_limits=True)
    except TypeError:  # older API
        return TabPFNClassifier(seed=random_state)


def make_classifier(name: str, random_state: int = 0, **kwargs):
    rs = random_state
    if name == "logreg":
        return _pipe(LogisticRegression(C=kwargs.get("C", 0.5), max_iter=2000, class_weight="balanced"))
    if name == "logreg_en":
        import sklearn

        major, minor = (int(v) for v in sklearn.__version__.split(".")[:2])
        lr_kwargs = dict(solver="saga", l1_ratio=kwargs.get("l1_ratio", 0.5), C=kwargs.get("C", 0.5), max_iter=5000, class_weight="balanced", random_state=rs)
        if (major, minor) < (1, 8):
            lr_kwargs["penalty"] = "elasticnet"
        return _pipe(LogisticRegression(**lr_kwargs))
    if name == "rf":
        return _pipe(RandomForestClassifier(n_estimators=kwargs.get("n_estimators", 500), class_weight="balanced_subsample", random_state=rs, n_jobs=-1), scale=False)
    if name == "extratrees":
        return _pipe(ExtraTreesClassifier(n_estimators=kwargs.get("n_estimators", 500), class_weight="balanced_subsample", random_state=rs, n_jobs=-1), scale=False)
    if name == "svc":
        return _pipe(SVC(C=kwargs.get("C", 1.0), kernel="rbf", probability=True, class_weight="balanced", random_state=rs))
    if name == "mlp":
        return _pipe(MLPClassifier(hidden_layer_sizes=kwargs.get("hidden", (64,)), alpha=1e-2, max_iter=1000, random_state=rs))
    if name == "xgb":
        try:
            from xgboost import XGBClassifier
        except ImportError as e:
            raise ImportError("xgb needs the optional dependency xgboost") from e
        return _pipe(XGBClassifier(n_estimators=400, max_depth=3, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8, random_state=rs, eval_metric="logloss"), scale=False)
    if name == "tabpfn":
        est = _tabpfn(rs)
        if est is None:
            raise ImportError("tabpfn needs the optional dependency tabpfn")
        return _pipe(est, scale=False)
    if name == "taco":
        members = [
            ("rf", make_classifier("rf", rs)),
            ("extratrees", make_classifier("extratrees", rs)),
            ("svc", make_classifier("svc", rs)),
            ("mlp", make_classifier("mlp", rs)),
        ]
        tp = _tabpfn(rs)
        if tp is not None:
            members.insert(0, ("tabpfn", _pipe(tp, scale=False)))
        return VotingClassifier(members, voting="soft", n_jobs=1)
    if name == "cross_study":
        return CrossStudyEnsemble(base=kwargs.get("base", "logreg"), random_state=rs)
    raise KeyError(f"unknown classifier {name}")


class CrossStudyEnsemble(BaseEstimator, ClassifierMixin):
    """Fit one learner per training study; weight each by its mean AUROC on the other
    training studies (zero weight when below 0.5). Predictions are the weighted mean of
    member probabilities. Requires ``studies`` passed to ``fit``."""

    def __init__(self, base: str = "logreg", random_state: int = 0):
        self.base = base
        self.random_state = random_state

    def fit(self, X, y, studies: Optional[pd.Series] = None):
        if studies is None:
            raise ValueError("CrossStudyEnsemble.fit needs studies")
        X = pd.DataFrame(X)
        y = pd.Series(np.asarray(y), index=X.index)
        studies = pd.Series(np.asarray(studies), index=X.index)
        self.members_: Dict[str, object] = {}
        for s in pd.unique(studies):
            idx = studies.index[studies == s]
            ys = y.loc[idx]
            if ys.nunique() < 2 or (ys == 1).sum() < 3 or (ys == 0).sum() < 3:
                continue
            m = make_classifier(self.base, self.random_state)
            m.fit(X.loc[idx], ys)
            self.members_[s] = m
        weights = {}
        for s, m in self.members_.items():
            aucs = []
            for t in self.members_:
                if t == s:
                    continue
                idx = studies.index[studies == t]
                if y.loc[idx].nunique() < 2:
                    continue
                aucs.append(roc_auc_score(y.loc[idx], m.predict_proba(X.loc[idx])[:, 1]))
            w = np.mean(aucs) - 0.5 if aucs else 0.5
            weights[s] = max(w, 0.0)
        tot = sum(weights.values())
        self.weights_ = {s: (w / tot if tot > 0 else 1.0 / len(weights)) for s, w in weights.items()}
        self.classes_ = np.array([0, 1])
        return self

    def predict_proba(self, X):
        X = pd.DataFrame(X)
        p = np.zeros(X.shape[0])
        for s, m in self.members_.items():
            p += self.weights_[s] * m.predict_proba(X)[:, 1]
        return np.column_stack([1 - p, p])

    def predict(self, X):
        return (self.predict_proba(X)[:, 1] >= 0.5).astype(int)


CLASSIFIERS = ["logreg", "logreg_en", "rf", "extratrees", "svc", "mlp", "xgb", "tabpfn", "taco", "cross_study"]
