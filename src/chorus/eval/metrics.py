"""Metrics with bootstrap confidence intervals."""

from __future__ import annotations

from typing import Callable, Dict, Tuple

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    matthews_corrcoef,
    roc_auc_score,
)


def expected_calibration_error(y: np.ndarray, p: np.ndarray, n_bins: int = 10) -> float:
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    ece = 0.0
    for b in range(n_bins):
        m = idx == b
        if m.any():
            ece += m.mean() * abs(y[m].mean() - p[m].mean())
    return float(ece)


def binary_metrics(y: np.ndarray, p: np.ndarray, threshold: float = 0.5) -> Dict[str, float]:
    y = np.asarray(y).astype(int)
    p = np.asarray(p).astype(float)
    yhat = (p >= threshold).astype(int)
    out = {
        "auroc": float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else float("nan"),
        "auprc": float(average_precision_score(y, p)) if len(np.unique(y)) == 2 else float("nan"),
        "balanced_accuracy": float(balanced_accuracy_score(y, yhat)),
        "f1": float(f1_score(y, yhat, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, yhat)) if len(np.unique(yhat)) > 1 else 0.0,
        "brier": float(brier_score_loss(y, p)),
        "ece": expected_calibration_error(y, p),
        "n": int(len(y)),
        "prevalence": float(y.mean()),
    }
    return out


def bootstrap_ci(y: np.ndarray, p: np.ndarray, metric: Callable = roc_auc_score, n_boot: int = 1000, alpha: float = 0.05, seed: int = 0) -> Tuple[float, float, float]:
    """Point estimate and percentile CI by stratified bootstrap over samples."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    p = np.asarray(p)
    pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
    vals = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        try:
            vals.append(metric(y[idx], p[idx]))
        except ValueError:
            continue
    est = metric(y, p)
    lo, hi = np.percentile(vals, [100 * alpha / 2, 100 * (1 - alpha / 2)]) if vals else (np.nan, np.nan)
    return float(est), float(lo), float(hi)
