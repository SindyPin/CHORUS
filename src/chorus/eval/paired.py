"""Within-cohort evaluation of several omics layers measured on the same people (Test 2).

For one paired cohort (every sample has every layer) and one binary contrast, every model is
fitted on the same repeated stratified K fold splits, so that per fold metrics are comparable
across models:

  layer:<name>   one layer alone: per fold, features standardised on the training fold, the
                 ``k`` features with the largest absolute training effect (DifferentialSelector)
                 kept, then the classifier
  early_fusion   the per layer selections of the fold concatenated; each layer's features are
                 standardised on the training fold and divided by the square root of the number of
                 that layer's selected features, so each layer carries the same total variance
  late_fusion    chorus.integrate.mosaic.LateFusion: the meta learner is fitted on inner out of fold
                 probabilities of the layer models (inner stratified K fold on the training fold),
                 then applied to the probabilities of the layer models refitted on the whole
                 training fold
  mofa           MOFA+ (mofapy2) factors learned on the training fold (one group, layers as views,
                 all features of the fold), test samples projected onto the training loadings
                 (chorus.integrate.mosaic.mofa_factors), then the classifier on the factors
  diablo         mixOmics block.splsda (R) on the training fold, ``keepX`` features per block and
                 component, full design 0.1; the score is the positive class column of the
                 averaged weighted prediction (predict(...)$WeightedPredict), clipped to [0, 1]; it
                 is not a calibrated probability

Missing values are filled inside the fold from the training samples of that fold (``impute``
per layer: ``min`` or ``median``). Features observed in fewer than ``min_observed`` of the
cohort's samples are dropped before the splits (label free).

Cross layer agreement:
  effect_correlation       per pair of layers, Spearman correlation over the shared features of
                           the Hedges g of the contrast computed in each layer on all samples
  probability_correlation  per repeat, Spearman and Pearson correlation across subjects of the
                           out of fold probabilities of two models (each subject is in exactly one
                           test fold per repeat)
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold

from ..features import make_selector
from ..features.invariant import per_study_effects
from ..integrate.mosaic import LateFusion, mofa_factors
from ..models import make_classifier
from .metrics import binary_metrics

R_LIB = os.environ.get("CHORUS_R_LIB", str(Path("~/Library/R/chorus-mixomics").expanduser()))

DIABLO_R = r"""
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]; lib <- args[2]; ncomp <- as.integer(args[3]); keepx <- as.integer(args[4])
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(mixOmics))
blocks <- strsplit(readLines(file.path(dir, "blocks.txt")), "\t")[[1]]
folds <- read.csv(file.path(dir, "folds.csv"), stringsAsFactors = FALSE)
y <- read.csv(file.path(dir, "y.csv"), stringsAsFactors = FALSE)
X <- lapply(blocks, function(b) as.matrix(read.csv(file.path(dir, paste0(b, ".csv")), row.names = 1, check.names = FALSE)))
names(X) <- blocks
out <- list()
for (f in unique(folds$fold_id)) {
  te <- folds$sample[folds$fold_id == f]
  tr <- setdiff(y$sample, te)
  Xtr <- lapply(X, function(m) m[tr, , drop = FALSE]); Xte <- lapply(X, function(m) m[te, , drop = FALSE])
  ytr <- factor(y$y[match(tr, y$sample)], levels = c(0, 1))
  design <- matrix(0.1, length(blocks), length(blocks), dimnames = list(blocks, blocks)); diag(design) <- 0
  kx <- lapply(Xtr, function(m) rep(min(keepx, ncol(m)), ncomp))
  res <- tryCatch({
    fit <- block.splsda(Xtr, ytr, ncomp = ncomp, keepX = kx, design = design, scale = TRUE)
    pr <- predict(fit, Xte)
    wp <- pr$WeightedPredict[, "1", ncomp]
    data.frame(fold_id = f, sample = te, p = as.numeric(wp), error = NA_character_)
  }, error = function(e) data.frame(fold_id = f, sample = te, p = NA_real_, error = conditionMessage(e)))
  out[[length(out) + 1]] <- res
}
write.csv(do.call(rbind, out), file.path(dir, "pred.csv"), row.names = FALSE)
"""


def diablo_available(r_lib: str = R_LIB) -> bool:
    """True when Rscript runs and mixOmics loads from ``r_lib`` or the default libraries."""
    if shutil.which("Rscript") is None:
        return False
    code = f'.libPaths(c("{r_lib}", .libPaths())); quit(status = if (requireNamespace("mixOmics", quietly = TRUE)) 0 else 1)'
    return subprocess.run(["Rscript", "-e", code], capture_output=True).returncode == 0


def mofa_available() -> bool:
    try:
        import mofapy2  # noqa: F401
    except ImportError:
        return False
    return True


def _fill(train: pd.DataFrame, other: pd.DataFrame, how: str) -> tuple:
    ref = train.min(axis=0) if how == "min" else train.median(axis=0)
    return train.fillna(ref), other.fillna(ref)


def _standardise(train: pd.DataFrame, other: pd.DataFrame) -> tuple:
    mu, sd = train.mean(axis=0), train.std(axis=0, ddof=1)
    sd = sd.where(sd > 0, 1.0)
    return (train - mu) / sd, (other - mu) / sd


class PairedCV:
    """Repeated stratified K fold over the subjects of one paired cohort, shared by every model."""

    def __init__(
        self,
        views: Dict[str, pd.DataFrame],
        y: pd.Series,
        impute: Optional[Dict[str, str]] = None,
        k: int = 50,
        classifier: str = "logreg",
        n_splits: int = 5,
        n_repeats: int = 10,
        inner_splits: int = 5,
        min_observed: float = 0.5,
        mofa_factors_n: int = 10,
        diablo_ncomp: int = 2,
        diablo_keepx: int = 20,
        random_state: int = 0,
        r_lib: str = R_LIB,
    ):
        y = y.dropna().astype(int)
        common = y.index
        for V in views.values():
            common = common.intersection(V.index)
        if len(common) < len(y):
            raise ValueError(f"{len(y) - len(common)} labelled samples lack a layer; pass paired samples only")
        self.y = y.loc[common]
        self.views = {}
        self.dropped_features = {}
        for name, V in views.items():
            V = V.loc[common].apply(pd.to_numeric, errors="coerce")
            keep = V.notna().mean(axis=0) >= min_observed
            self.dropped_features[name] = int((~keep).sum())
            self.views[name] = V.loc[:, keep]
        self.impute = {n: (impute or {}).get(n, "median") for n in self.views}
        self.k, self.classifier, self.n_splits, self.n_repeats = k, classifier, n_splits, n_repeats
        self.inner_splits, self.min_observed, self.random_state = inner_splits, min_observed, random_state
        self.mofa_factors_n, self.diablo_ncomp, self.diablo_keepx, self.r_lib = mofa_factors_n, diablo_ncomp, diablo_keepx, r_lib

    # ----------------------------------------------------------------------------- folds
    def folds(self) -> pd.DataFrame:
        """One row per subject, repeat and fold: the test fold of each subject in each repeat.
        Repeat r uses the split seed ``random_state + r``."""
        rows = []
        for r in range(self.n_repeats):
            skf = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=self.random_state + r)
            for f, (_, te) in enumerate(skf.split(np.zeros(len(self.y)), self.y.values)):
                rows += [{"sample": s, "repeat": r, "fold": f} for s in self.y.index[te]]
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------- per fold pieces
    def _layer_fold(self, name: str, tr: pd.Index, te: pd.Index) -> tuple:
        """Imputed, standardised and selected training and test blocks of one layer."""
        V = self.views[name]
        a, b = _fill(V.loc[tr], V.loc[te], self.impute[name])
        a, b = _standardise(a, b)
        sel = make_selector("differential", k=min(self.k, a.shape[1])).fit(a, self.y.loc[tr], pd.Series("s", index=tr))
        cols = list(sel.selected_)
        return a[cols], b[cols]

    def _fit_predict(self, Xa: pd.DataFrame, ya: pd.Series, Xb: pd.DataFrame) -> np.ndarray:
        clf = make_classifier(self.classifier, self.random_state).fit(Xa, ya)
        return clf.predict_proba(Xb)[:, 1]

    def _layer_probs(self, tr: pd.Index, te: pd.Index) -> Dict[str, np.ndarray]:
        out = {}
        for name in self.views:
            a, b = self._layer_fold(name, tr, te)
            out[name] = self._fit_predict(a, self.y.loc[tr], b)
        return out

    def _late_fusion(self, tr: pd.Index, te: pd.Index, test_probs: Dict[str, np.ndarray]) -> np.ndarray:
        ytr = self.y.loc[tr]
        inner = StratifiedKFold(n_splits=self.inner_splits, shuffle=True, random_state=self.random_state)
        rows = []
        for itr, ite in inner.split(np.zeros(len(tr)), ytr.values):
            p = self._layer_probs(tr[itr], tr[ite])
            for name, v in p.items():
                rows.append(pd.DataFrame({"subject_id": tr[ite], "layer": name, "p": v, "y": ytr.values[ite]}))
        lf = LateFusion(list(self.views)).fit(pd.concat(rows, ignore_index=True))
        test = pd.concat([pd.DataFrame({"subject_id": te, "layer": n, "p": v}) for n, v in test_probs.items()], ignore_index=True)
        return lf.predict_proba(test).reindex(te).values

    def _mofa(self, tr: pd.Index, te: pd.Index) -> np.ndarray:
        views = {}
        for name, V in self.views.items():
            a, b = _fill(V.loc[tr], V.loc[te], self.impute[name])
            views[name] = pd.concat([a, b]).add_prefix(f"{name}|")  # MOFA+ needs names unique across views
        idx = tr.append(te)
        Z = mofa_factors(views, pd.Series("cohort", index=idx), n_factors=self.mofa_factors_n, seed=self.random_state, train_mask=pd.Series(idx.isin(tr), index=idx))
        Z = Z.astype(float)
        return self._fit_predict(Z.loc[tr], self.y.loc[tr], Z.loc[te])

    # --------------------------------------------------------------------------- DIABLO
    def _diablo(self, folds: pd.DataFrame) -> pd.DataFrame:
        """All folds in one R session. Imputation uses the training samples of each fold, so the
        imputed blocks are written per fold; scaling is done by block.splsda on the training fold."""
        d = Path(tempfile.mkdtemp(prefix="chorus_diablo_"))
        try:
            preds = []
            for (r, f), g in folds.groupby(["repeat", "fold"]):
                te = pd.Index(g["sample"])
                tr = self.y.index.difference(te)
                fd = d / f"r{r}_f{f}"
                fd.mkdir()
                for name, V in self.views.items():
                    a, b = _fill(V.loc[tr], V.loc[te], self.impute[name])
                    pd.concat([a, b]).loc[self.y.index].to_csv(fd / f"{name}.csv")
                (fd / "blocks.txt").write_text("\t".join(self.views))
                pd.DataFrame({"sample": self.y.index, "y": self.y.values}).to_csv(fd / "y.csv", index=False)
                pd.DataFrame({"sample": te, "fold_id": 0}).to_csv(fd / "folds.csv", index=False)
                (d / "diablo.R").write_text(DIABLO_R)
                res = subprocess.run(["Rscript", str(d / "diablo.R"), str(fd), self.r_lib, str(self.diablo_ncomp), str(self.diablo_keepx)], capture_output=True, text=True)
                if res.returncode != 0:
                    preds.append(pd.DataFrame({"sample": te, "repeat": r, "fold": f, "p": np.nan, "error": res.stderr[-200:]}))
                    continue
                p = pd.read_csv(fd / "pred.csv")
                preds.append(pd.DataFrame({"sample": p["sample"], "repeat": r, "fold": f, "p": p["p"].clip(0, 1), "error": p["error"]}))
            return pd.concat(preds, ignore_index=True)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    # ------------------------------------------------------------------------------ run
    def run(self, models: Optional[List[str]] = None, verbose: bool = True) -> tuple:
        """Returns (per fold metrics, out of fold predictions). ``models`` defaults to every
        layer, early_fusion, late_fusion, and mofa and diablo when their dependencies load."""
        layer_models = [f"layer:{n}" for n in self.views]
        if models is None:
            models = layer_models + ["early_fusion", "late_fusion"] + (["mofa"] if mofa_available() else []) + (["diablo"] if diablo_available(self.r_lib) else [])
        folds = self.folds()
        preds = []
        for (r, f), g in folds.groupby(["repeat", "fold"]):
            te = pd.Index(g["sample"])
            tr = self.y.index.difference(te)
            blocks = {n: self._layer_fold(n, tr, te) for n in self.views}
            layer_p = {n: self._fit_predict(a, self.y.loc[tr], b) for n, (a, b) in blocks.items()}
            fold_p: Dict[str, np.ndarray] = {}
            for n in self.views:
                fold_p[f"layer:{n}"] = layer_p[n]
            if "early_fusion" in models:
                A = pd.concat([a / np.sqrt(a.shape[1]) for a, _ in blocks.values()], axis=1, keys=list(blocks))
                B = pd.concat([b / np.sqrt(b.shape[1]) for _, b in blocks.values()], axis=1, keys=list(blocks))
                A.columns = [f"{l}|{c}" for l, c in A.columns]
                B.columns = A.columns
                fold_p["early_fusion"] = self._fit_predict(A, self.y.loc[tr], B)
            if "late_fusion" in models:
                fold_p["late_fusion"] = self._late_fusion(tr, te, layer_p)
            if "mofa" in models:
                try:
                    fold_p["mofa"] = self._mofa(tr, te)
                except Exception as e:  # record and continue
                    preds.append(pd.DataFrame({"sample": te, "repeat": r, "fold": f, "model": "mofa", "y": self.y.loc[te].values, "p": np.nan, "error": str(e)[:200]}))
            for m, p in fold_p.items():
                if m in models:
                    preds.append(pd.DataFrame({"sample": te, "repeat": r, "fold": f, "model": m, "y": self.y.loc[te].values, "p": p, "error": None}))
            if verbose:
                print(f"repeat {r} fold {f}: " + ", ".join(f"{m} {binary_metrics(self.y.loc[te].values, p)['auroc']:.3f}" for m, p in fold_p.items() if m in models), flush=True)
        if "diablo" in models:
            dp = self._diablo(folds)
            dp["model"] = "diablo"
            dp["y"] = self.y.loc[dp["sample"]].values
            preds.append(dp)
        P = pd.concat(preds, ignore_index=True)
        return fold_metrics(P), P


def fold_metrics(P: pd.DataFrame) -> pd.DataFrame:
    """Metrics of each model on each test fold; failed folds keep their error message."""
    rows = []
    for (m, r, f), g in P.groupby(["model", "repeat", "fold"], sort=False):
        err = g["error"].dropna()
        if g["p"].isna().any():
            rows.append({"model": m, "repeat": r, "fold": f, "n": len(g), "error": err.iloc[0] if len(err) else "missing prediction"})
            continue
        rows.append({"model": m, "repeat": r, "fold": f, "n": len(g), "prevalence": float(g["y"].mean()), **binary_metrics(g["y"].values, g["p"].values), "error": None})
    return pd.DataFrame(rows)


def summarise_folds(F: pd.DataFrame) -> pd.DataFrame:
    """Mean and standard deviation over folds of each metric, per model."""
    ok = F[F["error"].isna()]
    cols = [c for c in ["auroc", "auprc", "balanced_accuracy", "f1", "mcc", "brier", "ece"] if c in ok.columns]
    out = ok.groupby("model", sort=False)[cols].mean()
    out["auroc_sd"] = ok.groupby("model", sort=False)["auroc"].std(ddof=1)
    out["n_folds"] = ok.groupby("model", sort=False).size()
    out["n_failed_folds"] = F.groupby("model", sort=False)["error"].apply(lambda e: int(e.notna().sum())).reindex(out.index)
    return out.reset_index()


def effect_correlation(views: Dict[str, pd.DataFrame], y: pd.Series) -> tuple:
    """Hedges g of the contrast per feature in each layer (all paired samples), and per pair of
    layers the Spearman correlation over their shared features. Returns (effects, correlations)."""
    y = y.dropna()
    eff = {}
    for name, V in views.items():
        V = V.loc[y.index].apply(pd.to_numeric, errors="coerce")
        e = per_study_effects(V, y, pd.Series("s", index=y.index))
        eff[name] = e["s"]["g"] if "s" in e else pd.Series(np.nan, index=V.columns)
    E = pd.DataFrame(eff)
    return E, layer_pairs(E)


def layer_pairs(E: pd.DataFrame, min_shared: int = 10) -> pd.DataFrame:
    """Per pair of layers (columns of ``E``) the Spearman correlation over their shared features,
    with ``fewer_than_<min_shared>_shared`` set for pairs below ``min_shared`` shared features."""
    rows = []
    names = list(E.columns)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            ok = E[a].notna() & E[b].notna()
            rho, p = stats.spearmanr(E.loc[ok, a], E.loc[ok, b]) if ok.sum() > 2 else (np.nan, np.nan)
            rows.append({"layer_a": a, "layer_b": b, "n_shared": int(ok.sum()), "spearman": float(rho), "p": float(p),
                         f"fewer_than_{min_shared}_shared": bool(ok.sum() < min_shared)})
    return pd.DataFrame(rows)


def probability_correlation(P: pd.DataFrame) -> pd.DataFrame:
    """Per repeat and pair of models, Spearman and Pearson correlation across subjects of the
    out of fold probabilities (subjects with a prediction from both models)."""
    rows = []
    models = list(pd.unique(P["model"]))
    for r, g in P.groupby("repeat"):
        W = g.pivot_table(index="sample", columns="model", values="p")
        for i, a in enumerate(models):
            for b in models[i + 1 :]:
                if a not in W or b not in W:
                    continue
                ok = W[a].notna() & W[b].notna()
                if ok.sum() < 3:
                    continue
                rows.append(
                    {
                        "repeat": r,
                        "model_a": a,
                        "model_b": b,
                        "n": int(ok.sum()),
                        "spearman": float(stats.spearmanr(W.loc[ok, a], W.loc[ok, b])[0]),
                        "pearson": float(np.corrcoef(W.loc[ok, a], W.loc[ok, b])[0, 1]),
                    }
                )
    return pd.DataFrame(rows)
