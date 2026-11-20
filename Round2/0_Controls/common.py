"""Shared code of the Round2 known answer controls (Round2/PROTOCOL.md with Amendment 1).

Controls:
  sex  classification, recorded sex (female = 1, male = 0); every participant of a cohort, whatever the LC status
  age  regression, recorded age in years

Used by the step scripts 01_one_cohort, 02_horizontal_one_plus_one, 03_vertical_one_plus_one and 04_summary. It
holds: cohort loading (Round1 loader), one sample per participant, the choice of anchor and partner, the folds,
preprocessing and feature selection inside each training fold, the models with their fixed grids, the failure
wrapper (failures.csv), the metrics, the corrected resampled t test, the resumable runner and the timing pass.
"""

from __future__ import annotations

import glob
import hashlib
import json
import os
import signal
import subprocess
import tempfile
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

CODE = Path(__file__).resolve().parents[2]
CTRL = Path(__file__).resolve().parent
CACHE = CTRL / "cache"  # arrays per control and layer (data: not committed)
SEED = 20261120
N_SPLITS, N_REPEATS, INNER = 5, 10, 3  # the 10 seeded repeats of PROTOCOL.md generate the folds
N_REPEATS_USED = 3  # Amendment 3, item 1: the first 3 repeats are run and summarised
N_FOLDS = N_SPLITS * N_REPEATS_USED  # folds 0 to 14 (RepeatedStratifiedKFold lists repeat by repeat)
K_FEATURES, K_VERTICAL = 200, 100
FIT_LIMIT_S = 15 * 60
BUDGET_H = 12.0
WORKERS = int(os.environ.get("ROUND2_WORKERS", "10"))
MIN_PER_CLASS = 10
MIN_AGE_PARTICIPANTS = 20  # regression: "10 per class" read as at least 20 participants with the outcome
os.environ.setdefault("CHORUS_DISABLE_TABPFN", "1")  # the Round1 TACO ensemble without TabPFN (four members)

LAYER = {"bulk_rna": "transcriptomics", "olink": "proteomics", "somascan": "proteomics", "ms_proteomics": "proteomics",
         "metabolomics": "metabolomics"}
# not used as anchor or partner (PROTOCOL.md, Amendment 1 item 1), plus GSE300129: Round1 verification failed because
# its recorded sex agrees with sex inferred from expression no better than chance
EXCLUDED = {"GSE270045": "flagged: labels confounded with collection", "GSE228320": "flagged: non standard case definition",
            "CORSIP2026": "flagged: gene panel only", "Buggert2025": "flagged: split by country (Amendment 1, item 1)",
            "GSE300129": "Round1 verification failed: recorded sex no better than chance against inferred sex"}
SPLITS = {"Buggert2025": ("collection_country", {"United Kingdom": "Buggert2025_UK", "Sweden": "Buggert2025_SE"})}
PARENT = {k: p for p, (_, kids) in SPLITS.items() for k in kids.values()}
# Amendment 1, item 5 (reading recorded in 0_Labels/README.md): same pipeline and reference = same scale
SAME_PIPELINE = {"salmon 1.10.3, GENCODE v47 index (Round1 re-quantification)":
                 {"GSE224615", "PRJNA1149196_SA", "PRJNA1149196_QLD", "PRJNA1184005"}}

CLASSIFICATION_METRICS = ["AUROC", "AUPRC", "balanced_accuracy", "F1", "MCC", "sensitivity", "specificity", "Brier", "log_loss", "ECE"]
REGRESSION_METRICS = ["MAE", "RMSE", "R2", "Spearman"]
LOWER_IS_BETTER = {"Brier", "log_loss", "ECE", "MAE", "RMSE"}
PRIMARY = {"classification": "AUROC", "regression": "MAE"}
TASK = {"sex": "classification", "age": "regression"}


# ----------------------------------------------------------------------------------------------- cohorts

# On disk with a matrix and the control variables but in no Round1 config (Round1 used them for no task): added here
# so that the controls choose among every cohort that records the variable. Captur2022 records neither sex nor age;
# GSE300129 is excluded above.
ON_DISK = {
    "Holmes2021": dict(study_id="Holmes2021", matrix="data/metabolomics/Holmes2021/Holmes2021_ids_matrix.csv",
                       meta="data/meta/Holmes2021_meta.csv", feature_space="metabolite_id"),
    "Zoodsma2022_Nijmegen": dict(study_id="Zoodsma2022_Nijmegen", matrix="data/proteomics/Zoodsma2022_Nijmegen/Zoodsma2022_Nijmegen_matrix.csv",
                                 meta="data/meta/Zoodsma2022_Nijmegen_meta.csv", feature_space="uniprot"),
}


def config_entries() -> dict:
    seen = {}
    for f in sorted(glob.glob(str(CODE / "Round1/configs/*.yaml"))):
        for c in yaml.safe_load(open(f))["cohorts"]:
            if (CODE / c["matrix"]).exists():
                seen.setdefault(c["study_id"], c)
    for k, c in ON_DISK.items():
        seen.setdefault(k, c)
    return seen


def matrix_columns(path: Path) -> list:
    import gzip

    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        head = f.readline().rstrip("\n")
    sep = "\t" if (".tsv" in path.name or ".txt" in path.name) else ","
    return [c.strip('"') for c in head.split(sep)[1:]]


def read_meta(study: str) -> pd.DataFrame:
    parent = PARENT.get(study, study)
    m = pd.read_csv(CODE / f"data/meta/{parent}_meta.csv", low_memory=False)
    if study in PARENT:
        col, kids = SPLITS[parent]
        val = {v: k for k, v in kids.items()}[study]
        m = m[m[col].astype(str) == val]
    return m


def outcome(m: pd.DataFrame, control: str) -> pd.Series:
    if control == "sex":
        s = m["sex"].astype(str).str.lower() if "sex" in m else pd.Series(np.nan, index=m.index)
        return s.map({"female": 1.0, "male": 0.0})
    return pd.to_numeric(m["age"], errors="coerce") if "age" in m else pd.Series(np.nan, index=m.index)


def one_per_participant(m: pd.DataFrame) -> pd.DataFrame:
    """The earliest sample of each participant: smallest days_since_infection, then numeric timepoint, then sample_id."""
    key = pd.DataFrame({"d": pd.to_numeric(m.get("days_since_infection"), errors="coerce") if "days_since_infection" in m else np.nan,
                        "t": pd.to_numeric(m.get("timepoint"), errors="coerce") if "timepoint" in m else np.nan,
                        "s": m["sample_id"].astype(str)}, index=m.index)
    order = key.sort_values(["d", "t", "s"], na_position="last").index
    return m.loc[order].drop_duplicates("subject_id", keep="first")


def candidates(control: str) -> pd.DataFrame:
    """Every cohort that records the control variable: participants with it (one sample each, sample in the matrix)."""
    E = config_entries()
    studies = sorted(set(E) | set(PARENT))
    rows = []
    for s in studies:
        parent = PARENT.get(s, s)
        if parent not in E or not (CODE / f"data/meta/{parent}_meta.csv").exists():
            continue
        m = read_meta(s)
        if "omics_layer" not in m or m["omics_layer"].dropna().empty:
            continue
        layer = LAYER.get(str(m["omics_layer"].dropna().iloc[0]))
        if layer is None:
            continue
        cols = set(matrix_columns(CODE / E[parent]["matrix"]))
        m = m[m["sample_id"].astype(str).isin(cols)]
        y = outcome(m, control)
        m = one_per_participant(m[y.notna()])
        y = outcome(m, control)
        if len(m) == 0:
            continue
        r = dict(control=control, layer=layer, study_id=s, n=len(m), excluded=EXCLUDED.get(s, ""))
        if control == "sex":
            r.update(n_female=int((y == 1).sum()), n_male=int((y == 0).sum()))
            r["eligible"] = r["n_female"] >= MIN_PER_CLASS and r["n_male"] >= MIN_PER_CLASS and not r["excluded"]
        else:
            r.update(age_min=float(y.min()), age_median=float(y.median()), age_max=float(y.max()))
            r["eligible"] = len(m) >= MIN_AGE_PARTICIPANTS and not r["excluded"]
        rows.append(r)
    return pd.DataFrame(rows)


def choose(C: pd.DataFrame) -> pd.DataFrame:
    """Anchor = largest eligible cohort of the layer; partner = next largest (PROTOCOL.md); anchor sensitivity
    (Amendment 1, item 3): the three largest eligible within fewer than five participants."""
    rows = []
    for (control, layer), d in C.groupby(["control", "layer"], sort=True):
        ok = d[d["eligible"]].sort_values(["n", "study_id"], ascending=[False, True])
        top = ok.head(3)
        sens = list(top["study_id"]) if len(top) > 1 and top["n"].max() - top["n"].min() < 5 else []
        rows.append(dict(control=control, layer=layer, anchor=ok["study_id"].iloc[0] if len(ok) else "",
                         partner=ok["study_id"].iloc[1] if len(ok) > 1 else "",
                         anchor_n=int(ok["n"].iloc[0]) if len(ok) else 0, partner_n=int(ok["n"].iloc[1]) if len(ok) > 1 else 0,
                         sensitivity_anchors="; ".join(sens), same_scale=same_scale(*(list(ok["study_id"].iloc[:2]) + ["", ""])[:2])))
    return pd.DataFrame(rows)


def same_scale(a: str, b: str) -> bool:
    return any(a in g and b in g for g in SAME_PIPELINE.values())


_LOADED: dict = {}


def load(study: str):
    """(X, meta, is_counts) through the Round1 loader; one sample per participant comes later."""
    if study in _LOADED:
        return _LOADED[study]
    from chorus.io import Cohort

    parent = PARENT.get(study, study)
    c = config_entries()[parent]
    m = read_meta(study)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        co = Cohort.from_files(parent, str(m["omics_layer"].dropna().iloc[0]), str(CODE / c["matrix"]), str(CODE / c["meta"]),
                               features_in_rows=c.get("features_in_rows", True), is_counts=c.get("is_counts", False),
                               feature_space=c.get("feature_space", "ensembl_gene"))
    X = co.X.loc[:, co.X.notna().any(axis=0)]
    X.index = X.index.astype(str)
    m = m[m["sample_id"].astype(str).isin(X.index)].copy()
    _LOADED[study] = (X, m, bool(c.get("is_counts", False)), c["matrix"])
    return _LOADED[study]


def transform(X: pd.DataFrame, is_counts: bool) -> pd.DataFrame:
    """RNA counts: log2 CPM (per sample, so no fold statistics are needed). Other layers: as loaded (log scale: the
    Round1 loader applies log2(x + 1) to linear scale matrices)."""
    if not is_counts:
        return X.astype(float)
    lib = X.sum(axis=1).replace(0, np.nan)
    return np.log2(X.div(lib, axis=0) * 1e6 + 1)


def control_data(study: str, control: str, samples: list | None = None):
    """(X transformed, y, meta) for one cohort and control: one sample per participant with the outcome recorded."""
    X, m, is_counts, _ = load(study)
    if samples is not None:
        m = m[m["sample_id"].astype(str).isin(samples)]
    y = outcome(m, control)
    m = one_per_participant(m[y.notna()])
    ids = m["sample_id"].astype(str).tolist()
    return transform(X.loc[ids], is_counts), outcome(m, control).to_numpy(), m


def folds(y: np.ndarray, task: str) -> list:
    from sklearn.model_selection import RepeatedStratifiedKFold

    strat = y if task == "classification" else pd.qcut(pd.Series(y).rank(method="first"), 5, labels=False).to_numpy()
    rk = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=SEED)
    return [(tr, te) for tr, te in rk.split(np.zeros(len(y)), strat)]


def file_sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# --------------------------------------------------------------------------------------- preprocessing

def zscore_within(X: np.ndarray, ref: np.ndarray | None = None) -> np.ndarray:
    """z score each feature with the statistics of ``ref`` (default: X itself); no labels used."""
    ref = X if ref is None else ref
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mu, sd = np.nanmean(ref, axis=0), np.nanstd(ref, axis=0)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
    return (X - np.nan_to_num(mu)) / sd


def select(Xtr: np.ndarray, ytr: np.ndarray, task: str, k: int) -> np.ndarray:
    """Indices of the top k features: |Welch t| (classification) or |Spearman| (regression), training fold only."""
    if task == "classification":
        a, b = Xtr[ytr == 1], Xtr[ytr == 0]
        with np.errstate(all="ignore"):
            t = (a.mean(0) - b.mean(0)) / np.sqrt(a.var(0, ddof=1) / len(a) + b.var(0, ddof=1) / len(b))
    else:
        R = stats.rankdata(Xtr, axis=0)
        r = stats.rankdata(ytr)
        Rc, rc = R - R.mean(0), r - r.mean()
        with np.errstate(all="ignore"):
            t = (Rc * rc[:, None]).sum(0) / np.sqrt((Rc ** 2).sum(0) * (rc ** 2).sum())
    s = np.nan_to_num(np.abs(t), nan=-1.0)
    k = min(k, int((s >= 0).sum()) or len(s))
    return np.argsort(-s, kind="stable")[:k]


def prepare(Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, task: str, k: int):
    """Training fold median imputation, feature selection, standardisation with training fold statistics."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        med = np.nanmedian(Xtr, axis=0)
    keep = np.isfinite(med)
    Xtr, Xte, med = Xtr[:, keep], Xte[:, keep], med[keep]
    Xtr = np.where(np.isnan(Xtr), med, Xtr)
    Xte = np.where(np.isnan(Xte), med, Xte)
    idx = select(Xtr, ytr, task, k)
    Xtr, Xte = Xtr[:, idx], Xte[:, idx]
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    return (Xtr - mu) / sd, (Xte - mu) / sd, np.flatnonzero(keep)[idx]


# ---------------------------------------------------------------------------------------------- models

CLASSICAL = ["logistic_regression", "elastic_net", "linear_svm", "rbf_svm", "knn", "gaussian_nb", "random_forest",
             "extra_trees", "hist_gradient_boosting", "xgboost", "lightgbm", "catboost", "mlp"]
ENSEMBLE = ["taco_ensemble"]
FOUNDATION = ["tabpfn_v2", "tabpfn_2_5", "tabpfn_3", "tabpfn_3_5", "tabpfn_3_5_fast", "tabicl_v2", "tabdpt", "mitra"]
# TabPFN versions, each pinned explicitly (tabpfn 9.2.0; Amendment 5, A). The checkpoint loaded is written to every
# result row of the TabPFN 3 and 3.5 jobs; the v2 and 2.5 checkpoints are listed in Round2/environment.txt.
TABPFN_VERSIONS = {"tabpfn_v2": "V2", "tabpfn_2_5": "V2_5", "tabpfn_3": "V3", "tabpfn_3_5": "V3_5", "tabpfn_3_5_fast": "V3_5_FAST"}
EXPRESSION = ["bulkformer", "bulkrnabert"]  # gene expression only: embeddings, then (logistic or ridge) regression
# Amendment 2, item 1: removed from the model grid on this hardware. Its jobs are still listed, so that every task and
# arm gets a failure row; they are never run. The step 01 fold 0 fits that measured it are kept in
# 01_one_cohort/output/mitra_timed/ (out of raw/, so no summary treats them as results).
REMOVED = {"mitra": "time budget exceeded on this hardware"}
MITRA_TIMED = CTRL / "01_one_cohort" / "output" / "mitra_timed"


def models_for(layer: str) -> list:
    return CLASSICAL + ENSEMBLE + FOUNDATION + (EXPRESSION if layer == "transcriptomics" else [])


def classical(name: str, task: str):
    """(estimator, grid) of the fixed grids in PROTOCOL.md; regression forms (Amendment 1, item 7) with analogous grids."""
    from sklearn import ensemble as E, linear_model as L, neighbors as N, neural_network as NN, svm as S
    from sklearn.naive_bayes import GaussianNB

    rs = SEED
    c = task == "classification"
    if name == "logistic_regression":
        return (L.LogisticRegression(max_iter=5000), {"C": [0.01, 0.1, 1, 10]}) if c else (L.Ridge(), {"alpha": [0.1, 1, 10, 100]})
    if name == "elastic_net":
        if c:
            return L.LogisticRegression(solver="saga", l1_ratio=0.5, max_iter=5000, random_state=rs), {"C": [0.01, 0.1, 1], "l1_ratio": [0.2, 0.5, 0.8]}
        return L.ElasticNet(max_iter=5000, random_state=rs), {"alpha": [0.01, 0.1, 1], "l1_ratio": [0.2, 0.5, 0.8]}
    if name == "linear_svm":
        return (S.SVC(kernel="linear", probability=True, random_state=rs), {"C": [0.01, 0.1, 1, 10]}) if c else (S.SVR(kernel="linear"), {"C": [0.01, 0.1, 1, 10]})
    if name == "rbf_svm":
        g = {"C": [0.1, 1, 10], "gamma": ["scale", 0.001, 0.01]}
        return (S.SVC(kernel="rbf", probability=True, random_state=rs), g) if c else (S.SVR(kernel="rbf"), g)
    if name == "knn":
        g = {"n_neighbors": [3, 5, 11, 21], "weights": ["uniform", "distance"]}
        return (N.KNeighborsClassifier(), g) if c else (N.KNeighborsRegressor(), g)
    if name == "gaussian_nb":
        if not c:
            raise NotImplementedError("Gaussian naive Bayes has no regression form")
        return GaussianNB(), {"var_smoothing": [1e-9, 1e-6, 1e-3]}
    if name in ("random_forest", "extra_trees"):
        g = {"max_features": ["sqrt", 0.3], "min_samples_leaf": [1, 5]}
        cls = {("random_forest", True): E.RandomForestClassifier, ("random_forest", False): E.RandomForestRegressor,
               ("extra_trees", True): E.ExtraTreesClassifier, ("extra_trees", False): E.ExtraTreesRegressor}[(name, c)]
        return cls(n_estimators=500, random_state=rs, n_jobs=1), g
    if name == "hist_gradient_boosting":
        g = {"learning_rate": [0.05, 0.1], "max_depth": [3, None]}
        return (E.HistGradientBoostingClassifier(random_state=rs), g) if c else (E.HistGradientBoostingRegressor(random_state=rs), g)
    if name == "xgboost":
        import xgboost as X

        g = {"max_depth": [2, 4], "learning_rate": [0.05, 0.1]}
        kw = dict(n_estimators=300, random_state=rs, n_jobs=1, verbosity=0)
        return (X.XGBClassifier(eval_metric="logloss", **kw), g) if c else (X.XGBRegressor(**kw), g)
    if name == "lightgbm":
        import lightgbm as G

        g = {"num_leaves": [7, 15], "learning_rate": [0.05, 0.1]}
        kw = dict(n_estimators=300, min_child_samples=5, random_state=rs, n_jobs=1, verbose=-1)
        return (G.LGBMClassifier(**kw), g) if c else (G.LGBMRegressor(**kw), g)
    if name == "catboost":
        import catboost as CB

        g = {"depth": [4, 6], "learning_rate": [0.05, 0.1]}
        kw = dict(iterations=300, random_seed=rs, verbose=0, thread_count=1, allow_writing_files=False)
        return (CB.CatBoostClassifier(**kw), g) if c else (CB.CatBoostRegressor(**kw), g)
    if name == "mlp":
        g = {"hidden_layer_sizes": [(64,), (64, 32)], "alpha": [1e-4, 1e-2]}
        return (NN.MLPClassifier(max_iter=1000, random_state=rs), g) if c else (NN.MLPRegressor(max_iter=1000, random_state=rs), g)
    raise KeyError(name)


def taco(task: str):
    """The Round1 TACO ensemble: random forest, extra trees, RBF SVC and MLP, soft voting (chorus.models.classifiers,
    without TabPFN); for regression the same members as regressors with the same settings, averaged."""
    if task == "classification":
        from chorus.models.classifiers import make_classifier

        return make_classifier("taco", random_state=SEED)
    from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor, VotingRegressor
    from sklearn.neural_network import MLPRegressor
    from sklearn.svm import SVR

    return VotingRegressor([("rf", RandomForestRegressor(n_estimators=500, random_state=SEED, n_jobs=1)),
                            ("extratrees", ExtraTreesRegressor(n_estimators=500, random_state=SEED, n_jobs=1)),
                            ("svc", SVR(kernel="rbf", C=1.0)),
                            ("mlp", MLPRegressor(hidden_layer_sizes=(64,), alpha=1e-2, max_iter=1000, random_state=SEED))])


def foundation(name: str, task: str):
    c = task == "classification"
    if name in TABPFN_VERSIONS:
        from tabpfn import TabPFNClassifier, TabPFNRegressor
        from tabpfn.constants import ModelVersion

        cls = TabPFNClassifier if c else TabPFNRegressor
        return cls.create_default_for_version(getattr(ModelVersion, TABPFN_VERSIONS[name]), device="cpu", random_state=SEED)
    if name == "tabicl_v2":
        import tabicl

        return (tabicl.TabICLClassifier if c else tabicl.TabICLRegressor)(random_state=SEED, device="cpu")
    if name == "tabdpt":
        from tabdpt import TabDPTClassifier, TabDPTRegressor

        return (TabDPTClassifier if c else TabDPTRegressor)(device="cpu")
    raise KeyError(name)


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout(f"fit exceeded {FIT_LIMIT_S // 60} minutes")


def fit_predict(name: str, task: str, Xtr, ytr, Xte):
    """Fit one model on a prepared training fold and predict the held out fold. Returns (pred, info)."""
    from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold

    info = {}
    if name in CLASSICAL:
        est, grid = classical(name, task)
        cv = StratifiedKFold(INNER, shuffle=True, random_state=SEED) if task == "classification" else KFold(INNER, shuffle=True, random_state=SEED)
        gs = GridSearchCV(est, grid, cv=cv, scoring="roc_auc" if task == "classification" else "neg_mean_absolute_error",
                          n_jobs=1, error_score=np.nan, refit=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gs.fit(Xtr, ytr)
        failed = int(np.isnan(gs.cv_results_["mean_test_score"]).sum())
        info = {"best_params": json.dumps(gs.best_params_, default=str), "grid_points_failed": failed}
        m = gs.best_estimator_
    elif name == "taco_ensemble":
        m = taco(task)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(Xtr, ytr)
    elif name == "mitra":
        return mitra(task, Xtr, ytr, Xte), info
    else:
        m = foundation(name, task)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m.fit(Xtr, ytr)
        if name in TABPFN_VERSIONS and name not in ("tabpfn_v2", "tabpfn_2_5"):
            info = {"checkpoint": Path(str(m.model_path)).name}
    pred = m.predict_proba(Xte)[:, 1] if task == "classification" else m.predict(Xte)
    return np.asarray(pred, dtype=float), info


def mitra(task: str, Xtr, ytr, Xte):
    py = CODE / ".venv-mitra/bin/python"
    if not py.exists():
        raise ImportError("Mitra environment Code/.venv-mitra not installed")
    with tempfile.TemporaryDirectory() as tmp:
        np.savez(f"{tmp}/in.npz", Xtr=Xtr, ytr=ytr, Xte=Xte, task=task)
        r = subprocess.run([str(py), "-I", str(CTRL / "mitra_worker.py"), f"{tmp}/in.npz", f"{tmp}/out.npz"],
                           capture_output=True, text=True, timeout=FIT_LIMIT_S)
        if r.returncode != 0:
            raise RuntimeError(r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "mitra worker failed")
        return np.load(f"{tmp}/out.npz")["pred"]


def classify_failure(e: BaseException) -> str:
    s = f"{type(e).__name__} {e}".lower()
    if isinstance(e, (Timeout, subprocess.TimeoutExpired)):
        return "time"
    if isinstance(e, MemoryError) or "out of memory" in s:
        return "memory"
    if "license" in s or "licence" in s or "token" in s or "gated" in s:
        return "licence or token"
    if isinstance(e, (ImportError, ModuleNotFoundError)):
        return "install"
    if "cuda" in s or "gpu" in s:
        return "GPU"
    return "error"


# --------------------------------------------------------------------------------------------- metrics

def ece(y, p, bins=10):
    idx = np.minimum((p * bins).astype(int), bins - 1)
    return float(sum(abs(y[idx == b].mean() - p[idx == b].mean()) * (idx == b).mean() for b in range(bins) if (idx == b).any()))


def metrics(y, pred, task: str) -> dict:
    from sklearn import metrics as M

    if task == "classification":
        p = np.clip(pred, 1e-15, 1 - 1e-15)
        yh = (p >= 0.5).astype(int)
        tn, fp, fn, tp = M.confusion_matrix(y, yh, labels=[0, 1]).ravel()
        return {"AUROC": M.roc_auc_score(y, p), "AUPRC": M.average_precision_score(y, p),
                "balanced_accuracy": M.balanced_accuracy_score(y, yh), "F1": M.f1_score(y, yh, zero_division=0),
                "MCC": M.matthews_corrcoef(y, yh), "sensitivity": tp / (tp + fn) if tp + fn else np.nan,
                "specificity": tn / (tn + fp) if tn + fp else np.nan, "Brier": M.brier_score_loss(y, p),
                "log_loss": M.log_loss(y, p, labels=[0, 1]), "ECE": ece(y, p)}
    return {"MAE": M.mean_absolute_error(y, pred), "RMSE": float(np.sqrt(M.mean_squared_error(y, pred))),
            "R2": M.r2_score(y, pred), "Spearman": stats.spearmanr(y, pred).statistic}


def corrected_t(d: np.ndarray, n_test: float, n_train: float) -> tuple:
    """Nadeau and Bengio corrected resampled t: mean, 95% interval (t with J - 1 degrees of freedom)."""
    d = d[np.isfinite(d)]
    J = len(d)
    if J < 2:
        return np.nan, np.nan, np.nan
    var = (1 / J + n_test / n_train) * d.var(ddof=1)
    h = stats.t.ppf(0.975, J - 1) * np.sqrt(var)
    return d.mean(), d.mean() - h, d.mean() + h


def verdict(lo: float, hi: float) -> str:
    if not np.isfinite(lo):
        return "not evaluable"
    return "improves" if lo > 0 else ("gets worse" if hi < 0 else "does not change")


# ---------------------------------------------------------------------------------------------- runner

def _job(args):
    """One (job, fold) in a worker process: per fit time limit with SIGALRM."""
    job, fold, cache = args
    d = np.load(cache, allow_pickle=True)
    task, name = job["task"], job["model"]
    tr, te = d[f"tr{fold}"], d[f"te{fold}"]
    t0 = time.time()
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(FIT_LIMIT_S)
    try:
        pred, info = run_fold(job, d, tr, te)
        signal.alarm(0)
        y = d["y"][te]
        return dict(ok=True, fold=fold, seconds=time.time() - t0, pred=pred, info=info, metrics=metrics(y, pred, task))
    except BaseException as e:  # noqa: BLE001 (every failure is recorded)
        signal.alarm(0)
        return dict(ok=False, fold=fold, seconds=time.time() - t0, reason=classify_failure(e),
                    message=f"{type(e).__name__}: {str(e)[:300]}", trace=traceback.format_exc()[-800:])


def run_fold(job, d, tr, te):
    """Arms of steps 01 and 02 (one matrix, optional extra training rows) and the single layer arms of step 03."""
    if job.get("module"):  # another folder's arms (Round2/5_Clinical): its module is on sys.path
        import importlib

        return importlib.import_module(job["module"]).run_fold(job, d, tr, te)
    task, name, k = job["task"], job["model"], job["k"]
    X, y = d["X"], d["y"]
    Xtr, ytr, Xte = X[tr], y[tr], X[te]
    if "Xp" in d.files:  # partner rows added to every training fold (B_naive, B_scaled)
        Xtr, ytr = np.vstack([Xtr, d["Xp"]]), np.concatenate([ytr, d["yp"]])
    if job.get("scaled"):  # B_scaled: anchor rows scaled with the anchor training fold, partner with its own statistics
        Xa = X[tr]
        Xtr = np.vstack([zscore_within(Xa), zscore_within(d["Xp"])])
        Xte = zscore_within(X[te], ref=Xa)
    if name in EXPRESSION:  # precomputed embeddings: no feature selection; standardised, then tuned linear model
        E, Ep = d[f"E_{name}"], d[f"Ep_{name}"] if f"Ep_{name}" in d.files else None
        Etr, Ete = E[tr], E[te]
        if Ep is not None:
            Etr = np.vstack([zscore_within(E[tr]), zscore_within(Ep)]) if job.get("scaled") else np.vstack([Etr, Ep])
            if job.get("scaled"):
                Ete = zscore_within(E[te], ref=E[tr])
        mu, sd = Etr.mean(0), Etr.std(0)
        sd = np.where(sd > 0, sd, 1.0)
        return fit_predict("logistic_regression", task, (Etr - mu) / sd, ytr, (Ete - mu) / sd)
    if job.get("vertical"):  # step 03: same participants, two layers (X and X2), 100 features per layer
        arm = job["arm"]
        if arm == "layer_2":
            Xtr, Xte = d["X2"][tr], d["X2"][te]
        if arm == "joined":
            A1, B1, _ = prepare(X[tr], ytr, X[te], task, k)
            A2, B2, _ = prepare(d["X2"][tr], ytr, d["X2"][te], task, k)
            return fit_predict(name, task, np.hstack([A1, A2]), ytr, np.hstack([B1, B2]))
    A, B, _ = prepare(Xtr, ytr, Xte, task, k)
    return fit_predict(name, task, A, ytr, B)


def bulkformer_embed(counts: pd.DataFrame) -> np.ndarray:
    """BulkFormer-37M (the smallest released model), sample level mean embedding from raw counts (log1p TPM inside)."""
    import torch
    from bulkformer import BulkFormerModel, extract_features
    from bulkformer.preprocessing import align_to_vocab, normalize_data

    torch.manual_seed(SEED)
    m = BulkFormerModel.from_pretrained("zhfanrui/BulkFormer-37M")
    al, mp = align_to_vocab(normalize_data(counts), m.config.gene_ids)
    dev = "mps" if torch.backends.mps.is_available() else "cpu"
    e = extract_features(m, al.values, output_feature_type="sample_level", aggregate_type="mean", mask_prob=mp, device=dev, batch_size=16)
    return np.asarray(e.detach().cpu() if hasattr(e, "detach") else e, dtype=float)


def save_cache(path: Path, X, y, fold_list, **extra):
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"X": np.asarray(X, dtype=float), "y": np.asarray(y, dtype=float)}
    for i, (tr, te) in enumerate(fold_list):
        arrays[f"tr{i}"], arrays[f"te{i}"] = tr, te
    arrays.update({k: np.asarray(v, dtype=float) for k, v in extra.items() if v is not None})
    np.savez(path, **arrays)


class Runner:
    """Resumable: every finished (job, fold) is appended to raw/<job>.csv (metrics) and raw/<job>_pred.csv
    (held out predictions); a rerun skips them. Failures go to failures.csv (task, arm, model, reason, message,
    device, memory) and the job continues with the next fold."""

    def __init__(self, out: Path, step: str, folder: str = "0_Controls"):
        self.out, self.step, self.folder = out, step, folder
        (out / "raw").mkdir(parents=True, exist_ok=True)
        self.failures = out / "failures.csv"

    def done(self, jid: str) -> set:
        f = self.out / "raw" / f"{jid}.csv"
        return set(pd.read_csv(f)["fold"]) if f.exists() else set()

    def retry_licence(self, jobs: list):
        """Amendment 3, item 3: TabPFN 2.5 jobs that failed on the licence are tried again. Their old rows leave
        failures.csv; a new failure writes a new row, a success leaves none."""
        if not self.failures.exists():
            return
        F = pd.read_csv(self.failures)
        ids = {j["id"] for j in jobs if j["model"] == "tabpfn_2_5"}
        drop = F["job"].isin(ids) & (F["reason"] == "licence or token")
        if drop.any():
            F[~drop].to_csv(self.failures, index=False)
            print(f"retrying {int(drop.sum())} TabPFN 2.5 jobs that failed on the licence", flush=True)

    def failed(self, jid: str) -> bool:
        if not self.failures.exists():
            return False
        F = pd.read_csv(self.failures)
        return bool((F["job"] == jid).any())

    def record(self, job, r, n_train, n_test):
        jid = job["id"]
        if r["ok"]:
            row = dict(job=jid, step=self.step, control=job["control"], layer=job["layer"], arm=job["arm"], model=job["model"],
                       fold=r["fold"], seconds=round(r["seconds"], 2), n_train=n_train, n_test=n_test, **r["metrics"], **r["info"])
            f = self.out / "raw" / f"{jid}.csv"
            pd.DataFrame([row]).to_csv(f, mode="a", header=not f.exists(), index=False)
            p = self.out / "raw" / f"{jid}_pred.csv"
            pd.DataFrame({"fold": r["fold"], "pred": r["pred"]}).to_csv(p, mode="a", header=not p.exists(), index=False)
        else:
            import psutil

            row = dict(task=f"{self.folder}/{self.step}/{job['control']}/{job['layer']}", arm=job["arm"], model=job["model"],
                       reason=r["reason"], message=r["message"], device="cpu (Apple M4 Pro, MPS not used)",
                       memory=f"{psutil.virtual_memory().total / 2**30:.0f} GB", job=jid, fold=r["fold"])
            pd.DataFrame([row]).to_csv(self.failures, mode="a", header=not self.failures.exists(), index=False)

    def run(self, jobs: list, fold_ids: list, workers: int = WORKERS):
        """Run the given folds of every job; a job that failed is not retried on later folds."""
        todo = []
        for job in jobs:
            if self.failed(job["id"]):
                continue
            have = self.done(job["id"])
            todo += [(job, f) for f in fold_ids if f not in have]
        heavy = [t for t in todo if t[0]["model"] in FOUNDATION + EXPRESSION + ["taco_ensemble"]]
        light = [t for t in todo if t not in heavy]
        stopped = set()
        for batch, w in ((light, workers), (heavy, 2)):
            if not batch:
                continue
            with ProcessPoolExecutor(max_workers=w) as ex:
                futs = {ex.submit(_job, (j, f, j["cache"])): (j, f) for j, f in batch}
                for fu in as_completed(futs):
                    j, f = futs[fu]
                    r = fu.result()
                    if j["id"] in stopped:
                        continue
                    d = np.load(j["cache"])
                    n_te = len(d[f"te{f}"])
                    n_tr = len(d[f"tr{f}"]) + (len(d["yp"]) if "yp" in d.files else 0)
                    self.record(j, r, n_tr, n_te)
                    if not r["ok"]:
                        stopped.add(j["id"])
                        for g, (jj, ff) in futs.items():
                            if jj["id"] == j["id"]:
                                g.cancel()
                    print(f"{j['id']} fold {f}: {'ok' if r['ok'] else 'FAILED ' + r['reason']} {r['seconds']:.1f}s", flush=True)


def mitra_measured(control: str, layer: str) -> tuple:
    """(seconds per fit, source) for Mitra: the step 01 fold 0 fit of the same control and layer (anchor, all
    features); for a layer step 01 did not time (step 03 cohort pairs), the largest step 01 fit of the control."""
    T = [pd.read_csv(f) for f in sorted(MITRA_TIMED.glob("*__mitra.csv"))]
    if not T:
        return np.nan, "not measured"
    T = pd.concat(T, ignore_index=True)
    hit = T[(T["control"] == control) & (T["layer"] == layer)]
    if len(hit):
        return float(hit["seconds"].iloc[0]), f"measured: step 01, {control} {layer}, A_all fold 0 (n_train {int(hit['n_train'].iloc[0])})"
    c = T[T["control"] == control]
    what = f"largest step 01 {control} fit"
    if c.empty:  # another folder's task (Round2/5_Clinical): the largest step 01 fit of any control
        c, what = T, "largest 0_Controls step 01 fit"
    return float(c["seconds"].max()), f"not timed for this pair; {what} used (n_train {int(c.loc[c['seconds'].idxmax(), 'n_train'])})"


def estimate(out: Path, jobs: list) -> pd.DataFrame:
    """Estimated wall clock of the full run from the timed fold 0 of every job (light models run WORKERS at a
    time, heavy ones two at a time). Removed models (Amendment 2) keep their rows, marked not run, and are left out
    of the step total."""
    rows = []
    for j in jobs:
        f = out / "raw" / f"{j['id']}.csv"
        if j["model"] in REMOVED:
            t, src = mitra_measured(j["control"], j["layer"])
            status = f"not run: {REMOVED[j['model']]} (Amendment 2); {src}"
        elif f.exists() and (pd.read_csv(f)["fold"] == 0).any():
            t, status = float(pd.read_csv(f).query("fold == 0")["seconds"].iloc[0]), "run"
        else:
            continue
        par = 2 if j["model"] in FOUNDATION + EXPRESSION + ["taco_ensemble"] else WORKERS
        rows.append(dict(job=j["id"], model=j["model"], fold0_seconds=t, cpu_hours=t * N_FOLDS / 3600,
                         wall_hours=t * N_FOLDS / 3600 / par, status=status))
    return pd.DataFrame(rows)


def run_step(runner: "Runner", jobs: list, out: Path, step: str, time_only: bool) -> bool:
    """Amendment 1, item 8 and Amendment 2, item 2: fold 0 of every job, the estimate (removed models left out), then
    the full run (N_FOLDS folds, Amendment 3) only when the step fits in BUDGET_H. Returns True when the full run was done."""
    runner.retry_licence(jobs)
    runner.run(jobs, [0])
    E = estimate(out, jobs)
    E.to_csv(out / "timing.csv", index=False)
    run = E[E["status"] == "run"]
    total = run["wall_hours"].sum()
    print(f"estimated wall clock of step {step} (models run): {total:.1f} h; budget {BUDGET_H:.0f} h per step")
    print(f"not counted (Amendment 2): {E.loc[E['status'] != 'run', 'wall_hours'].sum():.1f} h of removed models")
    if time_only:
        return False
    if total > BUDGET_H:
        print("over the budget: stopping before the full run. Estimate per model (h):")
        print(run.groupby("model")["wall_hours"].sum().sort_values(ascending=False).round(2).to_string())
        return False
    runner.run(jobs, list(range(N_FOLDS)))
    return True


# ------------------------------------------------------------------------------------- caches and jobs

ARMS_H = ["A_all", "A_shared", "B_naive", "B_scaled"]
ARMS_V = ["layer_1", "layer_2", "joined"]  # late fusion is the mean of the layer_1 and layer_2 predictions (step 03)


def embeddings_for(study: str, ids: list, genes=None) -> dict:
    """Expression model embeddings of the given samples (raw counts; genes restricted when given). Unsupervised and per
    sample, so computed once per cohort and gene set, outside the folds. A failure is returned as an exception object."""
    X, m, is_counts, _ = load(study)
    out = {}
    counts = X.loc[ids] if genes is None else X.loc[ids, [g for g in X.columns if g in set(genes)]]
    for name in EXPRESSION:
        t0 = time.time()
        try:
            if not is_counts:
                raise ValueError("needs raw counts; the matrix is not counts")
            if name == "bulkformer":
                out[name] = bulkformer_embed(counts)
            else:
                out[name] = bulkrnabert_embed(counts)
            if time.time() - t0 > FIT_LIMIT_S:
                raise Timeout(f"embedding took {time.time() - t0:.0f} s for {len(ids)} samples (limit {FIT_LIMIT_S} s)")
        except BaseException as e:  # noqa: BLE001
            out[name] = e
    return out


BULKRNABERT_SECONDS_PER_SAMPLE = 10.2  # measured 20 November 2026: one TCGA sample on the CPU; on MPS 2 samples ran out of memory


def bulkrnabert_embed(counts: pd.DataFrame) -> np.ndarray:
    """InstaDeepAI/BulkRNABert (transformers 4.51 in Code/.venv-bulkrnabert, the version of its config; it does not load
    under transformers 5). Full attention over 19,062 genes: measured 10.2 s per sample on the CPU and out of memory on
    MPS, so embedding a training fold exceeds the 15 minute limit per fit; recorded as a time failure without running."""
    n = len(counts)
    raise Timeout(f"measured {BULKRNABERT_SECONDS_PER_SAMPLE} s per sample on the CPU (MPS out of memory at 2 samples); "
                  f"{n} samples need about {n * BULKRNABERT_SECONDS_PER_SAMPLE / 60:.0f} min, over the {FIT_LIMIT_S // 60} minute limit per fit")


def build_horizontal(control: str, layer: str, anchor: str, partner: str) -> dict:
    """Caches of the four horizontal arms; returns {arm: cache path} and writes the participants used."""
    task = TASK[control]
    Xa, ya, ma = control_data(anchor, control)
    Xp, yp, mp = control_data(partner, control) if partner else (None, None, None)
    fl = folds(ya, task)
    paths, emb_fail = {}, {}
    shared = [c for c in Xa.columns if partner and c in set(Xp.columns)]
    E_all = embeddings_for(anchor, list(Xa.index)) if layer == "transcriptomics" else {}
    E_sh = embeddings_for(anchor, list(Xa.index), shared) if layer == "transcriptomics" and partner else {}
    Ep_sh = embeddings_for(partner, list(Xp.index), shared) if layer == "transcriptomics" and partner else {}
    for arm in ARMS_H:
        if arm != "A_all" and not partner:
            continue
        p = CACHE / f"{control}_{layer}_{arm}.npz"
        extra = {}
        if arm == "A_all":
            X = Xa
            for n, e in E_all.items():
                extra[f"E_{n}"] = None if isinstance(e, BaseException) else e
                if isinstance(e, BaseException):
                    emb_fail[(arm, n)] = e
        else:
            X = Xa[shared]
            for n in E_sh:
                bad = E_sh[n] if isinstance(E_sh[n], BaseException) else (Ep_sh[n] if isinstance(Ep_sh[n], BaseException) else None)
                if bad is not None:
                    emb_fail[(arm, n)] = bad
                    continue
                extra[f"E_{n}"] = E_sh[n]
                if arm != "A_shared":
                    extra[f"Ep_{n}"] = Ep_sh[n]
            if arm != "A_shared":
                extra.update(Xp=Xp[shared].to_numpy(), yp=yp)
        save_cache(p, X.to_numpy(), ya, fl, **extra)
        paths[arm] = p
    info = dict(control=control, layer=layer, anchor=anchor, partner=partner, n_anchor=len(ya), n_partner=0 if yp is None else len(yp),
                n_features_anchor=Xa.shape[1], n_shared=len(shared))
    return paths, emb_fail, info, (ma, mp)


def jobs_for(control: str, layer: str, paths: dict, arms: list, step: str) -> list:
    task = TASK[control]
    out = []
    for arm in arms:
        if arm not in paths:
            continue
        for model in models_for(layer):
            out.append(dict(id=f"{step}__{control}__{layer}__{arm}__{model}", control=control, layer=layer, arm=arm, model=model,
                            task=task, k=K_FEATURES, scaled=arm == "B_scaled", cache=str(paths[arm])))
    return out


def record_static_failures(runner: Runner, jobs: list, emb_fail: dict):
    """Failures known before any fold: removed models (Amendment 2), no regression form (Gaussian naive Bayes),
    embedding failures."""
    for j in jobs:
        if runner.failed(j["id"]) or runner.done(j["id"]):
            continue
        e = None
        if j["model"] in REMOVED:
            t, src = mitra_measured(j["control"], j["layer"])
            runner.record(j, dict(ok=False, fold=-1, seconds=0.0, reason=REMOVED[j["model"]],
                                  message=f"Amendment 2: removed from the model grid; {t:.0f} s per fit ({src}); estimated total "
                                          f"{t * N_FOLDS / 3600:.1f} h CPU for the {N_FOLDS} folds of this job "
                                          f"({t * N_FOLDS / 3600 / 2:.1f} h wall clock at two at a time)"), 0, 0)
            continue
        if j["model"] == "gaussian_nb" and j["task"] == "regression":
            e = NotImplementedError("Gaussian naive Bayes has no regression form")
        elif j["model"] in EXPRESSION and (j["arm"], j["model"]) in emb_fail:
            e = emb_fail[(j["arm"], j["model"])]
        if e is not None:
            runner.record(j, dict(ok=False, fold=-1, seconds=0.0, reason=classify_failure(e) if not isinstance(e, NotImplementedError) else "error",
                                  message=f"{type(e).__name__}: {str(e)[:300]}"), 0, 0)
