"""Build release v1 of the harmonised datasets into release/ (README.md is written by hand; this script writes the data,
the split files and manifest.csv). Existing results are only read.

  early_detection/  GSE215865 development set (earliest library per subject, timepoint 0 marked, later LC labels);
                    Su2022 and Cervia2024_SomaScan acute samples on their 386 shared proteins, without harmonisation and
                    with center_scale fitted per cohort without labels (frozen test sets); Su2022 acute Olink and
                    metabolomics on the same participants; split files; baselines from results/early_detection
  clustering/       LC samples per layer, without harmonisation and with center_scale and fsqn fitted on the LC samples of
                    all cohorts of the layer without labels; sample metadata (cohort, sex, age, days since infection)
  signature/        per layer, per cohort Hedges g of LC versus COVID_recovered on each cohort's own features, the random
                    effects pooled estimate (DerSimonian and Laird) with 95% interval, Q, I2, tau2, the number of cohorts
                    and the number agreeing in sign, and direction_inconsistent

Matrices are CSV (gzip) with samples in rows: row_id (the benchmark's pooled id, "<cohort>::<sample>"; the sample id for
the paired Su2022 files), then the features; values with 7 significant digits. tests/test_release.py checks that every
matrix reproduces the values the benchmark loads for the same cohort, samples and features.

  python scripts/build_release.py [--out release]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from concordance import feature_symbols  # noqa: E402
from paired_layers import load_layers, paired_samples, PAIRED  # noqa: E402

from chorus.cli import load_config  # noqa: E402
from chorus.eval.paired import PairedCV  # noqa: E402
from chorus.features.invariant import meta_analyse, per_study_effects  # noqa: E402
from chorus.harmonise import make_harmoniser  # noqa: E402
from chorus.schema import binary_label  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
HGNC = str(ROOT / "data/ref/hgnc_complete_set.txt")
FMT = "%.7g"
ROLE = {"PRJNA1149196_QLD": "training_only", "CORSIP2026": "test_only", "GSE228320": "test_only", "Saito2024": "concordance_only"}
FLAGS = {"GSE270045": "label_confounded_with_collection", "Cervia2024_MS": "same_participants_as_Cervia2024_SomaScan"}
NOT_POOLED = {"Cervia2024_MS"}  # same sera as Cervia2024_SomaScan: reported, not pooled twice
CLUSTER_LAYERS = {
    "transcriptomics_v4": "configs/lc_transcriptomics_blood_v4.yaml",
    "proteomics_affinity_v5": "configs/lc_proteomics_affinity_v5.yaml",
    "proteomics_explore_v5": "configs/lc_proteomics_explore_v5.yaml",
    "proteomics_ms_v4": "configs/lc_proteomics_ms_v4.yaml",
    "metabolomics": "configs/lc_metabolomics.yaml",
}
SIGNATURE_LAYERS = {
    "transcriptomics": "configs/concordance_rna_v5.yaml",
    "proteomics": "configs/concordance_proteomics_all_v5.yaml",
    "metabolomics": "configs/concordance_metabolomics_all.yaml",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Release:
    def __init__(self, out: Path, commit: str):
        self.out, self.commit, self.rows = out, commit, []

    def write(self, rel: str, df: pd.DataFrame, *, kind: str, layer: str, cohorts, harmoniser: str = "", config: str = "",
              n_features=None, y=None, roles: dict = None, index_label: str = "row_id", source: dict = None) -> None:
        path = self.out / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, float_format=FMT, index=index_label is not None, index_label=index_label)
        cohorts = list(dict.fromkeys(cohorts))
        roles = roles or {}
        self.rows.append({
            "path": rel, "sha256": sha256(path), "kind": kind, "layer": layer, "cohorts": ";".join(cohorts),
            "roles": ";".join(f"{c}:{roles.get(c, ROLE.get(c, 'training_and_held_out'))}" for c in cohorts),
            "n_rows": len(df), "n_features": n_features if n_features is not None else "", "harmoniser": harmoniser,
            "n_LC": int((y == 1).sum()) if y is not None else "", "n_COVID_recovered": int((y == 0).sum()) if y is not None else "",
            "config": config, "commit": self.commit,
            "flags": ";".join(f"{c}:{FLAGS[c]}" for c in cohorts if c in FLAGS),
            "source": json.dumps(source or {}, sort_keys=True),
        })
        print(f"{rel}: {df.shape}", flush=True)


def pooled_labelled(config: str, features=None):
    cs, lc, cfg = load_config(str(ROOT / config))
    X, M = cs.pooled(features if features is not None else cfg["filtered_features"])
    y = binary_label(M, lc.positive, lc.negative)
    keep = y.notna()
    return X[keep], M[keep], y[keep].astype(int), cs, lc, cfg


def harmonise(name: str, X: pd.DataFrame, s: pd.Series) -> pd.DataFrame:
    H = make_harmoniser(name).fit(X, s)  # label free: no y passed
    return H.transform(X, s)


def sample_table(M: pd.DataFrame, y: pd.Series = None, extra=()) -> pd.DataFrame:
    cols = ["study_id", "sample_id", "subject_id", "lc_status", "sex", "age", "days_since_infection", "time_anchor", "sample_matrix", "sample_phase", "timepoint", *extra]
    t = M.reindex(columns=[c for c in cols if c in M.columns or c in ("sex", "age", "days_since_infection")]).rename(columns={"study_id": "cohort"})
    if y is not None:
        t.insert(3, "y_LC", y.loc[t.index].values)
    return t


# ------------------------------------------------------------------------------------------------ early detection
def early_detection(R: Release) -> None:
    d = "early_detection"
    # development set
    cfg_dev, cfg_ref = "configs/release_v1_gse215865_dev.yaml", "configs/early_detection_reference.yaml"
    X, M, y, *_ = pooled_labelled(cfg_dev)
    Xr, Mr, yr, _, lcr, _ = pooled_labelled(cfg_ref)
    t0 = M.index.isin(Xr.index)
    S = sample_table(M, y)
    S.insert(4, "is_timepoint0", t0)
    R.write(f"{d}/gse215865_development_log2cpm.csv.gz", X, kind="matrix", layer="transcriptomics", cohorts=["GSE215865"], harmoniser="none",
            config=cfg_dev, n_features=X.shape[1], y=y, roles={"GSE215865": "training"}, source={"loader": "config", "config": cfg_dev})
    R.write(f"{d}/gse215865_development_samples.csv", S, kind="samples", layer="transcriptomics", cohorts=["GSE215865"], config=cfg_dev, y=y, roles={"GSE215865": "training"})
    # the folds of run_within_study_cv on the timepoint 0 samples (repeat r: StratifiedKFold(5, shuffle, random_state r))
    sp = pd.DataFrame(index=M.index)
    for r in range(10):
        f = pd.Series(np.nan, index=Xr.index)
        for k, (_, te) in enumerate(StratifiedKFold(n_splits=5, shuffle=True, random_state=lcr.random_state + r).split(Xr, yr)):
            f.iloc[te] = k
        sp[f"fold_repeat{r}"] = f.reindex(M.index).astype("Int64")
    sp.insert(0, "is_timepoint0", t0)
    R.write(f"{d}/gse215865_development_splits.csv", sp, kind="splits", layer="transcriptomics", cohorts=["GSE215865"], config=cfg_ref, y=y, roles={"GSE215865": "training"})

    # frozen acute protein test sets
    cfg_h = "configs/early_detection_horizontal.yaml"
    X, M, y, cs, lc, cfg = pooled_labelled(cfg_h)
    s = M["study_id"].astype(str)
    roles = {c: "test_only" for c in s.unique()}
    for h in ("none", "center_scale"):
        Z = X if h == "none" else harmonise(h, X, s)
        R.write(f"{d}/acute_proteins_386_{h}.csv.gz", Z, kind="matrix", layer="proteomics", cohorts=list(s.unique()), harmoniser=h, config=cfg_h,
                n_features=Z.shape[1], y=y, roles=roles, source={"loader": "config", "config": cfg_h, "harmoniser": h, "groups": "cohort"})
    R.write(f"{d}/acute_proteins_386_samples.csv", sample_table(M, y), kind="samples", layer="proteomics", cohorts=list(s.unique()), config=cfg_h, y=y, roles=roles)
    sp = pd.DataFrame({"cohort": s, "train_Su2022_test_Cervia2024_SomaScan": np.where(s == "Su2022", "train", "test"),
                       "train_Cervia2024_SomaScan_test_Su2022": np.where(s == "Su2022", "test", "train")}, index=M.index)
    R.write(f"{d}/acute_proteins_386_splits.csv", sp, kind="splits", layer="proteomics", cohorts=list(s.unique()), config=cfg_h, y=y, roles=roles)

    # Su2022 acute Olink and metabolomics on the same participants
    views, metas, impute, _ = load_layers("Su2022")
    samples = paired_samples(views, metas, "early", PAIRED["Su2022"]["timepoint_order"])
    yv = samples["y"].astype(int)
    cv = PairedCV({n: V.loc[yv.index] for n, V in views.items()}, yv, impute=impute)
    roles = {"Su2022": "training"}
    for n, V in cv.views.items():
        R.write(f"{d}/su2022_paired_{n}.csv.gz", V, kind="matrix", layer=n, cohorts=["Su2022"], harmoniser="none", config="scripts/paired_layers.py (Su2022, contrast early)",
                n_features=V.shape[1], y=yv, roles=roles, source={"loader": "paired", "cohort": "Su2022", "layer": n})
    st = samples.assign(y_LC=yv, lc_status=np.where(yv == 1, "LC", "COVID_recovered")).rename_axis("sample_id").reset_index()
    st = st.merge(metas["olink"][["sex", "age", "days_since_infection"]], left_on="sample_id", right_index=True, how="left").set_index("sample_id")
    R.write(f"{d}/su2022_paired_samples.csv", st, kind="samples", layer="olink;metabolomics", cohorts=["Su2022"], config="scripts/paired_layers.py", y=yv, roles=roles, index_label="sample_id")
    F = cv.folds().pivot(index="sample", columns="repeat", values="fold").rename(columns=lambda r: f"fold_repeat{r}").loc[yv.index]
    R.write(f"{d}/su2022_paired_splits.csv", F, kind="splits", layer="olink;metabolomics", cohorts=["Su2022"], config="scripts/paired_layers.py", y=yv, roles=roles, index_label="sample_id")

    # baselines
    E = ROOT / "results/early_detection"
    h = pd.read_csv(E / "horizontal/loso_results.csv")
    h = h[h["auroc"].notna()].merge(pd.read_csv(E / "horizontal/permutation_summary.csv"), on=["held_out", "harmoniser", "selector", "classifier"])
    base = [pd.DataFrame({"set": "acute_proteins_386", "test": h["held_out"], "trained_on": np.where(h["held_out"] == "Su2022", "Cervia2024_SomaScan", "Su2022"),
                          "model": h["harmoniser"] + " / " + h["selector"] + " / " + h["classifier"], "auroc": h["auroc"], "auroc_lo": h["auroc_lo"], "auroc_hi": h["auroc_hi"],
                          "p_permutation": h["p_permutation"], "evaluation": "train on one cohort, test on the other"})]
    v = pd.read_csv(E / "vertical/early/summary.csv")
    base.append(pd.DataFrame({"set": "su2022_paired", "test": "Su2022", "trained_on": "Su2022", "model": v["model"], "auroc": v["auroc"], "auroc_sd": v["auroc_sd"],
                              "evaluation": "repeated stratified 5 fold CV, 10 repeats, mean over 50 folds (su2022_paired_splits.csv)"}))
    w = pd.read_csv(E / "reference/within_study_cv.csv")
    base.append(pd.DataFrame({"set": "gse215865_development", "test": "GSE215865 timepoint 0", "trained_on": "GSE215865", "model": w["selector"] + " / " + w["classifier"],
                              "auroc": w["auroc"], "auroc_sd": w["auroc_sd"], "evaluation": "repeated stratified 5 fold CV on the timepoint 0 samples, mean over 50 folds (gse215865_development_splits.csv)"}))
    base.append(pd.DataFrame([{"set": "gse215865_development", "test": "GSE215865", "trained_on": "GSE215865", "model": "TACO (published, Thompson 2023)", "auroc": 0.675, "evaluation": "published value, not recomputed"}]))
    B = pd.concat(base, ignore_index=True)
    R.write(f"{d}/baselines.csv", B, kind="baselines", layer="", cohorts=["GSE215865", "Su2022", "Cervia2024_SomaScan"], index_label=None,
            roles={"GSE215865": "training", "Su2022": "test_only", "Cervia2024_SomaScan": "test_only"})


# ------------------------------------------------------------------------------------------------ clustering
def clustering(R: Release) -> None:
    for layer, config in CLUSTER_LAYERS.items():
        X, M, y, cs, lc, cfg = pooled_labelled(config)
        lcm = (y == 1).to_numpy()
        X, M, y = X[lcm], M[lcm], y[lcm]
        s = M["study_id"].astype(str)
        coh = list(dict.fromkeys(s))
        for h in ("none", "center_scale", "fsqn"):
            Z = X if h == "none" else harmonise(h, X, s)
            R.write(f"clustering/{layer}_{h}.csv.gz", Z, kind="matrix", layer=layer, cohorts=coh, harmoniser=h, config=config, n_features=Z.shape[1], y=y,
                    source={"loader": "config", "config": config, "samples": "LC", "harmoniser": h, "groups": "cohort"})
        R.write(f"clustering/{layer}_samples.csv", sample_table(M).drop(columns=["lc_status"], errors="ignore"), kind="samples", layer=layer, cohorts=coh, config=config, y=y)


# ------------------------------------------------------------------------------------------------ signature
def signature(R: Release) -> None:
    for layer, config in SIGNATURE_LAYERS.items():
        cs, lc, cfg = load_config(str(ROOT / config))
        X, M = cs.pooled(cfg["filtered_features"])
        y = binary_label(M, lc.positive, lc.negative)
        s = M["study_id"].astype(str)
        keep = y.notna() & ~s.map(lambda c: ROLE.get(c) == "test_only")
        X, M, y, s = X[keep], M[keep], y[keep].astype(int), s[keep]
        eff = per_study_effects(X, y, s)
        G = pd.DataFrame({c: e["g"] for c, e in eff.items()})
        V = pd.DataFrame({c: e["var"] for c, e in eff.items()})
        cohorts = list(G.columns)

        def pool(use, prefix):
            g, v = G[use], V[use]
            ok = g.notna().sum(axis=1) >= 2
            m = meta_analyse({c: pd.DataFrame({"g": g.loc[ok, c], "var": v.loc[ok, c]}) for c in use})
            agree = (np.sign(g.loc[ok]).eq(np.sign(m["effect"]), axis=0) & g.loc[ok].notna()).sum(axis=1)
            out = pd.DataFrame({f"{prefix}g": m["effect"], f"{prefix}se": m["se"], f"{prefix}ci_lo": m["effect"] - 1.959964 * m["se"], f"{prefix}ci_hi": m["effect"] + 1.959964 * m["se"],
                                f"{prefix}p": m["p"], f"{prefix}Q": m["Q"], f"{prefix}Q_p": m["Q_p"], f"{prefix}I2": m["I2"], f"{prefix}tau2": m["tau2"],
                                f"{prefix}n_cohorts": m["k"].astype(int), f"{prefix}n_agree_sign": agree.astype(int)}).reindex(G.index)
            out[f"{prefix}direction_inconsistent"] = (out[f"{prefix}n_agree_sign"] < out[f"{prefix}n_cohorts"]).where(out[f"{prefix}n_cohorts"].notna())
            return out

        pooled_cohorts = [c for c in cohorts if c not in NOT_POOLED]
        T = pd.concat([G.add_prefix("g_"), pool(pooled_cohorts, "pooled_")], axis=1)
        flagged = [c for c in pooled_cohorts if c in FLAGS]
        if flagged:
            T = pd.concat([T, pool([c for c in pooled_cohorts if c not in flagged], "pooled_unflagged_")], axis=1)
        sym = feature_symbols(T.index, HGNC, cs.cohorts[0].feature_space)
        T.insert(0, "symbol", sym.reindex(T.index).values)
        T = T.rename_axis("feature")
        R.write(f"signature/{layer}_signature.csv", T, kind="signature", layer=layer, cohorts=cohorts, config=config, n_features=len(T), y=y, index_label="feature",
                source={"pooled": pooled_cohorts, "not_pooled": [c for c in cohorts if c in NOT_POOLED], "unflagged_pool_excludes": flagged})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "release"))
    a = ap.parse_args()
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--untracked-files=no", "--", "src", "scripts", "configs"], capture_output=True, text=True).stdout.strip()
    if dirty:
        raise SystemExit(f"commit the code and configs first; uncommitted changes:\n{dirty}")
    R = Release(Path(a.out), commit)
    early_detection(R)
    clustering(R)
    signature(R)
    man = pd.DataFrame(R.rows)
    man.to_csv(Path(a.out) / "manifest.csv", index=False)
    print(man[["path", "kind", "layer", "n_rows", "n_features", "harmoniser", "n_LC", "n_COVID_recovered"]].to_string(index=False))


if __name__ == "__main__":
    main()
