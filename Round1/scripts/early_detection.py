"""Early detection baseline: the later outcome, LC versus COVID_recovered, from samples drawn in the acute phase.

Samples: sample_phase acute (GSE215865 has no sample_phase; every library was drawn during the hospital stay) with a
later lc_status of LC or COVID_recovered, one sample per subject, the earliest. Before any model is fitted the subjects
per cohort and group are written to cohorts.csv, and the run stops when a group has fewer than --min-subjects.

horizontal  configs/early_detection_horizontal.yaml: Su2022 (Olink) and Cervia2024_SomaScan on their shared proteins,
            trained on one cohort and tested on the other, both directions (chorus.eval.run_loso with one training
            study). Grid of the config, AUROC with a stratified bootstrap interval. Label permutation control: lc_status
            is shuffled within each cohort (class counts unchanged) and the same grid is rerun --n-perm times;
            permutation_summary.csv gives per configuration the observed AUROC, the mean and SD over permutations and
            p = (1 + permutations with AUROC >= observed) / (1 + n_perm).
reference   configs/early_detection_reference.yaml: GSE215865, the timepoint 0 library of each subject, within cohort
            repeated stratified 5 fold CV (chorus.eval.run_within_study_cv) with the classifiers of the horizontal
            analysis. The published TACO value (0.675) is written beside it as published_taco_auroc.

The vertical analysis (Su2022 acute Olink and metabolomics on the same subjects) is scripts/paired_layers.py with
--contrasts early.

  python scripts/early_detection.py horizontal configs/early_detection_horizontal.yaml --out results/early_detection/horizontal
  python scripts/early_detection.py reference configs/early_detection_reference.yaml --out results/early_detection/reference
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd

from chorus.cli import load_config
from chorus.eval import run_loso, run_within_study_cv, summarise
from chorus.io import Cohort, CohortSet

PUBLISHED_TACO_AUROC = 0.675  # Thompson et al. 2023 (TACO), GSE215865, LC prediction from the acute draw
KEYS = ["held_out", "harmoniser", "selector", "classifier"]


def subject_counts(cs: CohortSet, positive: str, negative: str) -> pd.DataFrame:
    rows = []
    for c in cs.cohorts:
        st = c.meta["lc_status"].astype(str)
        for g in (positive, negative):
            m = st == g
            rows.append({"study_id": c.study_id, "group": g, "n_samples": int(m.sum()), "n_subjects": int(c.meta.loc[m, "subject_id"].nunique())})
    return pd.DataFrame(rows)


def check_counts(cs: CohortSet, lc, out: Path, min_subjects: int) -> pd.DataFrame:
    counts = subject_counts(cs, lc.positive, lc.negative)
    out.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out / "cohorts.csv", index=False)
    print(counts.to_string(index=False), flush=True)
    small = counts[counts["n_subjects"] < min_subjects]
    if len(small):
        raise SystemExit(f"stopped: fewer than {min_subjects} subjects in {small[['study_id', 'group']].values.tolist()}")
    return counts


def permuted(cs: CohortSet, positive: str, negative: str, rng: np.random.Generator) -> CohortSet:
    cohorts = []
    for c in cs.cohorts:
        meta = c.meta.copy()
        m = meta["lc_status"].isin([positive, negative])
        meta.loc[m, "lc_status"] = rng.permutation(meta.loc[m, "lc_status"].values)
        cohorts.append(Cohort(c.study_id, c.layer, c.X, meta, c.feature_space, c.is_counts, dict(c.source)))
    return CohortSet(cohorts)


def permutation_summary(observed: pd.DataFrame, perms: pd.DataFrame) -> pd.DataFrame:
    obs = observed[observed["auroc"].notna()].set_index(KEYS)["auroc"] if "auroc" in observed else pd.Series(dtype=float)
    P = perms[perms["auroc"].notna()] if "auroc" in perms else perms.iloc[0:0]
    rows = []
    for key, g in P.groupby(KEYS):
        a = g["auroc"].to_numpy()
        o = obs.get(key, np.nan)
        rows.append(
            {
                **dict(zip(KEYS, key)),
                "auroc_observed": o,
                "auroc_perm_mean": float(a.mean()),
                "auroc_perm_sd": float(a.std(ddof=1)) if len(a) > 1 else np.nan,
                "n_permutations": len(a),
                "p_permutation": float((1 + (a >= o).sum()) / (1 + len(a))) if np.isfinite(o) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def horizontal(config: str, out: Path, n_perm: int, seed: int, min_subjects: int) -> None:
    cs, lc, cfg = load_config(config)
    print(cs.summary().to_string(index=False))
    check_counts(cs, lc, out, min_subjects)
    feats = pd.Index(cfg["filtered_features"])
    info = {"config": config, "n_features": len(feats), "n_permutations": n_perm, "seed": seed}
    res, preds, selected = run_loso(cs, lc)
    res.to_csv(out / "loso_results.csv", index=False)
    preds.to_csv(out / "loso_predictions.csv", index=False)
    (out / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in selected.items()}, indent=1))
    summarise(res, by=tuple(KEYS)).to_csv(out / "loso_summary.csv")
    rng = np.random.default_rng(seed)
    cfg_p = dataclasses.replace(lc, n_boot=0)
    runs = []
    for i in range(n_perm):
        r, _, _ = run_loso(permuted(cs, lc.positive, lc.negative, rng), cfg_p, verbose=False)
        runs.append(r.assign(permutation=i))
        print(f"permutation {i}: {int(r['auroc'].notna().sum()) if 'auroc' in r else 0} rows", flush=True)
    perms = pd.concat(runs, ignore_index=True)
    perms[["permutation", *KEYS, *[c for c in ("auroc", "error") if c in perms]]].to_csv(out / "permutation_results.csv", index=False)
    ps = permutation_summary(res, perms)
    ps.to_csv(out / "permutation_summary.csv", index=False)
    print(ps.to_string(index=False))
    (out / "run_info.json").write_text(json.dumps(info, indent=1))


def reference(config: str, out: Path, n_repeats: int, min_subjects: int) -> None:
    cs, lc, cfg = load_config(config)
    print(cs.summary().to_string(index=False))
    check_counts(cs, lc, out, min_subjects)
    w, folds = run_within_study_cv(cs, lc, n_repeats=n_repeats, return_folds=True)
    w.insert(3, "published_taco_auroc", PUBLISHED_TACO_AUROC)
    w.to_csv(out / "within_study_cv.csv", index=False)
    folds.to_csv(out / "within_study_cv_folds.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"config": config, "n_features": len(cfg["filtered_features"]), "n_repeats": n_repeats, "published_taco_auroc": PUBLISHED_TACO_AUROC}, indent=1))
    print(w.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("analysis", choices=["horizontal", "reference"])
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--n-repeats", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-subjects", type=int, default=10, help="stop when a cohort has fewer subjects in a group")
    a = ap.parse_args()
    if a.analysis == "horizontal":
        horizontal(a.config, Path(a.out), a.n_perm, a.seed, a.min_subjects)
    else:
        reference(a.config, Path(a.out), a.n_repeats, a.min_subjects)


if __name__ == "__main__":
    main()
