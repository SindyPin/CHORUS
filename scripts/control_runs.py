"""Positive and permutation controls for the primary leave-one-study-out benchmark.

matrix   Positive control: the label is sample_matrix (whole_blood = positive, PBMC = negative)
         instead of lc_status. GSE226260_WB and GSE226260_PBMC are merged into one held-out study,
         GSE226260, because a cohort with a single matrix gives a one-class test set; every other
         cohort is training only. Grid: identity, center_scale, combat_addon and chorus; variance
         and invariant selectors; logreg. sample_matrix is removed from the covariates and from
         protect, because it is the label here.

permute  Label permutation control: lc_status is shuffled within each cohort among the samples in
         the LC versus COVID_recovered contrast (class counts per cohort unchanged), then the
         primary leave-one-study-out benchmark is run with identity and chorus, the invariant
         selector and logreg, for each permutation. The summary gives the mean AUROC per held-out
         cohort and harmoniser over permutations.

Both start from the cohorts, preprocessing and covariates of the given config (load_config).

  python scripts/control_runs.py matrix configs/lc_transcriptomics_blood.yaml --out results/control_matrix_label
  python scripts/control_runs.py permute configs/lc_transcriptomics_blood.yaml --out results/control_label_permutation --n-perm 20
"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

import numpy as np
import pandas as pd

from chorus.cli import load_config
from chorus.eval import run_loso, summarise
from chorus.io import Cohort, CohortSet

MERGE = ("GSE226260_WB", "GSE226260_PBMC")


def matrix_control(config: str, out: Path) -> None:
    cs, lc, _ = load_config(config)
    merged, rest = [cs[s] for s in MERGE], [c for c in cs.cohorts if c.study_id not in MERGE]
    X = pd.concat([c.X for c in merged])
    meta = pd.concat([c.meta for c in merged]).assign(study_id="GSE226260")
    cohorts = [Cohort("GSE226260", cs.layer, X, meta, merged[0].feature_space, merged[0].is_counts, {"accession": "GSE226260"})]
    for c in rest:
        cohorts.append(Cohort(c.study_id, c.layer, c.X, c.meta.copy(), c.feature_space, c.is_counts, {**c.source, "train_only": "true"}))
    for c in cohorts:
        c.meta["lc_status"] = c.meta["sample_matrix"]  # the control label
    cs2 = CohortSet(cohorts)
    print(cs2.summary().to_string(index=False))
    cfg = dataclasses.replace(
        lc,
        positive="whole_blood",
        negative="PBMC",
        harmonisers=["identity", "center_scale", "combat_addon", "chorus"],
        selectors=[s for s in lc.selectors if s["name"] in ("variance", "invariant")],
        classifiers=["logreg"],
        covariates=[c for c in lc.covariates if c != "sample_matrix"],
        protect=[c for c in lc.protect if c != "sample_matrix"],
        nested_select_harmoniser=False,
    )
    res, preds, _ = run_loso(cs2, cfg)
    out.mkdir(parents=True, exist_ok=True)
    cs2.summary().to_csv(out / "cohorts.csv", index=False)
    res.to_csv(out / "loso_results.csv", index=False)
    preds.to_csv(out / "loso_predictions.csv", index=False)
    summarise(res).to_csv(out / "loso_summary.csv")
    print(res.drop(columns=["inner_scores"], errors="ignore").to_string(index=False))


def permutation_control(config: str, out: Path, n_perm: int, seed: int) -> None:
    cs, lc, _ = load_config(config)
    cfg = dataclasses.replace(
        lc,
        harmonisers=["identity", "chorus"],
        selectors=[s for s in lc.selectors if s["name"] == "invariant"],
        classifiers=["logreg"],
        nested_select_harmoniser=False,
    )
    rng = np.random.default_rng(seed)
    runs = []
    for i in range(n_perm):
        cohorts = []
        for c in cs.cohorts:
            meta = c.meta.copy()
            m = meta["lc_status"].isin([lc.positive, lc.negative])
            meta.loc[m, "lc_status"] = rng.permutation(meta.loc[m, "lc_status"].values)
            cohorts.append(Cohort(c.study_id, c.layer, c.X, meta, c.feature_space, c.is_counts, dict(c.source)))
        res, _, _ = run_loso(CohortSet(cohorts), cfg, verbose=False)
        res.insert(0, "permutation", i)
        runs.append(res)
        print(f"permutation {i}: " + ", ".join(f"{r.held_out} {r.harmoniser} {r.auroc:.3f}" for r in res.itertuples()), flush=True)
    res = pd.concat(runs, ignore_index=True)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "permutation_results.csv", index=False)
    summ = (
        res.groupby(["held_out", "harmoniser", "selector", "classifier"])["auroc"]
        .agg(["mean", "std", "min", "max", "count"])
        .rename(columns={"count": "n_permutations"})
        .reset_index()
    )
    summ.to_csv(out / "permutation_summary.csv", index=False)
    print(summ.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("control", choices=["matrix", "permute"])
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    if a.control == "matrix":
        matrix_control(a.config, Path(a.out))
    else:
        permutation_control(a.config, Path(a.out), a.n_perm, a.seed)


if __name__ == "__main__":
    main()
