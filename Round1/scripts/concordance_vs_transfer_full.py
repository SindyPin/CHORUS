"""Full range of concordance against transfer: the LC points of concordance_vs_transfer_v5 together with positive control
contrasts run through the same leave one study out pipeline.

Points: one per held-out cohort and label-free harmoniser (identity, center_scale, quantile_ref, fsqn, rank), invariant
selector, logreg, from loso_results.csv (no error rows): auroc, effect_spearman (agreement of the held-out cohort's effect
with the training effect over the selected features) and sign_agreement.
  LC        the pooled blocks of results/concordance_vs_transfer_v5 (transcriptomics v4, affinity v5, MS v4, metabolomics)
  controls  acute versus never infected, proteomics (results/control_acute_proteomics); inferred sex, transcriptomics
            (results/control_sex_label/all_genes). results/control_sex_label/no_xy is reported as a sensitivity set.
            The acute versus never infected transcriptomic control was never run (one cohort qualifies,
            results/control_acute_transcriptomics/cohorts.csv), so it contributes no point.

Outputs (--out):
  points.csv        every point with its contrast, block and group (LC or control)
  correlation.csv   Spearman of AUROC with each concordance metric for LC alone, controls alone, LC and controls, and LC
                    and controls with no_xy in place of all_genes: cohort bootstrap interval and permutation p
                    (scripts/concordance_vs_transfer.correlation, resampling within block)
  monotone.csv      one increasing isotonic fit of AUROC on the metric over LC and control points together: R2 of the
                    joint fit; mean residual of the LC and of the control points; the reduction in squared error when LC
                    and controls get separate isotonic fits, with a permutation p (group labels shuffled, 2,000 times);
                    the range of the metric in each group and the LC points inside the control range
  curve.csv         the joint isotonic fit evaluated at every point

  python scripts/concordance_vs_transfer_full.py --out results/concordance_vs_transfer_full
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression

sys.path.insert(0, str(Path(__file__).resolve().parent))
from concordance_vs_transfer import METRICS, correlation, load_rows  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LC_BLOCKS = {
    "transcriptomics": "results/lc_transcriptomics_blood_v4/loso_results.csv",
    "proteomics_affinity": "results/lc_proteomics_affinity_v5/loso_results.csv",
    "proteomics_ms": "results/lc_proteomics_ms_v4/loso_results.csv",
    "metabolomics": "results/lc_metabolomics/loso_results.csv",
}
CONTROLS = {
    "acute_proteomics": "results/control_acute_proteomics/loso_results.csv",
    "sex_transcriptomics": "results/control_sex_label/all_genes/loso_results.csv",
}
SENSITIVITY = {"sex_transcriptomics_no_xy": "results/control_sex_label/no_xy/loso_results.csv"}


def sse(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return 0.0
    f = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(x, y)
    return float(((y - f.predict(x)) ** 2).sum())


def monotone(P: pd.DataFrame, metric: str, n_perm: int, seed: int) -> tuple:
    d = P.dropna(subset=[metric, "auroc"])
    x, y, g = d[metric].to_numpy(), d["auroc"].to_numpy(), (d["group"] == "LC").to_numpy()
    iso = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(x, y)
    fit = iso.predict(x)
    res = y - fit
    joint = float((res**2).sum())
    split = sse(x[g], y[g]) + sse(x[~g], y[~g])
    stat = joint - split
    rng = np.random.default_rng(seed)
    null = np.array([(lambda gg: joint - sse(x[gg], y[gg]) - sse(x[~gg], y[~gg]))(rng.permutation(g)) for _ in range(n_perm)])
    lo_c, hi_c = x[~g].min(), x[~g].max()
    inside = g & (x >= lo_c) & (x <= hi_c)
    row = {
        "metric": metric, "n_LC": int(g.sum()), "n_control": int((~g).sum()),
        "r2_joint": 1 - joint / float(((y - y.mean()) ** 2).sum()),
        "mean_residual_LC": float(res[g].mean()), "mean_residual_control": float(res[~g].mean()),
        "sse_joint": joint, "sse_separate": split, "sse_reduction_separate": stat, "p_permutation_separate": float((1 + (null >= stat).sum()) / (1 + n_perm)),
        "LC_metric_min": float(x[g].min()), "LC_metric_max": float(x[g].max()), "control_metric_min": float(lo_c), "control_metric_max": float(hi_c),
        "n_LC_inside_control_range": int(inside.sum()), "mean_residual_LC_inside_control_range": float(res[inside].mean()) if inside.any() else np.nan,
    }
    curve = d.assign(metric=metric, fitted=fit, residual=res)[["contrast", "block", "group", "held_out", "harmoniser", "metric", metric, "auroc", "fitted", "residual"]].rename(columns={metric: "value"})
    return row, curve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-perm-curve", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    parts = [load_rows(LC_BLOCKS).assign(group="LC", contrast="LC_vs_recovered")]
    for name, path in {**CONTROLS, **SENSITIVITY}.items():
        parts.append(load_rows({name: path}).assign(group="control" if name in CONTROLS else "sensitivity", contrast=name))
    P = pd.concat(parts, ignore_index=True)
    P.to_csv(out / "points.csv", index=False)
    print(P.groupby(["group", "block"]).agg(n=("auroc", "size"), cohorts=("held_out", "nunique"), auroc=("auroc", "median"), es=("effect_spearman", "median")).round(3).to_string())
    analyses = {
        "LC": P[P["group"] == "LC"],
        "controls": P[P["group"] == "control"],
        "LC_and_controls": P[P["group"].isin(["LC", "control"])],
        "LC_and_controls_no_xy": P[(P["group"] == "LC") | (P["block"].isin(["acute_proteomics", "sex_transcriptomics_no_xy"]))],
    }
    rows = [{"analysis": n, **correlation(d, m, a.n_boot, a.n_perm, a.seed)} for n, d in analyses.items() for m in METRICS]
    C = pd.DataFrame(rows)
    C.to_csv(out / "correlation.csv", index=False)
    mono, curves = [], []
    for name in ("LC_and_controls", "LC_and_controls_no_xy"):
        d = analyses[name].assign(group=lambda t: np.where(t["group"] == "LC", "LC", "control"))
        for m in METRICS:
            r, c = monotone(d, m, a.n_perm_curve, a.seed)
            mono.append({"analysis": name, **r})
            curves.append(c.assign(analysis=name))
    Mo = pd.DataFrame(mono)
    Mo.to_csv(out / "monotone.csv", index=False)
    pd.concat(curves, ignore_index=True).to_csv(out / "curve.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"lc_blocks": LC_BLOCKS, "controls": CONTROLS, "sensitivity": SENSITIVITY, "not_available": {"acute_transcriptomics": "not run: one cohort qualifies (results/control_acute_transcriptomics/cohorts.csv)"},
                                                   "n_boot": a.n_boot, "n_perm": a.n_perm, "n_perm_curve": a.n_perm_curve, "seed": a.seed}, indent=1))
    print(C[["analysis", "metric", "n_rows", "n_cohorts", "spearman", "ci_cohort_lo", "ci_cohort_hi", "p_permutation"]].round(3).to_string(index=False))
    print(Mo.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
