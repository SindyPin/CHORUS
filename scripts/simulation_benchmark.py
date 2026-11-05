"""Run the CHORUS benchmark on simulated multi-study data and draw the summary figure.

python scripts/simulation_benchmark.py --out results/simulation
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from chorus.eval import LosoConfig, run_loso, run_within_study_cv, summarise
from chorus.simulate import simulate_cohorts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/simulation")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-studies", type=int, default=6)
    ap.add_argument("--n-features", type=int, default=2000)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cs = simulate_cohorts(n_studies=a.n_studies, n_features=a.n_features, n_causal=30, n_shortcut=60, seed=a.seed)
    cs.summary().to_csv(out / "cohorts.csv", index=False)
    cfg = LosoConfig(
        harmonisers=["identity", "center_scale", "quantile_ref", "fsqn", "rank", "combat_addon"],
        selectors=[{"name": "variance", "k": 300}, {"name": "differential", "k": 60}, {"name": "invariant", "k": 60}],
        classifiers=["logreg", "rf", "cross_study"],
        covariates=["age", "sex"],
        nested_select_harmoniser=True,
        n_boot=100,
        random_state=a.seed,
    )
    res, preds, selected = run_loso(cs, cfg, verbose=True)
    res.to_csv(out / "loso_results.csv", index=False)
    preds.to_csv(out / "loso_predictions.csv", index=False)
    summ = summarise(res)
    summ.to_csv(out / "loso_summary.csv")
    within = run_within_study_cv(cs, cfg)
    within.to_csv(out / "within_study_cv.csv", index=False)
    print(summ.to_string())
    print("\nwithin-study CV mean AUROC:")
    print(within.groupby(["selector", "classifier"])["auroc"].mean().round(3).to_string())
    draw(res, within, out / "figure_simulation.png")


def draw(res: pd.DataFrame, within: pd.DataFrame, path: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ok = res[res["classifier"].notna() & res.get("error", pd.Series(index=res.index)).isna()] if "error" in res.columns else res
    ok = ok[ok["classifier"] == "logreg"]
    order_h = ["identity", "center_scale", "quantile_ref", "fsqn", "rank", "combat_addon", "nested"]
    order_s = ["variance", "differential", "invariant"]
    piv = ok.pivot_table(index="harmoniser", columns="selector", values="auroc", aggfunc="mean").reindex(order_h)[order_s]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
    colours = {"variance": "#8da0cb", "differential": "#fc8d62", "invariant": "#66c2a5"}
    x = range(len(piv))
    w = 0.26
    for i, s in enumerate(order_s):
        axes[0].bar([xi + (i - 1) * w for xi in x], piv[s].values, width=w, label=s, color=colours[s])
    axes[0].axhline(0.5, color="grey", lw=0.8, ls="--")
    wm = within[within["classifier"] == "logreg"]["auroc"].mean()
    axes[0].axhline(wm, color="black", lw=0.8, ls=":", label=f"within-study CV mean ({wm:.2f})")
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(piv.index, rotation=30, ha="right")
    axes[0].set_ylabel("Leave-one-study-out AUROC (mean over held-out studies)")
    axes[0].set_ylim(0.4, 1.0)
    axes[0].set_title("A. Harmoniser x feature selector (L2 logistic regression)")
    axes[0].legend(frameon=False, fontsize=8)
    sub = ok[(ok["harmoniser"].isin(["identity", "center_scale"])) & (ok["selector"] == "invariant")]
    for h, c in [("identity", "#bdbdbd"), ("center_scale", "#66c2a5")]:
        d = sub[sub["harmoniser"] == h].sort_values("held_out")
        axes[1].errorbar(range(len(d)), d["auroc"], yerr=[d["auroc"] - d["auroc_lo"], d["auroc_hi"] - d["auroc"]], fmt="o", color=c, label=h, capsize=3)
    axes[1].set_xticks(range(len(d)))
    axes[1].set_xticklabels(d["held_out"])
    axes[1].axhline(0.5, color="grey", lw=0.8, ls="--")
    axes[1].set_ylim(0.3, 1.0)
    axes[1].set_title("B. Per held-out study, invariant features")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    print(f"figure written to {path}")


if __name__ == "__main__":
    main()
