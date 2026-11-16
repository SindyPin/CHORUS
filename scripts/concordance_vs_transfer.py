"""Concordance and transfer: does the held-out AUROC of the leave one study out benchmark go with the agreement of
the held-out cohort's effects with the training effects?

Rows: every held-out cohort of the blocks below with the invariant selector, logreg and the label free harmonisers
(identity, center_scale, quantile_ref, fsqn, rank: fitted without labels and applied to the held-out cohort in addon
mode; quantile_ref was run in the transcriptomic block only; combat_addon and chorus use the training labels and
combat_pooled is transductive, so they are left out). Error rows are dropped. Per row, from loso_results.csv:
auroc, effect_spearman (Spearman correlation over the selected features between the pooled training effect and the
held-out cohort's effect after harmonisation) and sign_agreement (fraction of the selected features with the same
effect sign); both use the held-out labels and are diagnostics, not inputs of the model.

  rows.csv         block, held_out, harmoniser, auroc, effect_spearman, sign_agreement, n_genes
  correlation.csv  per block and pooled, per concordance metric: n_rows, n_cohorts, Spearman rho of AUROC with the
                   metric, two bootstrap intervals (resampling held-out cohorts with all their rows, within block for
                   the pooled analysis: ci_cohort_*; resampling rows: ci_row_*) and a permutation p value (the
                   metric shuffled across rows, within block for the pooled analysis; two sided on |rho|)
  variance.csv     per block and pooled: two way ANOVA of AUROC on held-out cohort and harmoniser (type II sums of
                   squares); share = sum of squares / total sum of squares, with the residual

  python scripts/concordance_vs_transfer.py --out results/concordance_vs_transfer
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = {
    "transcriptomics": "results/lc_transcriptomics_blood_v4/loso_results.csv",
    "proteomics_affinity": "results/lc_proteomics_affinity_v4/loso_results.csv",
    "proteomics_ms": "results/lc_proteomics_ms_v4/loso_results.csv",
    "metabolomics": "results/lc_metabolomics/loso_results.csv",
}
LABEL_FREE = ["identity", "center_scale", "quantile_ref", "fsqn", "rank"]
METRICS = ["effect_spearman", "sign_agreement"]


def load_rows(blocks: dict = BLOCKS, root: Path = ROOT, selector: str = "invariant", classifier: str = "logreg") -> pd.DataFrame:
    parts = []
    for block, path in blocks.items():
        r = pd.read_csv(root / path)
        if "error" in r.columns:
            r = r[r["error"].isna()]
        r = r[(r["selector"] == selector) & (r["classifier"] == classifier) & r["harmoniser"].isin(LABEL_FREE)]
        parts.append(r.assign(block=block)[["block", "held_out", "harmoniser", "auroc", *METRICS, "n_features"]])
    out = pd.concat(parts, ignore_index=True)
    return out.dropna(subset=["auroc", *METRICS])


def _rho(a, b) -> float:
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(stats.spearmanr(a, b).statistic)


def correlation(df: pd.DataFrame, metric: str, n_boot: int = 2000, n_perm: int = 10000, seed: int = 0) -> dict:
    """Spearman rho of auroc with ``metric``; cohort and row bootstrap intervals and a permutation p value, every
    resampling stratified by block (one block: unstratified)."""
    rng = np.random.default_rng(seed)
    a, m = df["auroc"].to_numpy(), df[metric].to_numpy()
    rho = _rho(a, m)
    blocks = df["block"].to_numpy()
    by_block = {b: np.where(blocks == b)[0] for b in pd.unique(blocks)}
    cohort_rows = {b: [np.where((blocks == b) & (df["held_out"].to_numpy() == c))[0] for c in pd.unique(df.loc[blocks == b, "held_out"])] for b in by_block}
    boot_c, boot_r = [], []
    for _ in range(n_boot):
        idx = np.concatenate([np.concatenate([rows[i] for i in rng.integers(0, len(rows), len(rows))]) for rows in cohort_rows.values()])
        boot_c.append(_rho(a[idx], m[idx]))
        idx = np.concatenate([rng.choice(ix, len(ix)) for ix in by_block.values()])
        boot_r.append(_rho(a[idx], m[idx]))
    perm = np.empty(n_perm)
    for i in range(n_perm):
        mp = m.copy()
        for ix in by_block.values():
            mp[ix] = rng.permutation(m[ix])
        perm[i] = _rho(a, mp)
    ci = lambda v: np.nanpercentile(v, [2.5, 97.5]) if np.isfinite(v).any() else (np.nan, np.nan)  # noqa: E731
    (cl, ch), (rl, rh) = ci(np.array(boot_c)), ci(np.array(boot_r))
    p = (1 + np.sum(np.abs(perm[np.isfinite(perm)]) >= abs(rho))) / (1 + np.isfinite(perm).sum()) if np.isfinite(rho) else np.nan
    return {
        "metric": metric,
        "n_rows": len(df),
        "n_cohorts": int(df.groupby("block")["held_out"].nunique().sum()),
        "spearman": rho,
        "ci_cohort_lo": cl,
        "ci_cohort_hi": ch,
        "ci_row_lo": rl,
        "ci_row_hi": rh,
        "p_permutation": p,
        "n_boot": n_boot,
        "n_perm": n_perm,
    }


def variance_shares(df: pd.DataFrame) -> pd.DataFrame:
    """Type II two way ANOVA of auroc on held-out cohort and harmoniser; each sum of squares over the total."""
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm

    d = df.assign(cohort=df["block"] + ":" + df["held_out"])
    fit = smf.ols("auroc ~ C(cohort) + C(harmoniser)", data=d).fit()
    a = anova_lm(fit, typ=2)
    total = float(((d["auroc"] - d["auroc"].mean()) ** 2).sum())
    names = {"C(cohort)": "held_out_cohort", "C(harmoniser)": "harmoniser", "Residual": "residual"}
    return pd.DataFrame(
        [{"term": names[t], "df": float(a.loc[t, "df"]), "sum_sq": float(a.loc[t, "sum_sq"]), "share": float(a.loc[t, "sum_sq"]) / total, "F": a.loc[t, "F"], "p": a.loc[t, "PR(>F)"]} for t in a.index]
    ).assign(total_sum_sq=total, n_rows=len(d), n_cohorts=d["cohort"].nunique(), n_harmonisers=d["harmoniser"].nunique())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    rows.to_csv(out / "rows.csv", index=False)
    corr, var = [], []
    analyses = {**{b: rows[rows["block"] == b] for b in BLOCKS}, "pooled": rows}
    for name, df in analyses.items():
        for metric in METRICS:
            corr.append({"analysis": name, **correlation(df, metric, a.n_boot, a.n_perm, a.seed)})
        var.append(variance_shares(df).assign(analysis=name))
    C = pd.DataFrame(corr)
    V = pd.concat(var, ignore_index=True)
    V = V[["analysis", *[c for c in V.columns if c != "analysis"]]]
    C.to_csv(out / "correlation.csv", index=False)
    V.to_csv(out / "variance.csv", index=False)
    (out / "run_info.json").write_text(json.dumps({"blocks": BLOCKS, "harmonisers": LABEL_FREE, "selector": "invariant", "classifier": "logreg", "n_boot": a.n_boot, "n_perm": a.n_perm, "seed": a.seed}, indent=1))
    print(C.round(3).to_string(index=False))
    print(V.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
