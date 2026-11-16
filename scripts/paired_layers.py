"""Test 2: several omics layers measured on the same people, evaluated inside each paired cohort.

For each paired cohort (PAIRED below) and each contrast with both classes on the paired samples:

  samples.csv                 the paired samples used (sample_id, subject_id, label, timepoint)
  fold_metrics.csv            per model, repeat and fold: metrics on the test fold
                              (chorus.eval.paired.PairedCV; the same folds for every model)
  summary.csv                 per model: mean of each metric over folds, auroc_sd, n_folds,
                              n_failed_folds
  predictions.csv             out of fold predictions: sample, repeat, fold, model, y, p, error
  probability_correlation.csv per repeat and pair of models, Spearman and Pearson correlation
                              across subjects of the out of fold probabilities
  effects.csv, effect_correlation.csv   Hedges g of the contrast per feature and layer, and the
                              Spearman correlation between layers over shared features (only when
                              two layers share a feature space, Cervia2024), with
                              fewer_than_10_shared
  python scripts/paired_layers.py --cohort Cervia2024 --out results/paired_Cervia2024 --pairs-only
  (--pairs-only rewrites effect_correlation.csv from the effects.csv already in --out/<contrast>)
  features.csv                per layer: features in the cohort, dropped for missingness, used

Contrasts: ``lc`` is LC versus COVID_recovered on post-acute samples; ``acute`` is acute (any later
outcome) versus never_infected (the positive control). One sample per subject: the earliest by
``timepoint_order``. A contrast is skipped when the paired samples have fewer than ``n_splits``
samples in a class.

  python scripts/paired_layers.py --cohort Cervia2024 --out results/paired_Cervia2024
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from chorus.eval.paired import PairedCV, effect_correlation, layer_pairs, probability_correlation, summarise_folds

ROOT = Path(__file__).resolve().parents[1]

# layer name -> (matrix file, metadata file, imputation inside folds). Kovarik2023 is not listed: its
# eicosanoid, Biocrates and triglyceride tables can be linked to the proteomics runs only through the
# group_index column labels of the supplement, and those labels do not link two of the supplement's own
# tables (docs/paired_cohorts.md).
PAIRED = {
    "Cervia2024": {
        "layers": {
            "somascan": ("data/proteomics/Cervia2024_SomaScan/Cervia2024_SomaScan_matrix.csv", "data/meta/Cervia2024_SomaScan_meta.csv", "median"),
            "ms": ("data/proteomics/Cervia2024_MS/Cervia2024_MS_matrix.csv", "data/meta/Cervia2024_MS_meta.csv", "min"),
        },
        "shared_space": True,  # both keyed by UniProt
        "timepoint_order": ["acute", "6M"],
    },
    "Su2022": {
        "layers": {
            "olink": ("data/proteomics/Su2022/Su2022_matrix.csv", "data/meta/Su2022_meta.csv", "median"),
            "metabolomics": ("data/metabolomics/Su2022_metabolomics/Su2022_metabolomics_matrix.csv", "data/meta/Su2022_metabolomics_meta.csv", "median"),
        },
        "shared_space": False,
        "timepoint_order": ["T1", "T2", "T3"],
    },
}

CONTRASTS = {
    "lc": {"positive": "LC", "negative": "COVID_recovered"},
    "acute": {"positive": "acute", "negative": "never_infected"},
}


def load_layers(cohort: str, root: Path = ROOT) -> tuple:
    """Matrices (samples x features) and metadata of every layer of the cohort whose files exist."""
    views, metas, impute, missing = {}, {}, {}, []
    for name, (mpath, metapath, how) in PAIRED[cohort]["layers"].items():
        if not (root / mpath).exists():
            missing.append(name)
            continue
        views[name] = pd.read_csv(root / mpath, index_col=0).T
        metas[name] = pd.read_csv(root / metapath, low_memory=False).set_index("sample_id")
        impute[name] = how
    return views, metas, impute, missing


def contrast_label(meta: pd.DataFrame, contrast: str) -> pd.Series:
    """1 / 0 / NaN for one contrast from the curated metadata."""
    st = meta["lc_status"].astype(str)
    phase = meta.get("sample_phase", pd.Series("unknown", index=meta.index)).astype(str)
    y = pd.Series(np.nan, index=meta.index)
    if contrast == "lc":
        post = phase != "acute"
        y[post & (st == "LC")] = 1.0
        y[post & (st == "COVID_recovered")] = 0.0
    else:
        y[phase == "acute"] = 1.0
        y[st == "never_infected"] = 0.0
    return y


def paired_samples(views: dict, metas: dict, contrast: str, timepoint_order: list) -> pd.DataFrame:
    """Samples present in every layer with a label for the contrast (labels from the first layer's
    metadata; they must agree across layers), one per subject, the earliest by timepoint_order."""
    first = next(iter(metas))
    common = None
    for name, V in views.items():
        common = V.index if common is None else common.intersection(V.index)
    meta = metas[first].loc[metas[first].index.intersection(common)]
    y = contrast_label(meta, contrast)
    for name, m in metas.items():  # the labels of the layers must agree on the shared samples
        other = contrast_label(m.reindex(meta.index), contrast)
        if not other.fillna(-1).equals(y.fillna(-1)):
            raise ValueError(f"layer {name} disagrees with {first} on the {contrast} labels")
    meta = meta[y.notna()]
    rank = {t: i for i, t in enumerate(timepoint_order)}
    t = pd.DataFrame({"sample_id": meta.index, "subject_id": meta["subject_id"].astype(str).values, "y": y[meta.index].values,
                      "timepoint": meta.get("timepoint", pd.Series("", index=meta.index)).astype(str).values})
    t["order"] = t["timepoint"].map(rank).fillna(len(rank))
    t = t.sort_values(["order"], kind="stable").groupby("subject_id", sort=False).head(1)
    return t.drop(columns="order").set_index("sample_id").loc[[s for s in meta.index if s in set(t["sample_id"])]]


def run_cohort(cohort: str, out: Path, n_repeats: int = 10, n_splits: int = 5, k: int = 50, models=None) -> dict:
    views, metas, impute, missing = load_layers(cohort)
    out.mkdir(parents=True, exist_ok=True)
    info = {"cohort": cohort, "layers": list(views), "layers_missing": missing, "contrasts": {}}
    if len(views) < 2:
        info["skipped"] = "fewer than two layers on disk"
        (out / "run_info.json").write_text(json.dumps(info, indent=1))
        return info
    for contrast, spec in CONTRASTS.items():
        samples = paired_samples(views, metas, contrast, PAIRED[cohort]["timepoint_order"])
        counts = samples["y"].value_counts()
        n_pos, n_neg = int(counts.get(1.0, 0)), int(counts.get(0.0, 0))
        info["contrasts"][contrast] = {"positive": spec["positive"], "negative": spec["negative"], "n_positive": n_pos, "n_negative": n_neg}
        if min(n_pos, n_neg) < n_splits:
            info["contrasts"][contrast]["skipped"] = f"fewer than {n_splits} paired samples in a class"
            continue
        d = out / contrast
        d.mkdir(exist_ok=True)
        samples.assign(label=np.where(samples["y"] == 1, spec["positive"], spec["negative"])).to_csv(d / "samples.csv")
        y = samples["y"]
        sub = {n: V.loc[y.index] for n, V in views.items()}
        t0 = time.time()
        cv = PairedCV(sub, y, impute=impute, k=k, n_splits=n_splits, n_repeats=n_repeats)
        F, P = cv.run(models=models)
        info["contrasts"][contrast]["seconds"] = round(time.time() - t0, 1)
        F.to_csv(d / "fold_metrics.csv", index=False)
        summarise_folds(F).to_csv(d / "summary.csv", index=False)
        P.to_csv(d / "predictions.csv", index=False)
        probability_correlation(P[P["p"].notna()]).to_csv(d / "probability_correlation.csv", index=False)
        pd.DataFrame(
            [{"layer": n, "n_features_cohort": views[n].shape[1], "n_dropped_missing": cv.dropped_features[n], "n_used": cv.views[n].shape[1], "impute": impute[n]} for n in views]
        ).to_csv(d / "features.csv", index=False)
        if PAIRED[cohort]["shared_space"]:
            E, C = effect_correlation(cv.views, y)
            E.rename_axis("feature").to_csv(d / "effects.csv")
            C.to_csv(d / "effect_correlation.csv", index=False)
        print(summarise_folds(F).to_string(index=False), flush=True)
    (out / "run_info.json").write_text(json.dumps(info, indent=1))
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(PAIRED))
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-repeats", type=int, default=10)
    ap.add_argument("--k", type=int, default=50, help="features per layer kept by the differential selector inside each fold")
    ap.add_argument("--models", help="comma separated subset (default: all available)")
    ap.add_argument("--pairs-only", action="store_true", help="rewrite effect_correlation.csv from the existing effects.csv and stop")
    a = ap.parse_args()
    if a.pairs_only:
        for f in sorted(Path(a.out).glob("*/effects.csv")):
            layer_pairs(pd.read_csv(f, index_col=0)).to_csv(f.parent / "effect_correlation.csv", index=False)
        return
    info = run_cohort(a.cohort, Path(a.out), n_repeats=a.n_repeats, k=a.k, models=a.models.split(",") if a.models else None)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
