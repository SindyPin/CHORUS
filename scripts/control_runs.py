"""Positive and permutation controls for the primary leave-one-study-out benchmark.

sex      Positive control: the label is inferred_sex (female = positive, male = negative;
         ambiguous samples left out), from XIST and RPS4Y1 inside each cohort
         (scripts/infer_sex.py), on the samples of the primary LC versus COVID_recovered
         contrast. Grid: identity, center_scale, combat_addon and chorus; variance and
         invariant selectors; logreg. Run twice: on the filtered genes of the config
         (all_genes) and on the same genes without those on chrX or chrY (no_xy, chromosomes
         from the GENCODE v47 tx2gene table). This replaces the earlier sample_matrix control.

permute  Label permutation control: lc_status is shuffled within each cohort among the samples in
         the LC versus COVID_recovered contrast (class counts per cohort unchanged), then the
         primary leave-one-study-out benchmark is run with identity and chorus, the invariant
         selector and logreg, for each permutation. The summary gives the mean AUROC per held-out
         cohort and harmoniser over permutations.

Both start from the cohorts, preprocessing and covariates of the given config (load_config).

  python scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label
  python scripts/control_runs.py permute configs/lc_transcriptomics_blood.yaml --out results/control_label_permutation --n-perm 20
"""

from __future__ import annotations

import argparse
import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd

from chorus.cli import load_config
from chorus.eval import run_loso, summarise
from chorus.io import Cohort, CohortSet
from chorus.sex import RPS4Y1, XIST, gene_chromosomes, sex_chromosome_genes

GENE_TABLE = "data/quant/GSE224615/tx2gene_gencode_v47.tsv"


def sex_control(config: str, out: Path, gene_table: str) -> None:
    cs, lc, cfg = load_config(config)
    cohorts, counts = [], []
    for c in cs.cohorts:
        if "inferred_sex" not in c.meta.columns:
            raise ValueError(f"{c.study_id}: no inferred_sex column; run scripts/infer_sex.py first")
        meta = c.meta.copy()
        in_contrast = meta["lc_status"].isin([lc.positive, lc.negative])
        meta["lc_status"] = meta["inferred_sex"].where(in_contrast)  # the control label
        counts.append({"study_id": c.study_id, "n_contrast": int(in_contrast.sum()), **meta.loc[in_contrast, "inferred_sex"].value_counts().to_dict()})
        cohorts.append(Cohort(c.study_id, c.layer, c.X, meta, c.feature_space, c.is_counts, dict(c.source)))
    cs2 = CohortSet(cohorts)
    counts = pd.DataFrame(counts).fillna(0)
    print(counts.to_string(index=False))

    feats = pd.Index(cfg["filtered_features"])
    chroms = gene_chromosomes(gene_table)
    xy = set(sex_chromosome_genes(chroms, feats))
    feature_sets = {"all_genes": feats, "no_xy": feats[~feats.isin(xy)]}
    on = chroms.reindex(feats)
    info = pd.DataFrame(
        [
            {
                "feature_set": name,
                "n_features": len(f),
                "n_chrX": int(on.reindex(f).str.contains("chrX", na=False).sum()),
                "n_chrY": int(on.reindex(f).str.contains("chrY", na=False).sum()),
                "n_chromosome_unknown": int(on.reindex(f).isna().sum()),
                "XIST_in": XIST in f,
                "RPS4Y1_in": RPS4Y1 in f,
            }
            for name, f in feature_sets.items()
        ]
    )
    print(info.to_string(index=False))
    sel = [s for s in lc.selectors if s["name"] in ("variance", "invariant")]
    cfg2 = dataclasses.replace(
        lc,
        positive="female",
        negative="male",
        harmonisers=["identity", "center_scale", "combat_addon", "chorus"],
        selectors=sel,
        classifiers=["logreg"],
        nested_select_harmoniser=False,
    )
    out.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out / "cohorts.csv", index=False)
    info.to_csv(out / "feature_sets.csv", index=False)
    allres = []
    for name, f in feature_sets.items():
        res, preds, selected = run_loso(cs2, cfg2, features=f)
        d = out / name
        d.mkdir(exist_ok=True)
        res.to_csv(d / "loso_results.csv", index=False)
        preds.to_csv(d / "loso_predictions.csv", index=False)
        (d / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in selected.items()}, indent=1))
        if "auroc" in res.columns:
            summarise(res).to_csv(d / "loso_summary.csv")
        allres.append(res.assign(feature_set=name))
    res = pd.concat(allres, ignore_index=True)
    cols = ["feature_set", "held_out", "harmoniser", "selector", "classifier", "n_features", "auroc", "auroc_lo", "auroc_hi", "error"]
    res.reindex(columns=cols).to_csv(out / "auroc_by_held_out.csv", index=False)
    print(res.reindex(columns=cols).to_string(index=False))


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
    ap.add_argument("control", choices=["sex", "permute"])
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--gene-table", default=GENE_TABLE, help="tx2gene table or GTF with gene chromosomes")
    a = ap.parse_args()
    if a.control == "sex":
        sex_control(a.config, Path(a.out), a.gene_table)
    else:
        permutation_control(a.config, Path(a.out), a.n_perm, a.seed)


if __name__ == "__main__":
    main()
