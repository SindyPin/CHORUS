"""Positive and permutation controls for the primary leave-one-study-out benchmark.

sex      Positive control: the label is inferred_sex (female = positive, male = negative;
         ambiguous samples left out), from XIST and RPS4Y1 inside each cohort
         (scripts/infer_sex.py), on the samples of the primary LC versus COVID_recovered
         contrast. Grid: identity, center_scale, combat_addon, chorus and chorus_protected; variance and
         invariant selectors; logreg. Run twice: on the filtered genes of the config
         (all_genes) and on the same genes without those on chrX or chrY (no_xy, chromosomes
         from the GENCODE v47 tx2gene table). This replaces the earlier sample_matrix control.
         A third set (forced_sex_genes) adds XIST, RPS4Y1 and the expressed chrY genes to
         all_genes: a chrY gene (chrY only, not pseudoautosomal) counts as expressed in a cohort
         when its log2 CPM is at least log2(min_cpm) in at least min_fraction_samples of the
         cohort's inferred male samples (the expression filter of the config applied to males);
         elsewhere it is set to missing. Per fold a forced gene is used only when every training
         cohort and the held-out cohort have it (the held-out cohort must be scored), so
         LosoConfig.drop_incomplete_features is set for this set.

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


def forced_sex_genes(cs: CohortSet, chroms: pd.Series, min_cpm: float = 1.0, min_fraction: float = 0.5):
    """Mask the chrY only genes that are not expressed in a cohort's inferred males.

    Returns (cohort set with those values set to missing, candidate genes, presence table with
    one row per gene and cohort: in_matrix, expressed_in_males, fraction of males at or above
    the threshold, n_males)."""
    ychr = [g for g in chroms.index[chroms == "chrY"]]
    cand = [XIST, RPS4Y1] + [g for g in ychr if g != RPS4Y1]
    cand = [g for g in cand if any(g in c.features and c.X[g].notna().any() for c in cs.cohorts)]
    rows, cohorts = [], []
    for c in cs.cohorts:
        X = c.X.copy()
        males = c.meta.index[c.meta["inferred_sex"] == "male"]
        for g in cand:
            present = g in X.columns and bool(X[g].notna().all())
            frac = float((X.loc[males, g] >= np.log2(min_cpm)).mean()) if present and len(males) else np.nan
            expressed = g in (XIST, RPS4Y1) or (present and frac >= min_fraction)
            if g in X.columns and not (present and expressed):
                X[g] = np.nan
            rows.append({"study_id": c.study_id, "gene_id": g, "in_matrix": present, "frac_males_expressed": frac, "n_males": len(males), "used": present and expressed})
        cohorts.append(Cohort(c.study_id, c.layer, X, c.meta, c.feature_space, c.is_counts, dict(c.source)))
    return CohortSet(cohorts), cand, pd.DataFrame(rows)


def sex_control(config: str, out: Path, gene_table: str) -> None:
    chroms = gene_chromosomes(gene_table)
    ychr = list(chroms.index[chroms == "chrY"])
    cs, lc, cfg = load_config(config, extra_features=[XIST, RPS4Y1] + [g for g in ychr if g != RPS4Y1])
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
    xy = set(sex_chromosome_genes(chroms, feats))
    pre = cfg.get("preprocessing", {})
    cs2, forced, presence = forced_sex_genes(cs2, chroms, pre.get("min_cpm", 1.0), pre.get("min_fraction_samples", 0.5))
    forced = [g for g in forced if g not in set(feats)]
    feature_sets = {"all_genes": feats, "no_xy": feats[~feats.isin(xy)], "forced_sex_genes": feats.append(pd.Index(forced))}
    # per fold use of the forced genes: present in every training cohort and in the held-out cohort
    used = presence.pivot(index="gene_id", columns="study_id", values="used").reindex(forced)
    fold_rows = []
    for held_out in cs2.study_ids:
        if counts.set_index("study_id").reindex([held_out]).reindex(columns=["female", "male"]).fillna(0).min(axis=1).iloc[0] < lc.min_test_per_class:
            continue
        train = [s for s in cs2.study_ids if s != held_out]
        for g in forced:
            fold_rows.append({"held_out": held_out, "gene_id": g, "in_every_training_cohort": bool(used.loc[g, train].all()), "in_held_out": bool(used.loc[g, held_out]), "used": bool(used.loc[g, train].all() and used.loc[g, held_out])})
    folds = pd.DataFrame(fold_rows)
    on = chroms
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
                "n_forced": len([g for g in f if g in set(forced)]),
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
        harmonisers=["identity", "center_scale", "combat_addon", "chorus", "chorus_protected"],
        selectors=sel,
        classifiers=["logreg"],
        nested_select_harmoniser=False,
    )
    out.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out / "cohorts.csv", index=False)
    info.to_csv(out / "feature_sets.csv", index=False)
    sym = _symbols()
    presence.insert(2, "symbol", presence["gene_id"].map(sym))
    presence.insert(3, "chrom", presence["gene_id"].map(chroms))
    presence.to_csv(out / "forced_genes_by_cohort.csv", index=False)
    folds.insert(2, "symbol", folds["gene_id"].map(sym))
    folds.to_csv(out / "forced_genes_by_fold.csv", index=False)
    print(folds.groupby("held_out")["used"].sum().to_string())
    allres = []
    for name, f in feature_sets.items():
        res, preds, selected = run_loso(cs2, dataclasses.replace(cfg2, drop_incomplete_features=name == "forced_sex_genes"), features=f)
        d = out / name
        d.mkdir(exist_ok=True)
        res.to_csv(d / "loso_results.csv", index=False)
        preds.to_csv(d / "loso_predictions.csv", index=False)
        (d / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in selected.items()}, indent=1))
        if "auroc" in res.columns:
            summarise(res).to_csv(d / "loso_summary.csv")
        allres.append(res.assign(feature_set=name))
    res = pd.concat(allres, ignore_index=True)
    cols = ["feature_set", "held_out", "harmoniser", "selector", "classifier", "n_features", "n_features_dropped_incomplete", "auroc", "auroc_lo", "auroc_hi", "error"]
    res.reindex(columns=cols).to_csv(out / "auroc_by_held_out.csv", index=False)
    print(res.reindex(columns=cols).to_string(index=False))


def _symbols(path: str = "data/ref/hgnc_complete_set.txt") -> pd.Series:
    if not Path(path).exists():
        return pd.Series(dtype=str)
    t = pd.read_csv(path, sep="\t", usecols=["symbol", "ensembl_gene_id"], dtype=str).dropna()
    return t.drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id")["symbol"]


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
