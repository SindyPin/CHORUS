"""Command line interface.

  chorus simulate --out sim/ [--n-studies 6]
  chorus run config.yaml
  chorus summarise results.csv

A config file lists cohorts (matrix, metadata, layer, is_counts) and the benchmark grid:

  layer: bulk_rna
  cohorts:
    - study_id: GSE275334
      matrix: data/GSE275334_counts.tsv
      meta: data/GSE275334_meta.csv
      is_counts: true
      accession: GSE275334
  preprocessing:
    min_cpm: 1
    min_fraction_samples: 0.5
    matrices: [whole_blood]
  benchmark:
    harmonisers: [identity, center_scale, fsqn, rank, combat_addon]
    selectors:
      - {name: variance, k: 500}
      - {name: invariant, k: 200}
      - {name: prior, genes: priors/taco_411.txt, mode: union, k: 200}
    classifiers: [logreg_en, rf, taco, cross_study]
    covariates: [age, sex]
    nested_select_harmoniser: true
  out: results/
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from .eval import LosoConfig, run_loso, run_within_study_cv, summarise
from .features.priors import load_gene_list
from .harmonise.transcriptomics import filter_low_counts, select_tissue, to_continuous
from .io import Cohort, CohortSet
from .simulate import simulate_cohorts


def _cmd_simulate(a):
    cs = simulate_cohorts(n_studies=a.n_studies, seed=a.seed)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = {"layer": "bulk_rna", "cohorts": [], "preprocessing": {"matrices": ["whole_blood"]}, "benchmark": {}, "out": str(out / "results")}
    for c in cs.cohorts:
        c.X.T.to_csv(out / f"{c.study_id}_matrix.tsv", sep="\t")
        c.meta.to_csv(out / f"{c.study_id}_meta.csv", index=False)
        cfg["cohorts"].append({"study_id": c.study_id, "matrix": str(out / f"{c.study_id}_matrix.tsv"), "meta": str(out / f"{c.study_id}_meta.csv"), "is_counts": False, "accession": "simulated"})
    (out / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    print(f"wrote {len(cs.cohorts)} simulated cohorts and {out/'config.yaml'}")


def load_config(path: str) -> tuple:
    cfg = yaml.safe_load(Path(path).read_text())
    cohorts = []
    for c in cfg["cohorts"]:
        if not c.get("enabled", True):  # listed for the record, not loaded (missing matrix or labels)
            continue
        cohorts.append(
            Cohort.from_files(
                c["study_id"],
                cfg.get("layer", "bulk_rna"),
                c["matrix"],
                c["meta"],
                features_in_rows=c.get("features_in_rows", True),
                is_counts=c.get("is_counts", False),
                feature_space=c.get("feature_space", "ensembl_gene"),
                accession=str(c.get("accession", "")),
                provenance=str(c.get("provenance", "")),
                train_only=str(bool(c.get("train_only", False))).lower(),
            )
        )
    cs = CohortSet(cohorts)
    pre = cfg.get("preprocessing", {})
    if pre.get("matrices"):
        cs = select_tissue(cs, pre["matrices"])
    if pre.get("deduplicate", True):
        cs = cs.deduplicate_subjects()
    feats = filter_low_counts(cs, pre.get("min_cpm", 1.0), pre.get("min_fraction_samples", 0.5))
    cs = to_continuous(cs, feats)
    b = dict(cfg.get("benchmark", {}))
    for s in b.get("selectors", []):
        if s.get("name") == "prior" and isinstance(s.get("genes"), str):
            s["genes"] = load_gene_list(s["genes"])
    lc = LosoConfig(**b)
    return cs, lc, cfg


def _cmd_run(a):
    cs, lc, cfg = load_config(a.config)
    out = Path(cfg.get("out", "results"))
    out.mkdir(parents=True, exist_ok=True)
    print(cs.summary().to_string(index=False))
    results, preds, selected = run_loso(cs, lc)
    results.to_csv(out / "loso_results.csv", index=False)
    preds.to_csv(out / "loso_predictions.csv", index=False)
    (out / "selected_features.json").write_text(json.dumps({"|".join(k): v for k, v in selected.items()}, indent=1))
    summ = summarise(results)
    summ.to_csv(out / "loso_summary.csv")
    print("\nLOSO summary (mean AUROC across held-out studies):")
    print(summ.head(20).to_string())
    if a.within:
        w = run_within_study_cv(cs, lc)
        w.to_csv(out / "within_study_cv.csv", index=False)
        print("\nWithin-study CV AUROC (optimism reference):")
        print(w.groupby(["selector", "classifier"])["auroc"].mean().to_string())


def _cmd_summarise(a):
    r = pd.read_csv(a.results)
    print(summarise(r).to_string())


def main(argv=None):
    p = argparse.ArgumentParser(prog="chorus", description="CHORUS cross cohort harmonisation benchmark")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--out", default="sim")
    s.add_argument("--n-studies", type=int, default=6)
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(func=_cmd_simulate)
    r = sub.add_parser("run")
    r.add_argument("config")
    r.add_argument("--within", action="store_true", help="also run within-study CV")
    r.set_defaults(func=_cmd_run)
    m = sub.add_parser("summarise")
    m.add_argument("results")
    m.set_defaults(func=_cmd_summarise)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
