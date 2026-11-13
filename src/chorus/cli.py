"""Command line interface.

  chorus simulate --out sim/ [--n-studies 6]
  chorus run config.yaml
  chorus summarise results.csv
  chorus harmonisation config.yaml --out results/<name>_harmonisation

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
from .eval.loso import selector_label
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


def load_config(path: str, extra_features=()) -> tuple:
    """Load the cohorts and benchmark grid of a config. ``extra_features`` are kept after the
    expression filter when at least one cohort's matrix has them, as NaN in the cohorts without
    them (for diagnostics on marker genes; pass ``features`` explicitly to pooled views). The
    filtered gene set itself is returned as ``cfg["filtered_features"]``."""
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
    if pre.get("one_sample_per_subject", False):
        cs = cs.one_sample_per_subject(pre.get("acute_timepoints", ["acute"]))
    if pre.get("deduplicate", True):
        cs = cs.deduplicate_subjects()
    feats = filter_low_counts(cs, pre.get("min_cpm", 1.0), pre.get("min_fraction_samples", 0.5))
    cfg["filtered_features"] = list(feats)
    extra = [g for g in dict.fromkeys(extra_features) if g not in set(feats) and any(g in c.features for c in cs.cohorts)]
    cs = to_continuous(cs, feats.append(pd.Index(extra)))
    b = dict(cfg.get("benchmark", {}))
    for s in b.get("selectors", []):
        if s.get("name") == "prior" and isinstance(s.get("genes"), str):
            s.setdefault("label", selector_label(s))  # keeps the list name once the path is replaced by the genes
            s["genes"] = load_gene_list(s["genes"])
    lc = LosoConfig(**b)
    return cs, lc, cfg


def _cmd_run(a):
    cs, lc, cfg = load_config(a.config)
    out = Path(a.out or cfg.get("out", "results"))
    out.mkdir(parents=True, exist_ok=True)
    print(cs.summary().to_string(index=False))
    if not a.within_only:
        if a.only_harmonisers:
            lc.outer_harmonisers = a.only_harmonisers.split(",")
        results, preds, selected = run_loso(cs, lc)
        selected = {"|".join(k): v for k, v in selected.items()}
        if a.only_harmonisers:
            results, preds, selected = _merge_rerun(out, lc.outer_harmonisers, results, preds, selected, order=list(lc.harmonisers) + ["nested"])
        results.to_csv(out / "loso_results.csv", index=False)
        preds.to_csv(out / "loso_predictions.csv", index=False)
        (out / "selected_features.json").write_text(json.dumps(selected, indent=1))
        summ = summarise(results)
        summ.to_csv(out / "loso_summary.csv")
        print("\nLOSO summary (mean AUROC across held-out studies):")
        print(summ.head(20).to_string())
    if a.within or a.within_only:
        w, folds = run_within_study_cv(cs, lc, n_repeats=a.n_repeats, return_folds=True)
        w.to_csv(out / "within_study_cv.csv", index=False)
        folds.to_csv(out / "within_study_cv_folds.csv", index=False)
        print(f"\nWithin-study CV AUROC, mean over folds of 5-fold CV x {a.n_repeats} repeats (optimism reference):")
        print(w.groupby(["selector", "classifier"])[["auroc", "auroc_pooled"]].mean().to_string())


def _merge_rerun(out: Path, harmonisers, results: pd.DataFrame, preds: pd.DataFrame, selected: dict, order=None):
    """Replace the rows of ``harmonisers`` in the existing results of ``out`` with a rerun,
    keeping every other row, in the order of the existing files. A harmoniser new to the
    files is placed by ``order`` (the harmoniser loop order of run_loso), else last."""
    old_r = pd.read_csv(out / "loso_results.csv")
    old_p = pd.read_csv(out / "loso_predictions.csv")
    sel_path = out / "selected_features.json"
    old_s = json.loads(sel_path.read_text()) if sel_path.exists() else {}
    keep = {k: v for k, v in old_s.items() if k.split("|")[1] not in set(harmonisers)}
    r = pd.concat([old_r[~old_r["harmoniser"].isin(harmonisers)], results], ignore_index=True)
    p = pd.concat([old_p[~old_p["harmoniser"].isin(harmonisers)], preds], ignore_index=True)
    cols = ["held_out", "selector", "harmoniser", "classifier"]  # the loop order of run_loso
    levels = {c: list(pd.unique(pd.concat([old_r[c], results[c]]))) for c in cols}
    if order is not None:
        levels["harmoniser"] = [h for h in order if h in set(levels["harmoniser"])] + [h for h in levels["harmoniser"] if h not in set(order)]
    keys = {c: r[c].map({v: i for i, v in enumerate(levels[c])}) for c in cols}
    r = r.iloc[pd.DataFrame(keys).reset_index(drop=True).sort_values(cols, kind="stable").index].reset_index(drop=True)
    return r, p, {**keep, **selected}


def _cmd_harmonisation(a):
    from .eval.harmonisation import diagnostic_genes, run_harmonisation_diagnostics

    cs, lc, cfg = load_config(a.config, extra_features=diagnostic_genes())
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    run_harmonisation_diagnostics(cs, lc, cfg["filtered_features"], extra=diagnostic_genes(), out=out, n_pcs=a.n_pcs, k=a.k)


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
    r.add_argument("--within-only", action="store_true", help="run within-study CV only, without the leave-one-study-out benchmark")
    r.add_argument("--n-repeats", type=int, default=10, help="repeats of the stratified 5-fold within-study CV")
    r.add_argument("--out", help="output folder (default: out in the config)")
    r.add_argument("--only-harmonisers", help="comma separated rows to rerun (nested still chooses among all harmonisers); they replace the same rows of the existing results in the output folder")
    r.set_defaults(func=_cmd_run)
    d = sub.add_parser("harmonisation", help="mixing and conservation metrics of each harmoniser on every LOSO training pool")
    d.add_argument("config")
    d.add_argument("--out", required=True)
    d.add_argument("--n-pcs", type=int, default=50)
    d.add_argument("--k", type=int, default=20, help="nearest neighbours for the mixing fraction")
    d.set_defaults(func=_cmd_harmonisation)
    m = sub.add_parser("summarise")
    m.add_argument("results")
    m.set_defaults(func=_cmd_summarise)
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
