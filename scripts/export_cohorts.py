"""Export the analysis-ready LC versus COVID_recovered cohorts with a manifest of their roles and flags.

For each layer the cohorts and preprocessing of one config are used, with every cohort on its own features (the
concordance configs, ``features: per_cohort`` or ``union``):

  rna           configs/concordance_rna_v5.yaml: the eight lc_transcriptomics_blood_v4 cohorts, CORSIP2026, GSE228320;
                log2 CPM (library size over all genes) on the genes that pass the cohort's own expression filter
                (GSE297499: log2(TPM + 1))
  proteomics    configs/concordance_proteomics_all_v5.yaml: the version 5 protein cohorts (affinity and MS, with
                Woodruff2023 and Bansal2025); values as built, missing values filled per cohort as in the config
  metabolomics  configs/concordance_metabolomics_all.yaml

Samples: those the config keeps (post acute LC and COVID_recovered, one per subject). Written to --out (default
data/exports, not committed): <layer>/<cohort>_matrix.csv.gz (features x samples) and <layer>/<cohort>_meta.csv (the
ClinicalSchema metadata), and manifest.csv. The manifest is also written to exports/manifest.csv (committed).

Manifest columns: cohort, layer, role, flags, flag_notes, n_samples, n_LC, n_COVID_recovered, n_features,
feature_space, values, platform, sample_matrix, accession, overlaps_with_study, source_config, matrix_file, meta_file,
matrix_sha256, meta_sha256.

Roles:
  train_test   used for training and held out in the leave one study out benchmarks
  train_only   used for training, never held out (PRJNA1149196_QLD, 5 LC / 4 COVID_recovered; configs/lc_transcriptomics_blood_v4.yaml)
  test_only    never used for training: CORSIP2026 (NanoString panel; panel_transfer_CORSIP2026) and GSE228320 (case =
               DLCO below 60% of predicted, case_definition pulmonary_sequelae; sensitivity only)
  concordance_only   in no leave one study out benchmark, only in the concordance analyses: Saito2024 (shares 15 to 28
               metabolites with the other cohorts; configs/lc_metabolomics.yaml)
Flags:
  label_confounded_with_collection   GSE270045: LC and COVID_recovered were collected separately, so group is completely
                                     confounded with sample collection (docs/gse270045_check.md)

  python scripts/export_cohorts.py
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from chorus.cli import load_config

ROOT = Path(__file__).resolve().parents[1]
LAYERS = {
    "rna": "configs/concordance_rna_v5.yaml",
    "proteomics": "configs/concordance_proteomics_all_v5.yaml",
    "metabolomics": "configs/concordance_metabolomics_all.yaml",
}
ROLES = {"PRJNA1149196_QLD": "train_only", "CORSIP2026": "test_only", "GSE228320": "test_only", "Saito2024": "concordance_only"}  # default train_test
FLAGS = {
    "GSE270045": {
        "label_confounded_with_collection": "LC and COVID_recovered samples were collected separately; group is completely confounded with sample collection (docs/gse270045_check.md)",
    },
}
VALUES = {"log2cpm": "log2 CPM, library size over all genes", "as_provided": "as built"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def values_label(c) -> str:
    if c.study_id == "GSE297499":
        return "log2(TPM + 1), author TPM"
    v = VALUES.get(c.source.get("transform", "as_provided"), c.source.get("transform", ""))
    if c.source.get("impute"):
        v += f"; missing values filled per cohort ({c.source['impute']}, {c.source.get('n_values_imputed', 0)} values)"
    return v


def export(out: Path) -> pd.DataFrame:
    rows = []
    for layer, config in LAYERS.items():
        cs, lc, cfg = load_config(str(ROOT / config))
        (out / layer).mkdir(parents=True, exist_ok=True)
        for c in cs.cohorts:
            X = c.X.loc[:, c.X.notna().any(axis=0)]  # the cohort's own features
            meta = c.meta.loc[X.index]
            mf, tf = out / layer / f"{c.study_id}_matrix.csv.gz", out / layer / f"{c.study_id}_meta.csv"
            X.T.rename_axis(c.feature_space).to_csv(mf, float_format="%.6g")
            meta.to_csv(tf, index=False)
            st = meta["lc_status"].astype(str)
            flags = FLAGS.get(c.study_id, {})
            ov = meta["overlaps_with_study"].dropna().astype(str).unique() if "overlaps_with_study" in meta else []
            rows.append({
                "cohort": c.study_id, "layer": layer, "role": ROLES.get(c.study_id, "train_test"),
                "flags": ";".join(flags), "flag_notes": " | ".join(flags.values()),
                "n_samples": len(meta), "n_LC": int((st == lc.positive).sum()), "n_COVID_recovered": int((st == lc.negative).sum()),
                "n_features": X.shape[1], "feature_space": c.feature_space, "values": values_label(c),
                "platform": "; ".join(meta["platform"].dropna().astype(str).unique()) if "platform" in meta else "",
                "sample_matrix": ";".join(sorted(meta["sample_matrix"].dropna().astype(str).unique())),
                "accession": c.source.get("accession", ""), "overlaps_with_study": "; ".join(ov),
                "source_config": config, "matrix_file": str(mf.relative_to(out)), "meta_file": str(tf.relative_to(out)),
                "matrix_sha256": sha256(mf), "meta_sha256": sha256(tf),
            })
            print(f"{layer}/{c.study_id}: {len(meta)} samples x {X.shape[1]} features, role {rows[-1]['role']}" + (f", flags {rows[-1]['flags']}" if flags else ""), flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "data" / "exports"))
    a = ap.parse_args()
    out = Path(a.out)
    M = export(out)
    M.to_csv(out / "manifest.csv", index=False)
    (ROOT / "exports").mkdir(exist_ok=True)
    M.to_csv(ROOT / "exports" / "manifest.csv", index=False)
    print(M.groupby(["layer", "role"]).size().to_string())


if __name__ == "__main__":
    main()
