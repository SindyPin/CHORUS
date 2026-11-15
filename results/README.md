# Results

This file records how each result folder was produced. It does not interpret the results.

## Environment

- Every run below used the arm64 virtual environment `~/.venvs/chorus-arm64` (Homebrew Python 3.12.16), called as
  `arch -arm64 ~/.venvs/chorus-arm64/bin/...`.
- Package versions are pinned in `docs/requirements-chorus-arm64.txt`. Key versions: numpy 2.5.4, pandas 3.1.1,
  scikit-learn 1.9.1, torch 2.15.0, tabpfn 9.2.0, inmoose 0.9.1, xgboost 3.4.1.
- `docs/provenance.md` describes the environments and the GSE224615 re-quantification.

## TabPFN status

The TabPFN licence token file `~/.config/chorus/tabpfn_token` was absent when the real-data runs started. They ran with
`CHORUS_DISABLE_TABPFN=1`, so the `taco` classifier is the soft voting ensemble of **rf, extratrees, svc and mlp only,
without the TabPFN base learner**. The simulation runs do not use `taco`.

## Runs

Runtimes are wall-clock seconds from `results/runtimes.tsv`, all on the same laptop. The runs of 11 and 12 November were
made one after another. On 13 November `lc_transcriptomics_blood_within_v2` ran at the same time as the two morning
control runs (the label permutation and the sample_matrix control removed since), so those runtimes overlap. The later
13 November runs (from `infer_sex` on) ran one after another. `infer_sex` was run twice with identical output; the
runtime is that of the second run. The version 2 reruns of 13 November (from 11:19) ran one after another; see "Version 2
changes" for what they replaced. The version 3 reruns of 13 November (from 15:12) ran one after another; see "Version 3
changes". The runs of 15 November (from 06:49) ran one after another; see "Proteomics (15 November 2026)". The two
Test 2 runs (from 07:43) ran one after another; see "Test 2: layers on the same people (15 November 2026)". The Test 3
run (17:24) ran alone; see "Test 3: mosaic programmes (15 November 2026)". The metabolomics runs of 16 November ran one
after another (the all-cohort concordance was rerun into its final folder at 00:19:51); see "Metabolomics (16 November
2026)".

| folder | command | config | runtime (s) |
|---|---|---|---|
| `lc_transcriptomics_blood_v1` | `chorus run configs/lc_transcriptomics_blood.yaml --within` (written to `lc_transcriptomics_blood`, copied) | configs/lc_transcriptomics_blood.yaml | 1053 (started 2026-11-12 09:19:36) |
| `lc_transcriptomics_blood` | the v1 run, then `chorus run configs/lc_transcriptomics_blood.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_transcriptomics_blood.yaml | 397 for the version 2 rerun (started 2026-11-13 11:37:23), 404 for the version 3 rerun (started 15:12:49) |
| `lc_transcriptomics_blood_noGSE270045_v1` | `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --within` (copied) | configs/lc_transcriptomics_blood_noGSE270045.yaml | 805 (started 2026-11-12 09:37:09) |
| `lc_transcriptomics_blood_noGSE270045` | the v1 run, then `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_transcriptomics_blood_noGSE270045.yaml | 300 for the version 2 rerun (started 2026-11-13 11:44:00), 312 for the version 3 rerun (started 15:19:33) |
| `lc_acute_transfer_v1` | `chorus run configs/lc_acute_transfer.yaml --within` (copied) | configs/lc_acute_transfer.yaml | 231 (started 2026-11-12 09:50:34) |
| `lc_acute_transfer` | the v1 run, then `chorus run configs/lc_acute_transfer.yaml --only-harmonisers combat_addon,chorus,nested` (version 2), then `--only-harmonisers combat_addon,chorus_protected,nested` (version 3) | configs/lc_acute_transfer.yaml | 101 for the version 2 rerun (started 2026-11-13 11:49:00), 105 for the version 3 rerun (started 15:24:45) |
| `simulation_factors_v1` | `scripts/simulation_benchmark.py --out results/simulation_factors` (folder renamed) | built into the script | 231 (started 2026-11-11 09:09:10) |
| `simulation_factors` | `scripts/simulation_benchmark.py --out results/simulation_factors` | built into the script | 475 (version 3, started 2026-11-13 15:30:07) |
| `simulation_shifts_only_v1` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` (folder renamed) | built into the script | 346 (started 2026-11-11 09:13:14) |
| `simulation_shifts_only` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` | built into the script | 466 (version 3, started 2026-11-13 15:38:02) |
| `lc_transcriptomics_blood_within_v2` | `chorus run configs/lc_transcriptomics_blood.yaml --within-only --out results/lc_transcriptomics_blood_within_v2` | configs/lc_transcriptomics_blood.yaml | 1320 (started 2026-11-13 07:05:33) |
| `control_label_permutation` | `scripts/control_runs.py permute configs/lc_transcriptomics_blood.yaml --out results/control_label_permutation --n-perm 20` | configs/lc_transcriptomics_blood.yaml, changed by the script | 178 (version 3, started 2026-11-13 15:46:30) |
| `control_sex_label/inferred_sex.csv` | `scripts/infer_sex.py configs/lc_transcriptomics_blood.yaml --copy results/control_sex_label/inferred_sex.csv` | configs/lc_transcriptomics_blood.yaml (raw matrices) | 16 (started 2026-11-13 09:08:07) |
| `control_sex_label_v1` | `scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label` (folder renamed) | configs/lc_transcriptomics_blood.yaml, changed by the script | 84 (started 2026-11-13 08:56:20) |
| `control_sex_label` | `scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label` | configs/lc_transcriptomics_blood.yaml, changed by the script | 112 (rerun with `forced_features`, started 2026-11-15 06:57:24; version 3: 113) |
| `lc_transcriptomics_blood_harmonisation_v1` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` (folder renamed) | configs/lc_transcriptomics_blood.yaml | 100 (started 2026-11-13 09:01:09) |
| `lc_transcriptomics_blood_harmonisation` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` | configs/lc_transcriptomics_blood.yaml | 104 (version 3, started 2026-11-13 15:26:30) |
| `concordance` | `scripts/concordance.py configs/lc_transcriptomics_blood.yaml --out results/concordance` | configs/lc_transcriptomics_blood.yaml | 48 (started 2026-11-13 09:03:45) |
| `lc_proteomics_affinity` | `chorus run configs/lc_proteomics_affinity.yaml --within` | configs/lc_proteomics_affinity.yaml | 236 (started 2026-11-15 06:49:03) |
| `lc_proteomics_ms` | `chorus run configs/lc_proteomics_ms.yaml --within` | configs/lc_proteomics_ms.yaml | 226 (started 2026-11-15 06:52:59) |
| `control_acute_proteomics` | `chorus run configs/control_acute_proteomics.yaml --within` | configs/control_acute_proteomics.yaml | 23 (started 2026-11-15 06:56:45) |
| `lc_proteomics_affinity/harmonisation` | `chorus harmonisation configs/lc_proteomics_affinity.yaml --out results/lc_proteomics_affinity/harmonisation` | configs/lc_proteomics_affinity.yaml | 4 (started 2026-11-15 06:57:08) |
| `lc_proteomics_ms/harmonisation` | `chorus harmonisation configs/lc_proteomics_ms.yaml --out results/lc_proteomics_ms/harmonisation` | configs/lc_proteomics_ms.yaml | 1 (started 2026-11-15 06:57:12) |
| `lc_proteomics_affinity/concordance` | `scripts/concordance.py configs/lc_proteomics_affinity.yaml --out results/lc_proteomics_affinity/concordance --top 50` | configs/lc_proteomics_affinity.yaml | 2 (started 2026-11-15 06:57:13) |
| `lc_proteomics_ms/concordance` | `scripts/concordance.py configs/lc_proteomics_ms.yaml --out results/lc_proteomics_ms/concordance --top 50` | configs/lc_proteomics_ms.yaml | 1 (started 2026-11-15 06:57:16) |
| `concordance_proteomics_all` | `scripts/concordance.py configs/concordance_proteomics_all.yaml --out results/concordance_proteomics_all --top 50` | configs/concordance_proteomics_all.yaml | 7 (started 2026-11-15 06:57:17) |
| `control_acute_transcriptomics` | no run: fewer than three cohorts qualify (`cohorts.csv`) | none | none |
| `paired_Cervia2024` | `scripts/paired_layers.py --cohort Cervia2024 --out results/paired_Cervia2024` | built into the script | 2059 (started 2026-11-15 07:43:05) |
| `paired_Su2022` | `scripts/paired_layers.py --cohort Su2022 --out results/paired_Su2022` | built into the script | 1029 (started 2026-11-15 08:17:24) |
| `mosaic_programmes` | `scripts/mosaic_programmes.py --out results/mosaic_programmes` | configs/mosaic_lc_rna.yaml, configs/concordance_proteomics_all.yaml, configs/mosaic_acute_rna.yaml, configs/mosaic_acute_protein.yaml | 239 (started 2026-11-15 17:24:32; acute 28, lc 209) |
| `lc_metabolomics` | `chorus run configs/lc_metabolomics.yaml --within` | configs/lc_metabolomics.yaml | 205 (started 2026-11-16 00:16:10) |
| `lc_metabolomics/harmonisation` | `chorus harmonisation configs/lc_metabolomics.yaml --out results/lc_metabolomics/harmonisation --n-pcs 10` | configs/lc_metabolomics.yaml | 1 (started 2026-11-16 00:19:35) |
| `lc_metabolomics/concordance` | `scripts/concordance.py configs/lc_metabolomics.yaml --out results/lc_metabolomics/concordance --top 20` | configs/lc_metabolomics.yaml | 1 (started 2026-11-16 00:19:36) |
| `lc_metabolomics/concordance_all` | `scripts/concordance.py configs/concordance_metabolomics_all.yaml --out results/lc_metabolomics/concordance_all --top 50` | configs/concordance_metabolomics_all.yaml | 2 (started 2026-11-16 00:19:51) |
| `enrichment.csv` | `scripts/enrichment.py configs/lc_transcriptomics_blood.yaml results/lc_transcriptomics_blood/selected_features.json --out results/enrichment.csv` | configs/lc_transcriptomics_blood.yaml (expressed universe) | 42 (version 3, started 2026-11-13 15:45:48) |

- `lc_transcriptomics_blood`: primary leave-one-study-out benchmark, full grid (9 harmonisers, among them `chorus_protected` since version 3, plus nested selection,
  6 selector configurations, 4 classifiers), with the within-study cross-validation reference.
- `lc_transcriptomics_blood_noGSE270045`: identical except GSE270045 is disabled (sensitivity run;
  `docs/gse270045_check.md`).
- `lc_acute_transfer`: post-acute cohorts as training only, GSE215865 (acute admissions) as the single held-out study
  (`docs/acute_transfer.md`).
- `simulation_factors`: `scripts/simulation_benchmark.py` with its defaults (6 simulated studies, 3 sample-level
  technical factors).
- `simulation_shifts_only`: the same with `--n-factors 0` (study location and scale shifts only).
- `lc_transcriptomics_blood_within_v2`: within-study cross-validation only (`--within-only`) for the primary config, with
  the revised method: stratified 5-fold CV repeated 10 times, metrics per fold averaged over the 50 folds, and the pooled
  value kept for comparison. The leave-one-study-out benchmark was not rerun.
- `control_label_permutation`: label permutation control, 20 permutations. See "Control runs" below.
- `control_sex_label`: positive control with inferred sex as the label, on all filtered genes, without chrX and chrY
  genes, and with the sex genes forced in. It replaces the `sample_matrix` control (`control_matrix_label`, removed; it is in git history at commit
  aad458e). See "Inferred sex" and "Control runs" below.
- `lc_transcriptomics_blood_harmonisation`: mixing and conservation metrics of every harmoniser on every LOSO training
  pool of the primary config. See "Harmonisation diagnostics" below.
- `concordance`: within-cohort LC effects and their pairwise Spearman correlation across cohorts. See "Concordance" below.
- `enrichment.csv`: enrichment of each fold's invariant genes (from `lc_transcriptomics_blood/selected_features.json`) in
  the prior gene lists. See "Enrichment" below. Rerun in version 3 on the current `selected_features.json`.
- Folders ending in `_v1`: the results before the version 2 changes, kept for comparison. See "Version 2 changes" below.
  The version 2 results replaced in version 3 were not copied to new folders; they are in git at commit 44dffee.

The cohort table each real-data run loaded (samples and labels per study) is printed at the top of its log,
`results/<folder>.log`. The logs stay local (not committed). Simulation runs write it to `cohorts.csv`;
`control_sex_label/cohorts.csv` gives the sex counts per cohort.

**Two versions of within_study_cv.csv.** The `within_study_cv.csv` files in `lc_transcriptomics_blood`,
`lc_transcriptomics_blood_noGSE270045`, `lc_acute_transfer`, the simulation and smoke folders come from the earlier method:
one stratified 5-fold split with the out-of-fold predictions pooled into one AUROC per study. Their `auroc` equals
`auroc_pooled_repeat0` in `lc_transcriptomics_blood_within_v2` (checked: largest difference 0).

`smoke_transcriptomics` and `smoke_acute_transfer` are smoke tests from 12 November 2026 on the same cohorts with a reduced
grid: harmonisers identity and center_scale, selectors variance, differential and invariant, classifier logreg, no nested
selection, 200 bootstrap resamples.

## Cohorts in the real-data runs

| run | held out | training only |
|---|---|---|
| lc_transcriptomics_blood | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 | none |
| lc_transcriptomics_blood_noGSE270045 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849 | none |
| lc_acute_transfer | GSE215865 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 |
| control_label_permutation | as lc_transcriptomics_blood | none |
| control_sex_label (and _v1) | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE270045 (GSE251849 skipped: 4 males, below 5 per class) | none |
| lc_transcriptomics_blood_harmonisation | as lc_transcriptomics_blood (training pools of its folds) | none |
| lc_proteomics_affinity (and its harmonisation and concordance) | Buggert2025, Cervia2024_SomaScan, Su2022 | none |
| lc_proteomics_ms (and its harmonisation and concordance) | Cervia2024_MS, Chowdhury2024, Kovarik2023_plasma | none |
| control_acute_proteomics | Cervia2024_SomaScan, Su2022, Zoodsma2022_Hannover | none |
| concordance_proteomics_all | no folds; the six cohorts of the two blocks | none |

## GSE224615 samples flagged at QC

GSE224615 was re-quantified with salmon (`docs/provenance.md`). A sample is flagged when its mapping rate is below 40% or
more than 2.5 MAD (raw) below the cohort median: median 81.57%, MAD 1.77, lower bound 77.13%.

| sample | mapping rate (%) | mitochondrial fraction | reason |
|---|---|---|---|
| GSM7027478 | 52.53 | 0.580 | more than 2.5 MAD below median |
| GSM7027490 | 59.07 | 0.457 | more than 2.5 MAD below median |
| GSM7027501 | 56.84 | 0.052 | more than 2.5 MAD below median |

**All three are included in every run above.** Excluding them is left for a later sensitivity run. The full per-sample QC
table is `data/quant/GSE224615/salmon_qc.tsv` (local, not committed).

## Files in each folder

Only the CSV files, this README and `runtimes.tsv` are committed. `selected_features.json`, the figures and the logs stay
local.

- `loso_results.csv`: one row per held-out study, harmoniser, selector and classifier.
- `loso_summary.csv`: AUROC across held-out studies per harmoniser, selector and classifier.
- `loso_predictions.csv`: per-sample predicted probabilities.
- `within_study_cv.csv`: 5-fold stratified cross-validation inside each study (written with `--within`).
- `selected_features.json`: genes selected per held-out study, harmoniser and selector (real-data runs).
- `within_study_cv_folds.csv`: one row per study, selector, classifier, repeat and fold (`lc_transcriptomics_blood_within_v2`).
- `permutation_results.csv`, `permutation_summary.csv`: `control_label_permutation` only.
- `cohorts.csv`: simulation runs and `control_sex_label`. `figure_simulation.png`: simulation runs only.
- `control_sex_label`: `inferred_sex.csv`, `cohorts.csv`, `feature_sets.csv`, `forced_genes_by_cohort.csv`,
  `forced_genes_by_fold.csv`, `auroc_by_held_out.csv`, and the loso files in `all_genes/`, `no_xy/` and
  `forced_sex_genes/`. `control_sex_label_v1` has no forced set and no forced gene tables.
- `lc_transcriptomics_blood_harmonisation`: `<held_out>/<harmoniser>.csv` (one table per fold and harmoniser) and
  `all_metrics.csv` (all of them stacked).
- `concordance`: `effects.csv`, `spearman_all_genes.csv`, `spearman_all_genes_p.csv`, `spearman_top500.csv`,
  `spearman_top500_p.csv`, `top500_genes.csv`, `cohort_annotation.csv`.
- Proteomics folders (`lc_proteomics_affinity`, `lc_proteomics_ms`, `control_acute_proteomics`): the loso and
  within-study files above, plus `cohorts.csv` (samples, class counts and `n_values_imputed` per cohort) and
  `features_observed_fraction.csv` (per protein and cohort, the fraction of the cohort's samples with a value before
  imputation; one row per protein used, the proteins shared by all cohorts of the config).
  `lc_proteomics_affinity` and `lc_proteomics_ms` also hold `harmonisation/` (as `lc_transcriptomics_blood_harmonisation`)
  and `concordance/` (as `concordance`, with `top50` in place of `top500`, and `n_shared.csv`).
- `concordance_proteomics_all`: as `concordance`, with `top50` files and `n_shared.csv` (proteins with an effect in both
  cohorts of each pair).
- `control_acute_transcriptomics`: `cohorts.csv` only (see "Proteomics (15 November 2026)").
- `paired_<cohort>/<contrast>/` (`lc`, `acute`): `samples.csv`, `features.csv`, `fold_metrics.csv`, `summary.csv`,
  `predictions.csv`, `probability_correlation.csv`, and for Cervia2024 `effects.csv` and `effect_correlation.csv`;
  `paired_<cohort>/run_info.json` gives the class counts and why a contrast was skipped (see "Test 2").
- `mosaic_programmes/`: `genesets.csv`, `run_info.json` and, per contrast (`acute`, `lc`), `cohorts.csv`,
  `coverage.csv`, `effects.csv`, `effects_var.csv`, `concordance_<set>.csv` with `_p` and `n_shared_<set>.csv`
  (`all`, `hallmark`, `reactome`), `partition.csv`, `loso.csv`, `loso_predictions.csv`, and for `lc`
  `bridge_cervia.csv` (see "Test 3"). `scores.csv.gz` is not committed.
- `lc_metabolomics/`: as `lc_proteomics_affinity` (with `harmonisation/` and `concordance/` with `top20` files), plus
  `concordance_all/` (four cohorts, pairwise shared, `top50` files) and `mapping/` (`mapping_summary.csv`,
  `n_shared.csv`, `shared_by_k.csv`); see "Metabolomics".
- `concordance_epigenomics/`: `cohorts.csv` only (groups and why no effect was computed; see "Metabolomics").

## Column glossary

**loso_results.csv**

| column | meaning |
|---|---|
| held_out | study left out and used as the test set |
| harmoniser | harmoniser in the grid; `nested` means it was chosen by an inner leave-one-study-out over the training studies |
| harmoniser_chosen | harmoniser actually applied (differs from `harmoniser` only for `nested`) |
| leak_risk | true when the harmoniser is transductive (fitted on all studies, without labels) rather than on training studies only |
| n_factors | number of technical factors the `chorus` harmoniser removed in this fit (also for `nested` rows that chose it); empty for other harmonisers and in the `_v1` folders |
| selector | selector configuration: `variance`, `differential`, `invariant`, or `prior:<gene list>:<mode>[:k<k>]` |
| n_features | number of model inputs: selected genes plus covariate columns |
| n_features_dropped_incomplete | features left out of the fold because a training or held-out sample had no value (only `control_sex_label`, set `forced_sex_genes`; 0 otherwise) |
| n_strict_invariant | genes that passed the strict invariance screen (invariant selector; empty otherwise) |
| sign_agreement | over the selected genes, the fraction whose Hedges g between LC and COVID_recovered in the held-out study has the same sign as the pooled training effect; uses the held-out labels and is a diagnostic only. Empty when the selector has no effect table (variance, prior_only) or the held-out study has fewer than 3 samples per class |
| effect_spearman | Spearman correlation between the held-out Hedges g and the pooled training effect over the same genes; empty under the same conditions |
| classifier | classifier |
| n_train | training samples |
| n_train_studies | training studies |
| auroc, auprc | area under the ROC and precision-recall curves on the held-out study |
| balanced_accuracy, f1, mcc | at a probability threshold of 0.5; mcc is 0 when only one class is predicted |
| brier | Brier score |
| ece | expected calibration error |
| n | held-out samples |
| prevalence | fraction of LC among held-out samples |
| auroc_lo, auroc_hi | 95% percentile interval of the AUROC from a class-stratified bootstrap over held-out samples (`n_boot` resamples) |
| inner_scores | for `nested`: mean inner AUROC per candidate harmoniser |
| seconds | time for this harmoniser and selector fit plus this classifier |
| error | present only on rows where a fit failed; the message |

**loso_summary.csv**: per harmoniser, selector and classifier, the `mean`, `median`, `std`, `min` and `max` of AUROC
across held-out studies, `n_studies`, and the mean AUPRC (`auprc_mean`) and Brier score (`brier_mean`). Failed rows are
excluded.

**loso_predictions.csv**: `sample`, `held_out`, `harmoniser`, `selector`, `classifier`, true label `y` (1 = LC) and
predicted probability `p`.

**within_study_cv.csv** (`lc_transcriptomics_blood_within_v2`):

| column | meaning |
|---|---|
| study, selector, classifier | as above |
| auroc, auprc, balanced_accuracy, f1, mcc, brier, ece | computed on each test fold, then averaged over the `n_folds` folds |
| auroc_sd | standard deviation of the fold AUROCs |
| auroc_pooled | out-of-fold predictions pooled across the 5 folds of a repeat into one AUROC, averaged over repeats |
| auroc_pooled_repeat0 | the pooled AUROC of the first repeat, which uses the split of the earlier method |
| n, prevalence | samples in the study and fraction LC |
| n_folds, n_failed_folds | folds with metrics, and folds whose fit failed |
| n_splits, n_repeats | 5 and 10. Repeat r uses split seed `random_state + r` |

`within_study_cv_folds.csv` has the per-fold values with `repeat` and `fold`; there, `n` and `prevalence` refer to the
test fold.

**within_study_cv.csv** (earlier method, other folders): `study`, `selector`, `classifier` and the same metrics as in
loso_results.csv, from one pooled 5-fold stratified cross-validation inside the study. The invariant selector falls back to differential, because one study gives one
environment. `cross_study` is skipped. Studies with fewer than 5 samples in a class are skipped. Training-only studies
are included. This applies to both versions.

**permutation_results.csv**: the loso_results.csv columns plus `permutation` (0 to 19).

**permutation_summary.csv**: per `held_out`, `harmoniser`, `selector` and `classifier`, the `mean`, `std`, `min` and `max`
of AUROC over permutations, and `n_permutations`.

## Control runs

Both use `scripts/control_runs.py` and start from the cohorts, preprocessing and covariates of
`configs/lc_transcriptomics_blood.yaml`. TabPFN is not involved, because neither control uses `taco`.

**control_sex_label (positive control).**
- **Label:** `inferred_sex` (female positive, male negative; see "Inferred sex"), on the samples of the primary contrast
  (LC or COVID_recovered) after the primary preprocessing. Ambiguous samples are left out (one, GSM7069947).
- **Counts:** `cohorts.csv`. GSE251849 has 4 males in the contrast and is skipped as a held-out study
  (`min_test_per_class` 5); it is still used for training.
- **Grid:** harmonisers identity, center_scale, combat_addon, chorus and chorus_protected (since version 3); selectors variance (k 500) and invariant (k 200);
  classifier logreg; no nested selection; 1000 bootstrap resamples. Covariates and `protect` are `sample_matrix` as in
  the primary config.
- **Three feature sets** (`feature_sets.csv`; `_v1` has the first two):
  - `all_genes`: the 10,442 genes that pass the primary expression filter. These include 320 chrX genes and no chrY
    gene; XIST and RPS4Y1 are not among them, because the filter requires CPM of at least 1 in half the samples of every
    cohort.
  - `no_xy`: the same without chrX and chrY genes (10,122 genes). Chromosomes come from the GENCODE v47 tx2gene table
    (`data/quant/GSE224615/tx2gene_gencode_v47.tsv`); genes on both chrX and chrY (pseudoautosomal) count as both.
    Every filtered gene had a chromosome.
  - `forced_sex_genes` (version 2): `all_genes` plus XIST, RPS4Y1 and the chrY genes (chrY only, not
    pseudoautosomal) present in at least one cohort's matrix, 681 candidates in all. In each cohort a candidate is
    `used` when it has a value for every sample and, for chrY genes other than RPS4Y1, when its log2 CPM (prior count 2)
    is at least log2(1) in at least half of the cohort's inferred male samples (the config's expression thresholds
    applied to males); otherwise it is set to missing in that cohort (`forced_genes_by_cohort.csv`: `in_matrix`,
    `frac_males_expressed`, `n_males`, `used`). In each fold a candidate is a feature only when it is used in every
    training cohort and in the held-out cohort, because the held-out samples must be scored
    (`forced_genes_by_fold.csv`: `in_every_training_cohort`, `in_held_out`, `used`). The other candidates are dropped
    in that fold (`n_features_dropped_incomplete`, 673 in every fold).
  - In every fold the 8 forced genes used are RPS4Y1, KDM5D, DDX3Y, ZFY, USP9Y, TTTY14, UTY and EIF1AY. XIST is used in
    no fold: GSE270045 has no XIST row and is in every fold, as a training or the held-out cohort. In the GSE270045
    fold, XIST, TMSB4Y and LINC00278 are used in every training cohort but not in GSE270045.
  - In `selected_features.json` of `forced_sex_genes`, all 8 forced genes are among the selected genes for the variance
    selector with identity, center_scale and combat_addon in every fold, none for chorus with the variance selector,
    and none for the invariant selector. In the GSE224615 fold (identity), the 8 genes have I squared between 92 and 97
    and Q_p 0, so they fail the strict invariance screen (I squared at most 50), and 254 genes pass it; the invariant
    rows of `forced_sex_genes` therefore equal those of `all_genes`. This paragraph describes the version 3 files; since
    the 15 November rerun the forced genes are appended to every selector's output ("Proteomics (15 November 2026)", item 6).
- **auroc_by_held_out.csv**: `feature_set`, `held_out`, `harmoniser`, `selector`, `classifier`, `n_features`, `auroc`,
  `auroc_lo`, `auroc_hi`, `error`. The full rows are in `all_genes/loso_results.csv` and `no_xy/loso_results.csv`.

**control_label_permutation.**
- **Permutation:** within each cohort, `lc_status` is permuted among the samples labelled LC or COVID_recovered, so class
  counts per cohort are unchanged.
- **Benchmark:** the primary leave-one-study-out benchmark with harmonisers identity and chorus, the invariant selector
  (k 200) and logreg. Covariates and `protect` are `sample_matrix` as in the primary config. There is no nested selection.
- **Repeats:** 20 permutations (seed 0), 1000 bootstrap resamples per fit (the interval columns are not summarised).

## Inferred sex

`scripts/infer_sex.py` (module `chorus.sex`) infers sex inside each cohort of the primary config, on all samples of the
raw matrix (before one sample per subject and before the expression filter).

- **Input:** log2 CPM of XIST (ENSG00000229807) and RPS4Y1 (ENSG00000129824), library size over all genes, prior count 2.
- **Split:** a two-component Gaussian mixture (full covariance, `reg_covar` 0.25, 5 initialisations, seed 0). The
  component with the higher mean XIST minus RPS4Y1 is female.
- **Ambiguous:** posterior of the assigned component below 0.95 (`low_posterior`), or squared Mahalanobis distance to it
  above the chi-square 0.999 quantile (`outlier`).
- **GSE270045:** the author matrix has no XIST row, so the split uses RPS4Y1 alone (lower component female). The method
  string says so.
- **Stored:** columns `inferred_sex` and `inferred_sex_method` (method, thresholds, script and date) in
  `data/meta/<study_id>_meta.csv`, and the per-sample table `data/meta/inferred_sex.csv` (local), copied to
  `control_sex_label/inferred_sex.csv`. `scripts/build_cohort_meta.py` merges it back when it rebuilds the tables.

| cohort | samples | female | male | ambiguous | genes | agreement with recorded sex |
|---|---|---|---|---|---|---|
| GSE226260_WB | 103 | 67 | 36 | 0 | XIST, RPS4Y1 | no recorded sex |
| GSE226260_PBMC | 86 | 69 | 16 | 1 | XIST, RPS4Y1 | 85 of 85 non-ambiguous |
| GSE224615 | 36 | 20 | 16 | 0 | XIST, RPS4Y1 | 36 of 36 |
| GSE251849 | 23 | 18 | 5 | 0 | XIST, RPS4Y1 | no recorded sex |
| GSE270045 | 36 | 28 | 8 | 0 | RPS4Y1 only | no recorded sex |

- The ambiguous sample is GSM7069947 (GSE226260_PBMC, recorded female): XIST -0.22 and RPS4Y1 -0.36 log2 CPM, squared
  Mahalanobis distance 35.3. It is the sample with the malformed matrix tail (`docs/gse226260_check.md`).
- No subject with several samples has different inferred sex across its samples.

**inferred_sex.csv** columns: `study_id`, `sample_id`, `inferred_sex`, `posterior` (of the assigned component),
`mahalanobis2`, `reason` (empty, `low_posterior` or `outlier`), `xist_logcpm`, `rps4y1_logcpm`, `provenance`,
`recorded_sex` (the curated `sex` column; empty where not recorded).

## Harmonisation diagnostics

`chorus harmonisation` (module `chorus.eval.harmonisation`) uses the folds of the primary benchmark. For each held-out
study, the training pool is the labelled (LC or COVID_recovered) samples of the other studies. Each of the 9 harmonisers
(8 in `_v1`; `chorus_protected` added in version 3) is fitted as in the benchmark: on the training pool with its labels where the operator uses them, and for
`combat_pooled` on the pool plus the held-out samples without labels. The pool is compared before (log CPM as loaded)
and after.

- **Genes.** Mixing and LC effect metrics use the 10,442 filtered genes. XIST, RPS4Y1 and the marker genes are added to
  the harmonised matrix (none of the 10 markers needed adding; XIST and RPS4Y1 did). A gene is added in a fold only
  when every sample the harmoniser is fitted on has a value. GSE270045 has no XIST row, so XIST enters only the fold
  where GSE270045 is held out, and not for `combat_pooled` (which also fits on the held-out samples). Column
  `extra_genes` lists what was added.
- **Markers.** Neutrophil: FCGR3B, CSF3R, CXCR2, FPR1, MNDA. Lymphocyte: CD3E, CD3D, CD79A, MS4A1, CD2 (the lists in
  `chorus.harmonise.transcriptomics.composition_covariates`).
- **Sex labels** are `inferred_sex`; ambiguous samples are left out of the sex metrics.

| metric | family | definition |
|---|---|---|
| silhouette_study | technical_mixing | mean silhouette width with study as the label, Euclidean, on the top 50 principal components of the pool after each gene is standardised across the pool |
| knn_other_study | technical_mixing | mean over samples of the fraction of the 20 nearest neighbours (same space) from another study |
| knn_other_study_expected_random | technical_mixing | the same fraction under random mixing, mean of (N - n_study) / (N - 1) |
| lc_effect_spearman | biological_conservation | per training study, Spearman correlation over genes of the within-study Hedges g (LC minus COVID_recovered) before versus after; `mean_over_studies` is their mean |
| sex_auroc_XIST | biological_conservation | AUROC of XIST for female versus male over the pool |
| sex_auroc_RPS4Y1 | biological_conservation | AUROC of RPS4Y1 for male versus female over the pool |
| silhouette_sex | biological_conservation | mean silhouette width by inferred sex on XIST and RPS4Y1 (each standardised across the pool), on those of the two available |
| marker_score_spearman_neutrophil, marker_score_spearman_lymphocyte | biological_conservation | score = mean over the marker genes of the values standardised across the pool; Spearman correlation of the score before versus after, over the pool and within each study (`mean_over_studies` is the mean of the within-study values) |

**Table columns** (`<held_out>/<harmoniser>.csv` and `all_metrics.csv`): `held_out`, `harmoniser`, `family`, `metric`,
`scope` (`pool`, a study id or `mean_over_studies`), `before` and `after` (mixing and sex separation), `before_vs_after`
(LC effect and marker comparisons), `n` (samples, or studies for `mean_over_studies` of lc_effect_spearman), `n_genes`,
`n_pcs`, `k`, `leak_risk`, `n_factors` (`chorus` and `chorus_protected` only; not in `_v1`), `extra_genes`. A failed fit would give one row with `family` = `error` and the message in
`scope`; there were none.

**combat_addon with protect.** The `_v1` folder was made with the version 1 `AddonComBat`, in which a protected
covariate constant within a study moved the study location (see "Version 2 changes"). The current folder was made
with the version 3 `AddonComBat` (see "Version 3 changes").

## Concordance

`scripts/concordance.py` on the primary config, after its preprocessing (one sample per subject, deduplication).

- **Effect:** inside each cohort, Hedges g of LC minus COVID_recovered per gene on log CPM, over the 10,442 common
  expressed genes (the genes that pass the expression filter in every cohort).
- **Correlation:** pairwise Spearman correlation of g between cohorts, over all common genes and over the 500 genes
  with the largest mean |g| across the five cohorts.
- **Annotation** (`cohort_annotation.csv`), from the metadata: `case_definition` (of LC samples), `control_type`
  (`lc_status` of the negative class), `control_case_definition`, `time_anchor` (`missing` when the column is empty),
  `sample_matrix`, `n_LC`, `n_control`, `other_lc_status` (other statuses present in the cohort and not used), and
  median, minimum and maximum `days_since_infection` per class (empty when not recorded).

**Files:** `effects.csv` (`gene_id`, `symbol` from the HGNC complete set, one column of g per cohort);
`spearman_all_genes.csv` and `spearman_top500.csv` (cohort by cohort correlation matrices); `*_p.csv` (two-sided p
values of the same pairs); `top500_genes.csv` (`gene_id`, `symbol`, `mean_abs_g`).

## Enrichment

`scripts/enrichment.py` takes every `<held_out>|<harmoniser>|invariant` entry of
`results/lc_transcriptomics_blood/selected_features.json`: 5 folds by 10 harmonisers (including `nested` and
`chorus_protected`), 200 genes each.

- **Test:** one-sided hypergeometric (overlap at least as large as observed) over the expressed universe, the 10,442
  genes that pass the primary expression filter, which is the set the selector chose from. List genes outside the
  universe are dropped first.
- **Lists:**
  - `TACO_411`: `priors/taco_411_dce_genes.txt`, 411 genes, 333 in the universe.
  - `MRCONTROL_32`: `priors/mrcontrol_32_genes.txt`, 32 genes, 25 in the universe.
  - `review_complement`: `priors/review_complement_9.txt`. The file holds 11 UniProt accessions (C4A, C4BPB, C5, CFHR1,
    CFHR5, MBL2, FCN2, MASP2, C1QA, C1QB, C1QC). They map to 11 Ensembl genes through the `uniprot_ids` column of the
    HGNC complete set, and 2 are in the universe.
- **q_bh:** Benjamini-Hochberg over all 150 rows (135 before version 3).

**enrichment.csv** columns: `held_out`, `harmoniser`, `selector`, `gene_list`, `list_size` (genes after mapping),
`universe`, `list_in_universe`, `n_selected`, `overlap`, `expected` (n_selected x list_in_universe / universe),
`fold_enrichment` (overlap / expected), `p_hypergeom`, `q_bh`, `overlap_genes` (symbols).

## Version 2 changes (13 November 2026)

Three code changes, after which the affected results were rerun. The results before them are kept in the `_v1` folders.
The `_v1` results therefore differ from the current ones by both the `combat_addon` and the `chorus` change.

**1. AddonComBat (commit ee8248a).** When a protected covariate is constant within a training study, it is now dropped
from that study's design, so it no longer moves the study location. In version 1 the design had an intercept plus the
protected columns plus the label; a protected column constant within the study (as `sample_matrix` is in every cohort)
made the least-squares intercept non-unique, and the minimum-norm intercept of `numpy.linalg.lstsq` was used as the
study location. The label column is unchanged, so the location is still the intercept at label 0 when the label is in
the design. Test: two studies at means 0 and 5 with a one-hot protect by study; after the fix both study means are within
0.3 of 2.5 and equal to the result without protect (version 1 gave a difference of 6.8 between the study means in this test).

**2. Sex control (commit 564f9e5, with fixes f268a12 and c0cd0bd).** The third feature set
`forced_sex_genes` (see "Control runs"). `LosoConfig.drop_incomplete_features` drops, per fold, features with a missing
value in a training or held-out sample.

**3. chorus harmoniser (commit 34c59fb).**
- **Number of factors.** Chosen by parallel analysis on the null gene block of the training pool: squared singular
  values of the centred null block against the 0.95 quantile of those from 20 copies with each column permuted
  independently (seed 0); k is the number of leading components above it, and can be zero (no factor removed). Version 1
  used the number of training studies plus 4 (`k_method: studies_plus`, still available).
- **Protected marker scores** (option `protect_scores`, off by default and not used in the version 2 reruns; used by
  `chorus_protected` since version 3): named scores
  (`neutrophil`, `lymphocyte`, with the marker genes of the harmonisation diagnostics, or a dict of gene lists) enter
  the joint loading fit alongside the outcome and are dropped at transform time. A score is the mean of its genes
  standardised across the training pool. Its genes are left out of the null set, and the factors and the parallel
  analysis use the null block after the scores are projected out of it. Set it in a config with
  `benchmark: harmoniser_params: {chorus: {protect_scores: [neutrophil, lymphocyte]}}`.

**Reruns.**

| folder | what was rerun | kept as |
|---|---|---|
| `lc_transcriptomics_blood` | rows `combat_addon`, `chorus` and `nested` (`chorus run --only-harmonisers`); other rows, `within_study_cv.csv` and `within_study_cv_folds.csv` unchanged | `lc_transcriptomics_blood_v1` (full copy) |
| `lc_transcriptomics_blood_noGSE270045` | the same | `lc_transcriptomics_blood_noGSE270045_v1` |
| `lc_acute_transfer` | the same | `lc_acute_transfer_v1` |
| `simulation_factors`, `simulation_shifts_only` | whole run (no protect there, so only `chorus` and `nested` rows can change) | `simulation_factors_v1`, `simulation_shifts_only_v1` |
| `control_sex_label` | whole run, now with three feature sets | `control_sex_label_v1` |
| `lc_transcriptomics_blood_harmonisation` | whole run | `lc_transcriptomics_blood_harmonisation_v1` |

- `--only-harmonisers` reruns the named rows with the same folds and seeds; `nested` still chooses among all harmonisers
  of the config. The rerun rows replace the old ones in `loso_results.csv`, `loso_predictions.csv` and
  `selected_features.json`, and `loso_summary.csv` is recomputed. Rows of the other harmonisers were checked to be
  identical to `_v1` in all five rerun folders and in `control_sex_label` (identity and center_scale).
- Factors chosen by parallel analysis (`n_factors`): primary benchmark 5 to 8 per fold (version 1: 8); noGSE270045 4 to
  8 (version 1: 7); acute transfer 7 (version 1: 9); harmonisation diagnostics 5 to 8; sex control 4 to 8;
  `simulation_factors` 18 and 19 (the simulation has 3 sample-level factors and 6 studies); `simulation_shifts_only` 4.
- No fit failed in any rerun.
- **Not rerun:** `enrichment.csv` (it reads the version 1 `selected_features.json`, whose `combat_addon`, `chorus` and
  `nested` entries have since been replaced), `control_label_permutation` (uses `chorus`), `concordance` and
  `lc_transcriptomics_blood_within_v2` (no harmoniser), and the smoke runs. `enrichment.csv` and
  `control_label_permutation` were rerun in version 3.

## Version 3 changes (13 November 2026)

Two code changes, after which the affected results were rerun. The version 2 results they replace are in git at commit
44dffee. The `_v1` folders are unchanged.

**1. AddonComBat (commit 34b2eec).**
- **Diagnosis.** After version 2, the `combat_addon` training pool kept study structure in the harmonisation
  diagnostics: pool silhouette by study 0.21 to 0.26 and other-study neighbour fraction 0.07 to 0.12, against -0.03 to
  -0.02 and 0.57 to 0.63 for `center_scale`. On the fold with GSE226260_WB held out, the version 2 operator was refitted
  with parts of it switched off (scratch script, not kept). The columns are:
  - the maximum over training studies of |study mean minus pool mean|, median over genes, after fit and transform;
  - the pool silhouette by study;
  - the other-study neighbour fraction.

  | variant (version 2 code) | max abs study mean deviation | silhouette_study | knn_other_study |
  |---|---|---|---|
  | as run: protect `sample_matrix` and label in the design, `shrink` 0.5 | 1.378 | 0.229 | 0.102 |
  | label only in the design, `shrink` 0.5 | 1.378 | 0.229 | 0.102 |
  | protect and label, `shrink` 0 | 0.637 | 0.058 | 0.400 |
  | label only, `shrink` 0 | 0.637 | 0.058 | 0.400 |
  | no design, `shrink` 0.5 | 1.326 | 0.208 | 0.123 |
  | no design, `shrink` 0 | 0.000 | -0.033 | 0.605 |
  | `center_scale` | 0.000 | -0.033 | 0.605 |

  - **Standardisation:** the training studies were shifted and scaled to the grand mean and pooled variance per gene of
    the training pool.
  - **Shrinkage:** each study's per gene location deviation from the grand mean was replaced by 0.5 times itself plus
    0.5 times its average over genes. The log variance ratio was shrunk the same way. So half of each gene's study
    offset, less its across-gene average, was kept.
  - **Covariates:** `sample_matrix` is constant within every cohort and was already dropped from the design by the
    version 2 fix (the rows with and without it are identical). The label stayed in the design uncentred, so the
    location was the intercept at label 0 (the control mean), and studies with different case fractions kept different
    means.
- **Fix.**
  - The label and protected covariates are centred within each training study before the least-squares fit, so the
    intercept is the study mean. They still remove their variation from the variance estimate.
  - `shrink` defaults to 0 (it was 0.5). The option is kept.
  - The held-out study is unchanged: location and scale come from that study alone, without covariates.
  - On the same fold the fixed operator gives max deviation 0.000, silhouette -0.043 and other-study fraction 0.610.
- **Test** (`test_combat_addon_training_study_means_equal_after_transform`):
  - Setup: three studies with per gene offsets and scales, case fractions 0.2, 0.5 and 0.8, a label effect on half the
    genes, and a protected covariate that varies within studies.
  - After fit and transform, the per study gene means equal the pool means within 1e-8, with and without protect.
  - With `shrink` 0.5 they differ by more than 0.1.
  - It fails on the version 2 code (largest deviation 6.2).

**2. chorus_protected (commit bcf509d).** A second chorus variant, registered as harmoniser `chorus_protected`:
`InvarianceAnchored(protect_scores=["neutrophil", "lymphocyte"])`.
- Its other settings are those of `chorus` (factor count by parallel analysis).
- Its marker genes are those of the harmonisation diagnostics. All 10 are among the filtered genes of the three
  benchmarks.
- It was added to `harmonisers` in the three benchmark configs and in the sex control grid. The harmonisation
  diagnostics use the config list, so they include it.
- New config field `nested_candidates` (`LosoConfig.nested_candidates`): the harmonisers that `nested` chooses among,
  default `harmonisers`. The three configs set it to the version 2 list, so `chorus_protected` is not a nested
  candidate.
- It is not in the simulations. The simulated genes (`G00000` to `G01999`) contain no marker genes, so the protected
  scores cannot be formed (the operator raises an error when a named score has no genes). The simulation grid is
  unchanged.

**3. --only-harmonisers merge (commit 2b58d1b).** A harmoniser new to the results is placed in the `run_loso` loop order,
here `chorus_protected` between `chorus` and `nested`.

**Reruns.**

| folder | what was rerun |
|---|---|
| `lc_transcriptomics_blood`, `lc_transcriptomics_blood_noGSE270045`, `lc_acute_transfer` | rows `combat_addon`, `chorus_protected` (new) and `nested` (`--only-harmonisers combat_addon,chorus_protected,nested`); other rows and the within-study files unchanged |
| `simulation_factors`, `simulation_shifts_only` | whole run (only `combat_addon` and `nested` rows can change) |
| `control_sex_label` | whole run, now with `chorus_protected` |
| `lc_transcriptomics_blood_harmonisation` | whole run, now with `chorus_protected` |
| `enrichment.csv` | whole run on the current primary `selected_features.json` |
| `control_label_permutation` | whole run (identity and the current `chorus`) |

- **Unchanged rows checked against version 2:**
  - In the three benchmarks, both simulations and `control_sex_label`, the rows of the harmonisers not rerun are
    identical.
  - In the harmonisation diagnostics, every harmoniser other than `combat_addon` gives the same values.
  - In `control_label_permutation` the identity rows are identical.
- **enrichment.csv:** compared with the file computed before version 2, the `chorus` and `nested` rows changed (4 each),
  the `combat_addon` rows are unchanged, and the 15 `chorus_protected` rows are new.
- No fit failed in any rerun.
- **Factors chosen by parallel analysis for `chorus_protected` (`n_factors`):** primary benchmark 5 to 9 per fold;
  noGSE270045 4 to 10; acute transfer 5; harmonisation diagnostics 5 to 9.
- **Not rerun:** `concordance` and `lc_transcriptomics_blood_within_v2` (no harmoniser), and the smoke runs.

## Code changes after version 3 without reruns (14 November 2026)

- **Sex control, forced sex genes (commit efcf77b).**
  - **Why the rows were equal.** In `control_sex_label`, the `forced_sex_genes` rows equal the `all_genes` rows for the
    invariant selector only. The forced genes reached the candidate features, but the invariant selector dropped them
    (I squared 92 to 97 against its limit of 50), so the selected features were the same.
  - **Variance selector.** With identity, center_scale and combat_addon, the variance selector already kept all 8 forced
    genes. With chorus and chorus_protected it kept none.
  - **Fix.** The forced genes present in a fold are now appended to every selector's output (`LosoConfig.forced_features`,
    result column `n_forced_added`).
  - **Status.** `control_sex_label` was not rerun on 14 November. It was rerun on 15 November; see "Proteomics (15
    November 2026)", item 6.
- **Layer adapters (commit b5083c4).** The transcriptomics loader now goes through the `bulk_rna` adapter. The loaded
  matrices of `configs/lc_transcriptomics_blood.yaml` and `configs/lc_acute_transfer.yaml` were checked to be identical
  (hash of every cohort matrix) before and after the change.

## Proteomics (15 November 2026)

**1. Curation applied first** (in the builders; `docs/proteomics_cohorts.md`, "Curation decisions"):
- **Cervia2024:** the 5 healthy controls that match later COVID-19 subjects are linked to them
  (`linked_healthy_subject_id`), and their healthy samples are dropped (SomaScan: 39 to 34 healthy).
- **Chowdhury2024:** the authors' labels are kept, the controls have `control_type` convalescent_unverified, and the 4
  unlabelled runs are excluded.
- **Zoodsma2022:** the 16 Nijmegen patients stay excluded.
- **Captur2022:** not used.

The schema gains `control_type` convalescent_unverified and the pooled layer `proteomics`.

**Loader options added for these configs** (`chorus.cli.load_config`, `preprocessing`). Transcriptomics configs do not
use them.

| option | effect |
|---|---|
| `drop_sample_phases` | drops samples whose `sample_phase` is listed (Test 1: `acute`, because acute draws carry the later outcome in `lc_status`) |
| `keep_lc_status` | keeps only the listed `lc_status` values, before one sample per subject |
| `timepoint_order` | timepoint labels from earliest to latest, for one sample per subject. Chowdhury2024 keeps the pre-intervention sample, and Su2022 and Zoodsma keep the earliest acute draw. |
| `relabel_acute` | `lc_status` becomes `acute` for every sample with `sample_phase` acute (the curated value is kept in `lc_status_original`) |
| `impute` | per cohort and protein, from that cohort's own retained values, without labels. `study_median` for affinity (values masked at QC); `study_min` for MS (missing mostly below detection). A cohort entry can set its own `impute`. |
| `features: union` | every cohort keeps its own proteins (pairwise concordance); the default is the proteins shared by all cohorts |

- **Shared proteins.** For a cohort that is not counts, a protein counts as present when it has at least one value
  among the retained samples.
- **Other choices.** No covariates and no `protect` are used. Repeated runs of one Buggert2025 sample (20 subjects) are
  reduced to the first in file order.

**2. Test 1, LC versus COVID_recovered** (`lc_proteomics_affinity`, `lc_proteomics_ms`).
- **Samples:** one post-acute sample per subject. Leave-one-study-out within each block, on the proteins shared by all
  three cohorts of the block.
- **Grid:**
  - harmonisers identity, center_scale, rank, fsqn, combat_addon, chorus, and nested selection over these six;
  - selectors variance (k 100), differential (k 50), invariant (k 50);
  - classifiers logreg and rf;
  - 1000 bootstrap resamples;
  - within-study CV with the per-fold method (5 folds x 10 repeats, `--within`).
- **Selector sizes:** k was chosen for these blocks; the transcriptomics grid uses 500/200/200 on 10,442 genes.
- **Failed fits:** none, in 126 rows per block.

| block | cohort | samples (LC / COVID_recovered) | proteins used | values imputed |
|---|---|---|---|---|
| affinity | Buggert2025 | 197 (146 / 51) | 376 | 642 (median) |
| affinity | Cervia2024_SomaScan | 113 (40 / 73) | 376 | 0 |
| affinity | Su2022 | 125 (75 / 50) | 376 | 0 |
| MS | Cervia2024_MS | 113 (40 / 73) | 182 | 907 (minimum) |
| MS | Chowdhury2024 | 30 (21 / 9) | 182 | 56 (minimum) |
| MS | Kovarik2023_plasma | 26 (13 / 13) | 182 | 646 (minimum) |

The per-protein observed fractions are in `features_observed_fraction.csv`. In the MS block, 49 (Cervia2024_MS), 11
(Chowdhury2024) and 69 (Kovarik2023_plasma) of the 182 proteins have at least one missing value; the lowest observed
fractions are 0.44, 0.20 and 0.08.

**Su2022 and the harmonisers.** The Su2022 values are the authors' Z scores: batch corrected, adjusted for age, sex and
BMI, and scaled by the mean and standard deviation of the residuals of the matched healthy controls. The other two
affinity cohorts are NPX (Buggert2025) and log2 RFU (Cervia2024_SomaScan). What each harmoniser does to Su2022:

| harmoniser | on Su2022 | meaningful for Su2022 |
|---|---|---|
| identity | pools the Z scores with NPX and log2 RFU on their own scales | no: the units differ by construction |
| center_scale | per protein z-score within Su2022's retained samples, a linear map of the authors' values; the same operation as on the other cohorts | yes. The authors' age, sex and BMI adjustment stays, and is not applied to the other cohorts. |
| combat_addon | per protein location and scale of the study mapped to the training pool's grand mean and pooled variance | yes, as center_scale. The pool's grand mean and variance mix the three units. |
| fsqn | maps each protein's within-study ranks onto the training pool's quantiles of that protein | the within-study ranks are unaffected by the authors' per-protein transform. The reference quantiles are taken from the raw training values, so they mix units whenever cohorts with different units are in training (any affinity fold). |
| rank | within-sample rank across proteins | no: in Su2022 it orders proteins by deviation from healthy controls; in NPX and RFU cohorts, by measured level |
| chorus | centring by the training grand mean, factors from the null proteins of the pooled training values, then a per-study z-score | partly. The final per-study z-score is as center_scale. The factors are estimated on the pool before any per-study scaling, so they are fitted on mixed units. |
| nested | chooses one of the above per fold | as the harmoniser chosen (`harmoniser_chosen`) |

**3. Diagnostics and concordance.**
- **Harmonisation diagnostics:** `lc_proteomics_affinity/harmonisation` and `lc_proteomics_ms/harmonisation`, same
  metrics as for transcriptomics.
  - There are no sex rows (the cohorts have no `inferred_sex`).
  - The marker rows have `n_genes` 0 and empty values, because the marker genes are Ensembl ids, not UniProt.
  - n_pcs is 50.
  - No failed fit.
- **Concordance per block:** `lc_proteomics_affinity/concordance` and `lc_proteomics_ms/concordance`, on the 376 and
  182 shared proteins after the same preprocessing and imputation as the benchmark. `--top 50`.
  - Symbols come from the HGNC `uniprot_ids` column; a multi-protein key gets its symbols joined by "_", and an
    unmapped key has no symbol.
- **Concordance across all six cohorts:** `concordance_proteomics_all`, from `configs/concordance_proteomics_all.yaml`.
  - Same sample selection, with imputation per cohort (median for affinity, minimum for MS).
  - Each cohort is kept on its own proteins (`features: union`, 7,419 in all; 38 without a symbol), and each pair is
    correlated over the proteins both have.
  - Deduplication is off, because Cervia2024_SomaScan and Cervia2024_MS are the same people.
  - `n_shared.csv` gives the number of proteins per pair. Su2022 shares 4 or 5 proteins with each MS cohort, so those
    three correlations rest on 4 or 5 values.
  - Chowdhury2024 has 404 proteins here: one of its 405 has no value in the 30 retained samples.

**4. Positive control in proteomics, acute versus never_infected** (`control_acute_proteomics`).
- **Label:** acute = `sample_phase` acute, whatever the later outcome; never_infected = curated `lc_status`. The earliest
  acute draw per subject is kept.
- **Grid:** the six harmonisers plus nested; selectors variance (k 100) and invariant (k 50); logreg; within-study CV.
- **Cohorts:** those with both classes after one sample per subject and removal of subjects measured in two cohorts.
  - **Zoodsma2022_Breda is disabled.** Its 5 never infected samples are bridge samples of Hannover healthy donors
    already in Zoodsma2022_Hannover (TK_10, TK_33, TK_55, TK_63, TK_65), so it keeps acute samples only.
  - **Zoodsma2022_Nijmegen** has no never infected sample.
  - **Zoodsma2022_Hannover's acute samples** are 13 bridge samples (12 patients) of the Nijmegen hospitalised cohort.

| cohort | acute | never_infected | proteins used | values imputed |
|---|---|---|---|---|
| Cervia2024_SomaScan | 113 | 34 | 163 | 0 |
| Su2022 | 204 | 289 | 163 | 0 |
| Zoodsma2022_Hannover | 12 | 62 | 163 | 82 |

- **chorus failed in every fold** (6 error rows, both selectors): "only 0 null genes passed the screen" (Cervia2024_SomaScan
  held out), "only 6" (Su2022) and "only 4" (Zoodsma2022_Hannover). The default null screen of `InvarianceAnchored` was
  not changed.
- **Why the chorus rows are absent (checked 15 November).** They are present in `loso_results.csv` as error rows
  (column `error`, no metrics). They are absent from `loso_summary.csv` and `loso_predictions.csv`, which keep only
  fits that succeeded.
  - **The screen.** `InvarianceAnchored` needs at least 20 null proteins. A null protein has a meta-analytic p above
    0.2, an absolute Hedges g below 0.35 in every training cohort, I squared at most 50, and training variance at or
    above the 20th percentile.
  - **Per fold,** on the 163 proteins of the acute contrast:

    | held out | p above 0.2 | max abs g below 0.35 | I squared at most 50 | variance ok | all four | needed |
    |---|---|---|---|---|---|---|
    | Cervia2024_SomaScan | 40 | 5 | 40 | 130 | 0 | 20 |
    | Su2022 | 85 | 17 | 47 | 130 | 6 | 20 |
    | Zoodsma2022_Hannover | 73 | 14 | 29 | 130 | 4 | 20 |

  - **Decision.** The failure is the operator refusing to fit, not a bug, so it is documented and the screen
    thresholds were not relaxed. Relaxing them would make "null" proteins that are associated with the label in the
    training cohorts.
- `nested` chooses among all six candidates. A candidate that fails in an inner fold scores missing there.

**5. Positive control in transcriptomics, acute versus never_infected** (`control_acute_transcriptomics/cohorts.csv`).
- **Survey:** every transcriptomics cohort on disk was checked for acute COVID-19 samples and never infected or healthy
  samples. Only GSE215865 has both: 1,218 libraries (424 subjects) of COVID-19 positive admissions, and 151 libraries
  (64 subjects) of COVID-19 negative admissions curated as never_infected. Those are acutely ill patients without
  COVID-19, not healthy donors.
- **GSE226260 acute set (p20169):** 104 acute samples with a severity. The other 38 samples (36 subjects) have
  "pasc severity: N/A" and no other field. The paper reports 35 uninfected participants, but GEO labels no sample as
  uninfected.
- **Other cohorts:** GSE169687 (14 never infected) and GSE251849 (7) have no acute samples. GSE222253 (11 in hospital)
  has no uninfected samples.
- **Not run:** at most one cohort qualifies (two if the GSE226260 N/A samples were taken as uninfected), against the
  three required.

**6. Sex control rerun** (`control_sex_label`, with `LosoConfig.forced_features`).
- **Forced genes added:** `n_forced_added` is 8 in every `forced_sex_genes` row except the variance selector with
  identity, center_scale and combat_addon, where it is 0 because the selector already chose all 8. The forced genes
  therefore reach the model in every row.
- **Rows equal to all_genes:** 2 of the 40 `forced_sex_genes` rows have the same AUROC as `all_genes` (GSE224615 held
  out, invariant selector, center_scale and combat_addon).
- **all_genes and no_xy:** identical to the version 3 files except `seconds` and the new `n_forced_added` column (0).
- **Earlier statement superseded:** the statement under "Control runs" that the invariant rows of `forced_sex_genes`
  equal those of `all_genes` describes the version 3 files.

**Files not committed:** the logs (`results/<run>.log`, the diagnostics and concordance logs named with `_` in place of
`/`) and `selected_features.json`.

## Test 2: layers on the same people (15 November 2026)

**Paired sets** (`docs/paired_cohorts.md`). A sample is paired when its `sample_id` is in every layer and its labels
agree across the layers' metadata. One sample per subject is kept: the earliest by timepoint.

| cohort | layers (features used) | LC vs COVID_recovered | acute vs never_infected |
|---|---|---|---|
| Cervia2024 | SomaScan (6402), MS (301 of 308; 7 dropped as observed in fewer than half the samples) | 40 / 73 | not run: the MS layer has no acute or healthy sample |
| Su2022 | Olink (441), metabolomics (891, `data/metabolomics/Su2022_metabolomics`, Table S2.2) | 68 / 45 (T3) | 198 / 27 (earliest paired acute draw; healthy donors in both sheets) |
| Kovarik2023 | not built | | |

**Kovarik2023 not built.** The supplement was downloaded and checked, but its tables share only column labels
(`healthy_<n>` and the like). Those labels do not link the triglyceride table (S4) to the Biocrates table (S3) on the
126 triglycerides they share, so they cannot be relied on to link the eicosanoid and Biocrates tables to the proteomics
runs (`docs/paired_cohorts.md`).

**Models** (`chorus.eval.paired.PairedCV`). Repeated stratified 5-fold CV, 10 repeats, split seed `random_state + r`;
every model is fitted and scored on the same folds.
- **Missing values:** filled inside each fold from that fold's training samples (MS: per protein minimum; the other
  layers: median). Features observed in fewer than half of the cohort's paired samples are dropped first (label free).
- **Classifier:** logreg (the CHORUS `logreg`: standardised, C 0.5, balanced classes) for every model except DIABLO.

| model | per fold |
|---|---|
| `layer:<name>` | features standardised on the training fold, top 50 by absolute training effect (differential selector), logreg |
| `early_fusion` | the per-layer selections concatenated, each layer divided by the square root of its number of selected features |
| `late_fusion` | `chorus.integrate.mosaic.LateFusion` fitted on inner 5-fold out-of-fold layer probabilities of the training fold, applied to the layer models refitted on the whole training fold |
| `mofa` | MOFA+ (mofapy2 0.7.5), 10 factors, one group, layers as views (feature names prefixed by the layer), all features, fitted on the training fold; test samples projected onto the training loadings (`mofa_factors`); logreg on the factors |
| `diablo` | mixOmics 6.32.0 `block.splsda` in R 4.5.3, 2 components, keepX 20 per block and component, design 0.1 off the diagonal, `scale = TRUE`; the score is the positive-class column of `predict(...)$WeightedPredict` at component 2, clipped to [0, 1]. It is not a probability, so its Brier score and ECE are those of this score. |

**Agreement between layers.**
- `effect_correlation.csv` (Cervia2024 only, the one pair with a shared feature space): Spearman correlation of the
  Hedges g (LC minus COVID_recovered, all 113 paired samples) between SomaScan and MS over the UniProt keys in both
  layers (173).
- `probability_correlation.csv`: per repeat and pair of models, Spearman and Pearson correlation across subjects of
  the out-of-fold probabilities. Every subject is predicted once per repeat. The pairs include the layer models with
  each other and with the fusion models.

**Positive control:** acute versus never_infected, the same models and tables, in Su2022 only (the only paired cohort
with both classes).

**Failures:** none. Every model has 50 folds in each run (`n_failed_folds` 0).

**Column glossary.**
- `fold_metrics.csv`: `model`, `repeat`, `fold`, `n` (test samples), `prevalence`, and the metrics of
  loso_results.csv on that fold; `error` is filled on a failed fold.
- `summary.csv`: per model, the mean over folds of each metric, `auroc_sd`, `n_folds` and `n_failed_folds`.
- `predictions.csv`: `sample`, `repeat`, `fold`, `model`, `y`, `p`, `error`.
- `samples.csv`: `sample_id`, `subject_id`, `y`, `timepoint`, `label`.
- `features.csv`: `layer`, `n_features_cohort`, `n_dropped_missing`, `n_used`, `impute`.
- `effects.csv`: Hedges g per feature (rows) and layer (columns).
- `effect_correlation.csv`: `layer_a`, `layer_b`, `n_shared`, `spearman`, `p`.

**Not committed:** the logs (`results/paired_<cohort>.log`).

## Test 3: mosaic programmes (15 November 2026)

**Method:** `docs/mosaic.md`.
- **Programmes:** MSigDB 2026.2 Hallmark (50) and Reactome (1,867), CC BY 4.0, mapped by Entrez id through HGNC to Ensembl and UniProt.
- **Scores:** in each cohort, the mean of the members standardised within that cohort, with at least 5 measured members.
- **Environment:** cohort and layer.
- **Harmonisers and classifier:** identity and center_scale, with logreg.
- **Run order:** the positive control (acute versus never_infected) ran first.

**Gene set mapping** (`genesets.csv`). Hallmark: 7,322 gene memberships, 7,321 found in HGNC with an Ensembl id, 7,317 with a UniProt accession. Reactome: 105,477 memberships, 105,472 in HGNC, 105,464 with an Ensembl id, 105,307 with a UniProt accession.

### Positive control: acute versus never_infected (`acute/`)

**Cohorts and coverage** (`cohorts.csv`, `coverage.csv`):

| environment | role | positive / negative | features measured | Hallmark scored (of 50) | Reactome scored (of 1,867) | median fraction of members measured (scored programmes) |
|---|---|---|---|---|---|---|
| rna:GSE215865 | primary | 424 / 64 | 15173 | 50 | 1628 | 1.00 |
| protein:Cervia2024_SomaScan | primary | 113 / 34 | 6402 | 50 | 1540 | 1.00 |
| protein:Su2022 | primary | 204 / 289 | 441 | 33 | 177 | 0.11 |
| protein:Zoodsma2022_Hannover | primary | 12 / 62 | 359 | 26 | 122 | 0.12 |

**Concordance of programme effects** (`concordance_all.csv`; Spearman, with the number of programmes scored in both in parentheses; Hallmark and Reactome alone in `concordance_hallmark.csv` and `concordance_reactome.csv`):

| | rna:GSE215865 | protein:Cervia2024_SomaScan | protein:Su2022 | protein:Zoodsma2022_Hannover |
|---|---|---|---|---|
| rna:GSE215865 |  | -0.02 (1486) | 0.15 (210) | 0.13 (148) |
| protein:Cervia2024_SomaScan | -0.02 (1486) |  | 0.27 (210) | 0.31 (148) |
| protein:Su2022 | 0.15 (210) | 0.27 (210) |  | 0.37 (123) |
| protein:Zoodsma2022_Hannover | 0.13 (148) | 0.31 (148) | 0.37 (123) |  |

**Invariance partition** (`partition.csv`):

| analysis | environments | programmes partitioned (scored in every environment) | invariant | null | unstable | none | not_partitioned |
|---|---|---|---|---|---|---|---|
| primary | 4 | 123 | 0 | 0 | 123 | 0 | 1794 |

**Cross-layer leave one cohort out** (`loso.csv`, logreg; bootstrap 95% interval, 1000 resamples):

| scheme | test | role | positive / negative | programmes | AUROC identity (95% CI) | AUROC center_scale (95% CI) |
|---|---|---|---|---|---|---|
| protein_to_rna | rna:GSE215865 | primary | 424 / 64 | 123 | 0.637 (0.579 to 0.695) | 0.681 (0.626 to 0.733) |
| rna_plus_other_protein | protein:Cervia2024_SomaScan | primary | 113 / 34 | 123 | 0.472 (0.364 to 0.571) | 0.576 (0.471 to 0.683) |
| rna_plus_other_protein | protein:Su2022 | primary | 204 / 289 | 123 | 0.910 (0.880 to 0.937) | 0.829 (0.786 to 0.865) |
| rna_plus_other_protein | protein:Zoodsma2022_Hannover | primary | 12 / 62 | 123 | 0.996 (0.984 to 1.000) | 0.960 (0.909 to 0.993) |

### LC versus COVID_recovered (`lc/`)

**Cohorts and coverage:**

| environment | role | positive / negative | features measured | Hallmark scored (of 50) | Reactome scored (of 1,867) | median fraction of members measured (scored programmes) |
|---|---|---|---|---|---|---|
| rna:GSE226260_WB | primary | 17 / 20 | 12801 | 50 | 1617 | 0.89 |
| rna:GSE226260_PBMC | primary | 28 / 8 | 13866 | 50 | 1648 | 0.92 |
| rna:GSE224615 | primary | 23 / 13 | 15734 | 50 | 1660 | 0.93 |
| rna:GSE251849 | primary | 11 / 5 | 13872 | 50 | 1639 | 0.92 |
| rna:GSE270045 | primary | 19 / 17 | 15561 | 50 | 1739 | 0.97 |
| protein:Buggert2025 | primary | 146 / 51 | 2903 | 50 | 967 | 0.39 |
| protein:Cervia2024_SomaScan | primary | 40 / 73 | 6402 | 50 | 1540 | 0.91 |
| protein:Su2022 | primary | 75 / 50 | 441 | 33 | 177 | 0.09 |
| protein:Chowdhury2024 | primary | 21 / 9 | 404 | 22 | 134 | 0.14 |
| protein:Kovarik2023_plasma | primary | 13 / 13 | 314 | 19 | 144 | 0.12 |
| protein:Cervia2024_MS | sensitivity | 40 / 73 | 308 | 19 | 113 | 0.16 |

**Concordance of programme effects** (`concordance_all.csv`, as above; Cervia2024_MS is the sensitivity column, the same sera as Cervia2024_SomaScan):

| | rna:GSE226260_WB | rna:GSE226260_PBMC | rna:GSE224615 | rna:GSE251849 | rna:GSE270045 | protein:Buggert2025 | protein:Cervia2024_SomaScan | protein:Su2022 | protein:Chowdhury2024 | protein:Kovarik2023_plasma | protein:Cervia2024_MS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rna:GSE226260_WB |  | -0.06 (1662) | 0.25 (1665) | 0.21 (1657) | 0.19 (1661) | -0.04 (966) | 0.23 (1481) | 0.11 (210) | 0.38 (145) | -0.01 (153) | 0.35 (121) |
| rna:GSE226260_PBMC | -0.06 (1662) |  | -0.03 (1685) | -0.17 (1677) | -0.13 (1691) | -0.23 (976) | -0.08 (1501) | -0.13 (210) | 0.33 (148) | -0.26 (156) | 0.35 (124) |
| rna:GSE224615 | 0.25 (1665) | -0.03 (1685) |  | 0.24 (1681) | -0.08 (1701) | -0.07 (977) | 0.21 (1506) | -0.04 (210) | 0.50 (147) | -0.07 (156) | 0.60 (124) |
| rna:GSE251849 | 0.21 (1657) | -0.17 (1677) | 0.24 (1681) |  | -0.34 (1680) | -0.15 (974) | 0.16 (1496) | -0.07 (210) | -0.18 (149) | 0.04 (158) | -0.20 (126) |
| rna:GSE270045 | 0.19 (1661) | -0.13 (1691) | -0.08 (1701) | -0.34 (1680) |  | 0.07 (991) | 0.05 (1544) | 0.24 (210) | 0.30 (151) | 0.05 (159) | 0.31 (127) |
| protein:Buggert2025 | -0.04 (966) | -0.23 (976) | -0.07 (977) | -0.15 (974) | 0.07 (991) |  | -0.04 (1011) | -0.08 (210) | 0.16 (152) | -0.27 (160) | 0.17 (129) |
| protein:Cervia2024_SomaScan | 0.23 (1481) | -0.08 (1501) | 0.21 (1506) | 0.16 (1496) | 0.05 (1544) | -0.04 (1011) |  | 0.08 (210) | -0.01 (156) | 0.25 (163) | 0.21 (132) |
| protein:Su2022 | 0.11 (210) | -0.13 (210) | -0.04 (210) | -0.07 (210) | 0.24 (210) | -0.08 (210) | 0.08 (210) |  | 0.01 (86) | 0.03 (88) | -0.10 (73) |
| protein:Chowdhury2024 | 0.38 (145) | 0.33 (148) | 0.50 (147) | -0.18 (149) | 0.30 (151) | 0.16 (152) | -0.01 (156) | 0.01 (86) |  | -0.09 (139) | 0.72 (125) |
| protein:Kovarik2023_plasma | -0.01 (153) | -0.26 (156) | -0.07 (156) | 0.04 (158) | 0.05 (159) | -0.27 (160) | 0.25 (163) | 0.03 (88) | -0.09 (139) |  | -0.19 (124) |
| protein:Cervia2024_MS | 0.35 (121) | 0.35 (124) | 0.60 (124) | -0.20 (126) | 0.31 (127) | 0.17 (129) | 0.21 (132) | -0.10 (73) | 0.72 (125) | -0.19 (124) |  |

**Invariance partition** (`partition.csv`; `sensitivity_Cervia2024_MS` replaces Cervia2024_SomaScan with Cervia2024_MS):

| analysis | environments | programmes partitioned (scored in every environment) | invariant | null | unstable | none | not_partitioned |
|---|---|---|---|---|---|---|---|
| primary | 10 | 77 | 0 | 0 | 77 | 0 | 1840 |
| sensitivity_Cervia2024_MS | 10 | 67 | 0 | 0 | 67 | 0 | 1850 |

**Cross-layer leave one cohort out** (`loso.csv`). Cervia2024_MS is a test cohort only, under `rna_to_protein`:

| scheme | test | role | positive / negative | programmes | AUROC identity (95% CI) | AUROC center_scale (95% CI) |
|---|---|---|---|---|---|---|
| rna_to_protein | protein:Buggert2025 | primary | 146 / 51 | 961 | 0.568 (0.487 to 0.654) | 0.428 (0.348 to 0.513) |
| rna_to_protein | protein:Cervia2024_SomaScan | primary | 40 / 73 | 1470 | 0.536 (0.429 to 0.649) | 0.497 (0.382 to 0.608) |
| rna_to_protein | protein:Su2022 | primary | 75 / 50 | 210 | 0.659 (0.565 to 0.756) | 0.651 (0.561 to 0.747) |
| rna_to_protein | protein:Chowdhury2024 | primary | 21 / 9 | 144 | 0.614 (0.344 to 0.847) | 0.471 (0.238 to 0.704) |
| rna_to_protein | protein:Kovarik2023_plasma | primary | 13 / 13 | 153 | 0.562 (0.325 to 0.799) | 0.491 (0.260 to 0.728) |
| rna_to_protein | protein:Cervia2024_MS | sensitivity | 40 / 73 | 121 | 0.523 (0.410 to 0.638) | 0.468 (0.352 to 0.583) |
| protein_to_rna | rna:GSE226260_WB | primary | 17 / 20 | 77 | 0.682 (0.476 to 0.844) | 0.638 (0.453 to 0.806) |
| protein_to_rna | rna:GSE226260_PBMC | primary | 28 / 8 | 77 | 0.518 (0.295 to 0.746) | 0.473 (0.254 to 0.705) |
| protein_to_rna | rna:GSE224615 | primary | 23 / 13 | 77 | 0.548 (0.334 to 0.763) | 0.515 (0.294 to 0.739) |
| protein_to_rna | rna:GSE251849 | primary | 11 / 5 | 77 | 0.364 (0.109 to 0.655) | 0.309 (0.054 to 0.600) |
| protein_to_rna | rna:GSE270045 | primary | 19 / 17 | 77 | 0.851 (0.703 to 0.963) | 0.845 (0.709 to 0.954) |
| leave_one_cohort_out | rna:GSE226260_WB | primary | 17 / 20 | 77 | 0.688 (0.494 to 0.841) | 0.647 (0.459 to 0.812) |
| leave_one_cohort_out | rna:GSE226260_PBMC | primary | 28 / 8 | 77 | 0.464 (0.268 to 0.679) | 0.469 (0.254 to 0.705) |
| leave_one_cohort_out | rna:GSE224615 | primary | 23 / 13 | 77 | 0.552 (0.324 to 0.763) | 0.532 (0.311 to 0.746) |
| leave_one_cohort_out | rna:GSE251849 | primary | 11 / 5 | 77 | 0.400 (0.091 to 0.710) | 0.345 (0.073 to 0.655) |
| leave_one_cohort_out | rna:GSE270045 | primary | 19 / 17 | 77 | 0.789 (0.625 to 0.926) | 0.715 (0.542 to 0.870) |
| leave_one_cohort_out | protein:Buggert2025 | primary | 146 / 51 | 77 | 0.578 (0.485 to 0.667) | 0.551 (0.456 to 0.631) |
| leave_one_cohort_out | protein:Cervia2024_SomaScan | primary | 40 / 73 | 77 | 0.537 (0.427 to 0.642) | 0.554 (0.447 to 0.661) |
| leave_one_cohort_out | protein:Su2022 | primary | 75 / 50 | 77 | 0.512 (0.407 to 0.618) | 0.478 (0.377 to 0.583) |
| leave_one_cohort_out | protein:Chowdhury2024 | primary | 21 / 9 | 77 | 0.545 (0.338 to 0.730) | 0.529 (0.323 to 0.730) |
| leave_one_cohort_out | protein:Kovarik2023_plasma | primary | 13 / 13 | 77 | 0.704 (0.485 to 0.893) | 0.781 (0.568 to 0.953) |

**Bridge check in Cervia2024** (`bridge_cervia.csv`). SomaScan and MS programme scores of the same 113 sera, Spearman correlation per programme across the sera; `own`: each layer's own measured members, `shared`: only members measured in both layers.

| member set | collection | programmes scored in both layers | median Spearman | minimum | maximum | programmes with Spearman p < 0.05 |
|---|---|---|---|---|---|---|
| own | hallmark | 19 | 0.09 | -0.15 | 0.63 | 4 |
| own | reactome | 113 | 0.14 | -0.15 | 0.73 | 47 |
| shared | hallmark | 16 | 0.34 | 0.05 | 0.58 | 13 |
| shared | reactome | 84 | 0.34 | -0.02 | 0.71 | 62 |

**Column glossary.**
- `cohorts.csv`: `env`, `layer`, `cohort`, `role` (primary or sensitivity), `n_positive`, `n_negative`, `n_features_measured`, `n_scored_<collection>`, `n_programmes_<collection>`.
- `coverage.csv`: `layer`, `env`, `programme`, `collection`, `n_members` (members in the layer's feature space), `n_measured`, `fraction_measured`, `scored`.
- `effects.csv`, `effects_var.csv`: per programme, `collection` and Hedges g (or its variance) per environment; missing where the programme is not scored.
- `partition.csv`: `analysis`, `programme`, `collection`, the meta-analysis columns of `chorus.features.meta_analyse` (`effect`, `se`, `z`, `p`, `Q`, `Q_p`, `I2`, `tau2`, `k`, `sign_consistency`), `set`, `max_abs_g`, `n_env_scored`, `n_env`, `environments`; the meta-analysis columns are empty for `not_partitioned`.
- `loso.csv`: `scheme`, `test_role`, `train`, `test`, `n_train`, `n_test_pos`, `n_test_neg`, `n_programmes`, `harmoniser`, `classifier`, the metrics of loso_results.csv, `auroc_lo`, `auroc_hi`, `skipped` (empty: nothing was skipped).
- `loso_predictions.csv`: `scheme`, `test`, `harmoniser`, `classifier`, `sample` (`<cohort>::<sample_id>`), `y`, `p`.
- `bridge_cervia.csv`: `member_set`, `programme`, `collection`, `n_samples`, `spearman`, `spearman_p`, `pearson`, `pearson_p`, `n_members_somascan`, `n_members_ms`, `n_members_shared`.

**Code and config changes:**
- `chorus.integrate.programmes`, new.
- `preprocessing: features: per_cohort`: every cohort keeps the features that pass its own filter.
- `acute_lc_status`: per cohort, for `relabel_acute`.
- New configs: `configs/mosaic_lc_rna.yaml`, `configs/mosaic_acute_rna.yaml` and `configs/mosaic_acute_protein.yaml`.

**Not committed:** `scores.csv.gz` and the log (`results/mosaic_programmes.log`).

## Metabolomics (16 November 2026)

**Cohorts and identifiers** (`docs/metabolomics_cohorts.md`).
- **Built:** Martinez2024 (Metabolomics Workbench ST003103), Lopez2023 (Mendeley 8zfdjsypd8), Saito2024
  (Front Immunol supplementary Table 2) and Holmes2021 (Zenodo 4713391).
- **Not built:** Mikuteit2023: its supplement holds group summaries only.
- **Labels:** Holmes2021 is built but `lc_status` unknown, because the deposit does not explain its symptom composites.
- **Identifiers:** metabolite names are mapped to a shared identifier through RefMet: `HMDB:<id>`, else
  `RefMet:<name>` (`chorus.metabolites`, adapter `metabolomics_id`).
- **Mapping tables:** `lc_metabolomics/mapping/`:
  - `mapping_summary.csv`: per cohort, names mapped by each rule;
  - `n_shared.csv`: identifiers shared by every pair, Su2022_metabolomics included;
  - `shared_by_k.csv`: identifiers in at least k cohorts, and shared by each triple of labelled cohorts.

**Why these three cohorts for Test 1.** Of the labelled cohorts (Martinez2024, Lopez2023, Saito2024,
Su2022_metabolomics), one triple shares 30 or more identifiers: Martinez2024, Lopez2023 and Su2022_metabolomics (38).
- **Other triples:** Martinez2024 + Lopez2023 + Saito2024 share 15; Martinez2024 + Saito2024 + Su2022_metabolomics
  share 17; Lopez2023 + Saito2024 + Su2022_metabolomics share 25.
- **All four:** share 15.
- **So:** Test 1 runs on that triple and its 38 identifiers (`configs/lc_metabolomics.yaml`). Saito2024 enters the
  concordance only.

**Test 1 setup** (`lc_metabolomics`). Leave one study out, as for proteomics:
- **Harmonisers:** identity, center_scale, rank, fsqn, combat_addon and chorus, plus nested.
- **Selectors:** variance, differential and invariant. Their sizes are the proteomics fractions of the shared features
  (100 and 50 of 376) applied to 38: variance 10, differential and invariant 5.
- **Classifiers:** logreg and rf; nested selection on.
- **Within-study CV:** 5 folds × 10 repeats (per-fold method).
- **Samples:** post-acute, one per subject: Martinez2024 40 LC / 40 COVID_recovered, Lopez2023 30 / 18,
  Su2022_metabolomics 68 / 45 (T3).
- **Imputation:** `study_min`. 29 values were filled in Martinez2024 (one sample has no GC-MS values) and 1 in
  Lopez2023.
- **Su2022_metabolomics** holds the authors' adjusted Z scores, so the harmoniser table of "Su2022 and the
  harmonisers" (Proteomics, 2.) applies to it in the same way. The other two cohorts are log2 areas and log2 µM.

**chorus.** All 9 chorus rows (3 held-out studies × 3 selectors) carry an error: only 12 or 14 null metabolites
passed the screen, against the 20 required. As for the acute proteomics control, the thresholds were not relaxed.
The harmonisation diagnostics carry the same message for chorus. `nested` chooses among the other five candidates
when chorus fails in an inner fold. Choices over the 18 nested rows: rank 8, center_scale 4, combat_addon 2, identity 2, fsqn 2.

**LOSO, mean AUROC across the three held-out studies (minimum to maximum)** (`loso_summary.csv`; per held-out
rows with bootstrap intervals in `loso_results.csv`):

| harmoniser | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| identity | 0.594 (0.504 to 0.696) | 0.536 (0.425 to 0.626) | 0.586 (0.492 to 0.7) | 0.519 (0.5 to 0.534) | 0.528 (0.489 to 0.58) | 0.439 (0.369 to 0.547) |
| center_scale | 0.621 (0.577 to 0.667) | 0.596 (0.502 to 0.649) | 0.607 (0.575 to 0.644) | 0.501 (0.428 to 0.551) | 0.569 (0.538 to 0.593) | 0.586 (0.495 to 0.71) |
| rank | 0.554 (0.528 to 0.574) | 0.553 (0.479 to 0.639) | 0.539 (0.476 to 0.582) | 0.486 (0.479 to 0.491) | 0.567 (0.485 to 0.622) | 0.558 (0.444 to 0.71) |
| fsqn | 0.569 (0.498 to 0.622) | 0.544 (0.519 to 0.57) | 0.547 (0.51 to 0.594) | 0.475 (0.41 to 0.514) | 0.563 (0.469 to 0.668) | 0.505 (0.445 to 0.56) |
| combat_addon | 0.619 (0.581 to 0.666) | 0.593 (0.475 to 0.653) | 0.607 (0.575 to 0.643) | 0.501 (0.419 to 0.564) | 0.557 (0.455 to 0.721) | 0.519 (0.473 to 0.6) |
| nested | 0.606 (0.574 to 0.666) | 0.598 (0.502 to 0.652) | 0.574 (0.566 to 0.582) | 0.48 (0.428 to 0.534) | 0.525 (0.469 to 0.622) | 0.533 (0.444 to 0.71) |

**Within-study CV, mean AUROC over folds** (`within_study_cv.csv`; per fold in `within_study_cv_folds.csv`). With
k = 5, differential and invariant select the same metabolites within one study, so their rows are equal.

| study | differential logreg | differential rf | invariant logreg | invariant rf | variance logreg | variance rf |
|---|---|---|---|---|---|---|
| Lopez2023 | 0.661 | 0.559 | 0.661 | 0.559 | 0.621 | 0.59 |
| Martinez2024 | 0.804 | 0.765 | 0.804 | 0.765 | 0.779 | 0.821 |
| Su2022_metabolomics | 0.511 | 0.508 | 0.511 | 0.508 | 0.519 | 0.513 |

**Harmonisation diagnostics** (`lc_metabolomics/harmonisation`). These are the transcriptomics metrics, run with
`--n-pcs 10` because there are 38 features. `all_metrics.csv` has 213 rows. Chorus has one row per held-out study, with family `error` and the null screen message in `scope`.

**Concordance** (`lc_metabolomics/concordance`: the three Test 1 cohorts on their 38 identifiers, with `top20` files;
`lc_metabolomics/concordance_all`: the four labelled cohorts each on its own identifiers, pairwise over those shared,
`n_shared.csv`, with `top50` files). Spearman correlation of Hedges g (LC minus COVID_recovered):

Test 1 cohorts, 38 identifiers:

| | Martinez2024 | Lopez2023 | Su2022_metabolomics |
|---|---|---|---|
| Martinez2024 |  | 0.112 (38) | 0.532 (38) |
| Lopez2023 | 0.112 (38) |  | 0.345 (38) |
| Su2022_metabolomics | 0.532 (38) | 0.345 (38) |  |

All labelled cohorts, pairwise over shared identifiers (number shared in parentheses):

| | Martinez2024 | Lopez2023 | Saito2024 | Su2022_metabolomics |
|---|---|---|---|---|
| Martinez2024 |  | 0.156 (41) | -0.287 (17) | 0.058 (158) |
| Lopez2023 | 0.156 (41) |  | 0.007 (28) | 0.213 (61) |
| Saito2024 | -0.287 (17) | 0.007 (28) |  | 0.104 (66) |
| Su2022_metabolomics | 0.058 (158) | 0.213 (61) | 0.104 (66) |  |

- **`top50` files of `concordance_all`:** the 50 identifiers with the largest mean |g| are taken over the union of
  the four cohorts. So each pair of cohorts has 1 to 4 of them in common, and the top-50 correlations rest on those
  1 to 4 values.

**Epigenomics** (`concordance_epigenomics/cohorts.csv`, `docs/metabolomics_cohorts.md`). GSE174818, GSE197152 and
GSE247869 were downloaded (EPIC, whole blood); the groups are documented, but no per-CpG effect was computed.
- **GSE174818:** acute samples only.
- **GSE197152:** 15 one-year samples, all described with PASC symptoms, so there is no recovered group.
- **GSE247869:** 96 six-month samples. The paper's PASC subgroup (28) is not identifiable per sample in GEO or in
  the supplement.
- **So:** no series has per-sample LC and COVID_recovered labels.

**Code and config changes:**
- `chorus.metabolites` (RefMet matching with a cache, identifier rules, pairwise shared counts);
- the `metabolomics_id` adapter;
- `scripts/metabolomics/build_cohorts.py` and `scripts/metabolomics/map_ids.py`;
- `scripts/concordance.py`, which labels `metabolite_id` features with the RefMet name;
- `configs/lc_metabolomics.yaml` and `configs/concordance_metabolomics_all.yaml`;
- xlrd 2.0.2 added to the venv.

**Not committed:** the logs and `selected_features.json`.
