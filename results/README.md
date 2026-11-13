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
changes".

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
| `control_sex_label` | `scripts/control_runs.py sex configs/lc_transcriptomics_blood.yaml --out results/control_sex_label` | configs/lc_transcriptomics_blood.yaml, changed by the script | 113 (version 3, started 2026-11-13 15:28:14) |
| `lc_transcriptomics_blood_harmonisation_v1` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` (folder renamed) | configs/lc_transcriptomics_blood.yaml | 100 (started 2026-11-13 09:01:09) |
| `lc_transcriptomics_blood_harmonisation` | `chorus harmonisation configs/lc_transcriptomics_blood.yaml --out results/lc_transcriptomics_blood_harmonisation` | configs/lc_transcriptomics_blood.yaml | 104 (version 3, started 2026-11-13 15:26:30) |
| `concordance` | `scripts/concordance.py configs/lc_transcriptomics_blood.yaml --out results/concordance` | configs/lc_transcriptomics_blood.yaml | 48 (started 2026-11-13 09:03:45) |
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
    rows of `forced_sex_genes` therefore equal those of `all_genes`.
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
  - **Status.** `control_sex_label` was **not rerun** (no benchmark runs were requested), so its `forced_sex_genes`
    rows still come from the version 3 code.
- **Layer adapters (commit b5083c4).** The transcriptomics loader now goes through the `bulk_rna` adapter. The loaded
  matrices of `configs/lc_transcriptomics_blood.yaml` and `configs/lc_acute_transfer.yaml` were checked to be identical
  (hash of every cohort matrix) before and after the change.
