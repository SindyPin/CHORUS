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
made one after another. On 13 November `lc_transcriptomics_blood_within_v2` ran at the same time as the two control runs,
so the three runtimes overlap.

| folder | command | config | runtime (s) |
|---|---|---|---|
| `lc_transcriptomics_blood` | `chorus run configs/lc_transcriptomics_blood.yaml --within` | configs/lc_transcriptomics_blood.yaml | 1053 (started 2026-11-12 09:19:36) |
| `lc_transcriptomics_blood_noGSE270045` | `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --within` | configs/lc_transcriptomics_blood_noGSE270045.yaml | 805 (started 2026-11-12 09:37:09) |
| `lc_acute_transfer` | `chorus run configs/lc_acute_transfer.yaml --within` | configs/lc_acute_transfer.yaml | 231 (started 2026-11-12 09:50:34) |
| `simulation_factors` | `scripts/simulation_benchmark.py --out results/simulation_factors` | built into the script | 231 (started 2026-11-11 09:09:10) |
| `simulation_shifts_only` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` | built into the script | 346 (started 2026-11-11 09:13:14) |
| `lc_transcriptomics_blood_within_v2` | `chorus run configs/lc_transcriptomics_blood.yaml --within-only --out results/lc_transcriptomics_blood_within_v2` | configs/lc_transcriptomics_blood.yaml | 1320 (started 2026-11-13 07:05:33) |
| `control_matrix_label` | `scripts/control_runs.py matrix configs/lc_transcriptomics_blood.yaml --out results/control_matrix_label` | configs/lc_transcriptomics_blood.yaml, changed by the script | 49 (started 2026-11-13 07:06:20) |
| `control_label_permutation` | `scripts/control_runs.py permute configs/lc_transcriptomics_blood.yaml --out results/control_label_permutation --n-perm 20` | configs/lc_transcriptomics_blood.yaml, changed by the script | 192 (started 2026-11-13 07:07:17) |

- `lc_transcriptomics_blood`: primary leave-one-study-out benchmark, full grid (8 harmonisers plus nested selection,
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
- `control_matrix_label`: positive control. The label is `sample_matrix` (whole_blood positive, PBMC negative). See
  "Control runs" below.
- `control_label_permutation`: label permutation control, 20 permutations. See "Control runs" below.

The cohort table each real-data run loaded (samples and labels per study) is printed at the top of its log,
`results/<folder>.log`. The logs stay local (not committed). Simulation runs and `control_matrix_label` write it to
`cohorts.csv`.

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
| control_matrix_label | GSE226260 (WB and PBMC merged) | GSE224615, GSE251849, GSE270045 |
| control_label_permutation | as lc_transcriptomics_blood | none |

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
- `cohorts.csv`: simulation runs and `control_matrix_label`. `figure_simulation.png`: simulation runs only.

## Column glossary

**loso_results.csv**

| column | meaning |
|---|---|
| held_out | study left out and used as the test set |
| harmoniser | harmoniser in the grid; `nested` means it was chosen by an inner leave-one-study-out over the training studies |
| harmoniser_chosen | harmoniser actually applied (differs from `harmoniser` only for `nested`) |
| leak_risk | true when the harmoniser is transductive (fitted on all studies, without labels) rather than on training studies only |
| selector | selector configuration: `variance`, `differential`, `invariant`, or `prior:<gene list>:<mode>[:k<k>]` |
| n_features | number of model inputs: selected genes plus covariate columns |
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

**control_matrix_label (positive control).**
- **Label:** `lc_status` is replaced by `sample_matrix`, with whole_blood positive and PBMC negative.
- **Held-out study:** GSE226260_WB and GSE226260_PBMC are merged into one held-out study, GSE226260, because a single-matrix
  cohort would give a one-class test set. The other cohorts (GSE224615, GSE251849, GSE270045) are training only.
- **Grid:** harmonisers identity, center_scale, combat_addon and chorus; selectors variance (k 500) and invariant (k 200);
  classifier logreg; no nested selection; 1000 bootstrap resamples.
- **Covariates:** `sample_matrix` is dropped from the covariates and from `protect`, because it is the label.
- **Training studies hold one class each.** GSE224615 and GSE251849 are PBMC only, and GSE270045 is whole blood only.
  As a result, in `loso_results.csv`:
  - 5 of the 8 grid cells failed, and their `error` column says why.
  - Every invariant cell failed with "no study had enough cases and controls for effect estimation", because the invariant
    screen needs both classes within a training study.
  - chorus failed with "invariance anchored harmonisation needs at least two labelled training studies".
  - center_scale centres and scales each training study separately, which removes between-study differences, here the
    label difference. All its predictions are equal (AUROC 0.500, interval 0.50 to 0.50).
- **The merged held-out study differs in more than matrix.** In GSE226260 the matrix goes together with cohort,
  enrolment period and library chemistry (`docs/gse226260_check.md`).

**control_label_permutation.**
- **Permutation:** within each cohort, `lc_status` is permuted among the samples labelled LC or COVID_recovered, so class
  counts per cohort are unchanged.
- **Benchmark:** the primary leave-one-study-out benchmark with harmonisers identity and chorus, the invariant selector
  (k 200) and logreg. Covariates and `protect` are `sample_matrix` as in the primary config. There is no nested selection.
- **Repeats:** 20 permutations (seed 0), 1000 bootstrap resamples per fit (the interval columns are not summarised).
