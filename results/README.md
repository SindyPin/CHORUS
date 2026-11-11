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

Runtimes are wall-clock seconds from `results/runtimes.tsv`. The runs were made one after another on the same laptop.

| folder | command | config | runtime (s) |
|---|---|---|---|
| `lc_transcriptomics_blood` | `chorus run configs/lc_transcriptomics_blood.yaml --within` | configs/lc_transcriptomics_blood.yaml | 1053 (started 2026-11-12 09:19:36) |
| `lc_transcriptomics_blood_noGSE270045` | `chorus run configs/lc_transcriptomics_blood_noGSE270045.yaml --within` | configs/lc_transcriptomics_blood_noGSE270045.yaml | 805 (started 2026-11-12 09:37:09) |
| `lc_acute_transfer` | `chorus run configs/lc_acute_transfer.yaml --within` | configs/lc_acute_transfer.yaml | 231 (started 2026-11-12 09:50:34) |
| `simulation_factors` | `scripts/simulation_benchmark.py --out results/simulation_factors` | built into the script | 231 (started 2026-11-11 09:09:10) |
| `simulation_shifts_only` | `scripts/simulation_benchmark.py --out results/simulation_shifts_only --n-factors 0` | built into the script | 346 (started 2026-11-11 09:13:14) |

- `lc_transcriptomics_blood`: primary leave-one-study-out benchmark, full grid (8 harmonisers plus nested selection,
  6 selector configurations, 4 classifiers), with the within-study cross-validation reference.
- `lc_transcriptomics_blood_noGSE270045`: identical except GSE270045 is disabled (sensitivity run;
  `docs/gse270045_check.md`).
- `lc_acute_transfer`: post-acute cohorts as training only, GSE215865 (acute admissions) as the single held-out study
  (`docs/acute_transfer.md`).
- `simulation_factors`: `scripts/simulation_benchmark.py` with its defaults (6 simulated studies, 3 sample-level
  technical factors).
- `simulation_shifts_only`: the same with `--n-factors 0` (study location and scale shifts only).

The cohort table each real-data run loaded (samples and labels per study) is printed at the top of its log,
`results/<folder>.log`. The logs stay local (not committed). Simulation runs write it to `cohorts.csv`.

`smoke_transcriptomics` and `smoke_acute_transfer` are smoke tests from 12 November 2026 on the same cohorts with a reduced
grid: harmonisers identity and center_scale, selectors variance, differential and invariant, classifier logreg, no nested
selection, 200 bootstrap resamples.

## Cohorts in the real-data runs

| run | held out | training only |
|---|---|---|
| lc_transcriptomics_blood | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 | none |
| lc_transcriptomics_blood_noGSE270045 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849 | none |
| lc_acute_transfer | GSE215865 | GSE226260_WB, GSE226260_PBMC, GSE224615, GSE251849, GSE270045 |

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
- `cohorts.csv`, `figure_simulation.png`: simulation runs only.

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

**within_study_cv.csv**: `study`, `selector`, `classifier` and the same metrics as above, from 5-fold stratified
cross-validation inside the study. The invariant selector falls back to differential, because one study gives one
environment. `cross_study` is skipped. Studies with fewer than 5 samples in a class are skipped. Training-only studies
are included.
